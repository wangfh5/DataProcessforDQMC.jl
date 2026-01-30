using Statistics
using DataFrames
using Printf

export error, round_error, format_value_error
export iqr_fence_filter, outlier_filter, remove_outliers
export compute_cr_m2_covariance

## -------------------------------------------------------------------------- ##
##                              Outlier filtering                             ##
## -------------------------------------------------------------------------- ##

"""
    iqr_fence_filter(values; k=10.0, min_n=5)

Apply Tukey's IQR fence to remove extreme outliers.

Given samples `x₁,…,xₙ`, define the quartiles and fences; keep points in-range.

# Behaviour / edge cases
- Ignore `missing` / non-finite before quantiles.
- If finite samples `< min_n`, or `IQR <= 0`, keep all finite.
- If fence would remove all, fall back to keeping all finite.

# Returns
`(; keep::BitVector, filtered::Vector{Float64}, removed_min::Int, removed_max::Int)`
"""
function iqr_fence_filter(values; k::Real=10.0, min_n::Int=5)
    keep = falses(length(values))
    finite_idx = findall(x -> !ismissing(x) && isfinite(x), values)
    n_finite = length(finite_idx)
    removed_min = 0
    removed_max = 0

    if n_finite < min_n
        keep[finite_idx] .= true
        return (; keep, filtered=values[keep], removed_min, removed_max)
    end

    vals = Float64.(values[finite_idx])
    q1 = quantile(vals, 0.25)
    q3 = quantile(vals, 0.75)
    iqr = q3 - q1
    if !(isfinite(iqr) && iqr > 0)
        keep[finite_idx] .= true
        return (; keep, filtered=values[keep], removed_min, removed_max)
    end

    kk = Float64(k)
    kk <= 0 && throw(ArgumentError("iqrfence k must be > 0, got $kk"))
    lo = q1 - kk * iqr
    hi = q3 + kk * iqr

    inrange = (lo .<= vals) .& (vals .<= hi)
    if !any(inrange)
        keep[finite_idx] .= true
        removed_min = 0
        removed_max = 0
        return (; keep, filtered=values[keep], removed_min, removed_max)
    end

    keep[finite_idx[inrange]] .= true
    removed_min = count(<(lo), vals)
    removed_max = count(>(hi), vals)

    return (; keep, filtered=values[keep], removed_min, removed_max)
end

"""
    outlier_filter(values, mode, param; min_n=5)

统一的离群值处理入口（公开）。返回 keep mask 与筛选后的值，并附带删除计数。
"""
function outlier_filter(values::AbstractVector, mode, param; min_n::Int=5)
    m = Symbol(lowercase(String(mode)))
    dropnum = 0
    if m == :dropmaxmin
        dropnum = Int(param)
        dropnum < 0 && throw(ArgumentError("dropmaxmin must be >= 0, got $dropnum"))
    elseif m == :iqrfence
        k = Float64(param)
        k <= 0 && throw(ArgumentError("iqrfence k must be > 0, got $k"))
    else
        throw(ArgumentError("mode must be :dropmaxmin or :iqrfence, got $mode"))
    end

    keep = falses(length(values))
    finite_idx = findall(x -> !ismissing(x) && isfinite(x), values)
    n_finite = length(finite_idx)

    if n_finite < min_n
        keep[finite_idx] .= true
        return (; keep, filtered=values[keep], removed_min=0, removed_max=0)
    end

    if m == :dropmaxmin
        if n_finite <= 2 * dropnum
            keep[finite_idx] .= true
            return (; keep, filtered=values[keep], removed_min=0, removed_max=0)
        end
        vals = Float64.(values[finite_idx])
        order = sortperm(vals)
        kept = order[(dropnum + 1):(end - dropnum)]
        keep[finite_idx[kept]] .= true
        return (; keep, filtered=values[keep], removed_min=dropnum, removed_max=dropnum)
    elseif m == :iqrfence
        return iqr_fence_filter(values; k=param, min_n)
    end
end

"""
    remove_outliers(values, mode, param; min_n=5)

统一的离群值处理入口（返回筛选后数据 + 删除计数）。
"""
function remove_outliers(values::AbstractVector, mode, param; min_n::Int=5)
    res = outlier_filter(values, mode, param; min_n)
    return (; values=res.filtered, removed_min=res.removed_min, removed_max=res.removed_max)
end

## -------------------------------------------------------------------------- ##
##                           Statistical Quantities                           ##
## -------------------------------------------------------------------------- ##

"""
    error(data; sigma=1, bessel=true, auto_digits=false)
Calculate the standard error of the data.

Arguments:
- `data`: The data array
- `sigma`: Number of standard deviations (default: 1)
- `bessel`: Whether to use Bessel's correction (N-1) for sample standard deviation (default: true)
- `auto_digits`: Whether to automatically determine significant digits based on error of error (default: true)

Returns:
- Standard error multiplied by sigma, with appropriate precision if auto_digits=true
"""
function error(data; sigma=1, bessel=true, auto_digits=true)
    # Handle empty or single-element arrays
    if length(data) <= 1
        return zero(eltype(data))
    end

    n = length(data)

    # Calculate standard deviation with or without Bessel's correction
    divisor = bessel ? n - 1 : n
    mean_val = mean(data)
    s = sqrt(sum((data .- mean_val).^2) / divisor)

    # Calculate standard error
    std_err = s / sqrt(n)

    # Apply sigma multiplier
    result = sigma * std_err

    # Adjust significant digits based on error of error if requested
    if auto_digits
        # Error of error formula: std_err / sqrt(2*(n-1))
        err_of_err = std_err / sqrt(2 * (n - 1))

        # Use the round_error function to round the result
        result = round_error(result, err_of_err)[1]
    end

    return result
end


## -------------------------------------------------------------------------- ##
##                                Format error                                ##
## -------------------------------------------------------------------------- ##

"""
    round_error(error_val, error_of_error)
Round an error value based on its error of error.

The error is rounded to be precise to one digit before the first significant digit of error of error.
This ensures that the reported error has appropriate precision.

Arguments:
- `error_val`: The error value to be rounded
- `error_of_error`: The error of the error (uncertainty in the error estimate)

Returns:
- Rounded error value with appropriate precision
- The number of significant digits in the rounded error

Examples:
```julia
round_error(23456, 2345)    # Returns (30000,1)
round_error(2.3456, 0.002345)  # Returns (2.35,3)
```
"""
function round_error(error_val::AbstractFloat, error_of_error::Number)
    # Handle zero error of error case
    if error_of_error == 0
        return error_val, -1
    end

    # Calculate the number of significant digits in the error
    # For example:
    # - If error = 23456, error_of_error = 2345, significant_digits = 4 - 3 = 1
    # - If error = 2.3456, error_of_error = 0.002345, significant_digits = 0 - (-3) = 3
    err_of_err_order = floor(log10(abs(error_of_error)))
    error_order = floor(log10(abs(error_val)))
    significant_digits = Int(error_order - err_of_err_order)

    # Round the error to the appropriate precision using RoundUp mode
    # This ensures we never underestimate the error
    result = round(error_val, RoundUp; sigdigits=significant_digits)

    return result, significant_digits
end
function round_error(error_val::Int, error_of_error::Number)
    result, significant_digits = round_error(Float64(error_val), error_of_error)

    return Int(result), significant_digits
end
function round_error(error_val::Complex{T}, error_of_error::Complex{T}) where T <: AbstractFloat
    real_rounded, real_sig_digits = round_error(real(error_val), real(error_of_error))
    imag_rounded, imag_sig_digits = round_error(imag(error_val), imag(error_of_error))

    return complex(real_rounded, imag_rounded), (real_sig_digits, imag_sig_digits)
end

"""
    format_value_error(value, error, error_sig_digits=1; format=:scientific)
格式化数值及其误差，支持科学计数法或普通小数表示。

误差会按照给定的有效数字数目进行向上取整，数值会匹配同样的精度。

参数:
- `value`: 待格式化的数值
- `error`: 该数值的误差
- `error_sig_digits`: 误差保留的有效数字个数 (默认: 1)
- `format`: 输出格式，`:scientific`(默认) 返回科学计数法字符串，`:decimal` 返回普通小数

返回:
- `(formatted_value, formatted_error)` 字符串元组，格式由 `format` 决定

示例:
```julia
# 科学计数法 (默认)
format_value_error(2.36738, 0.0023)                      # ("2.367e+00", "0.003e0")
format_value_error(2367.38, 23, 2)                       # ("2.367e+03", "0.023e3")

# 普通小数
format_value_error(2.36738, 0.0023; format=:decimal)     # ("2.367", "0.003")
format_value_error(2367.38, 23; format=:decimal)         # ("2370", "30")
```
"""
function format_value_error(value::Number, error::Number, error_sig_digits::Int=1; format::Symbol=:scientific)
    @assert format in [:scientific, :decimal] "format must be :scientific or :decimal, got $format"
    # Step 0: Special-case zero error.
    if iszero(error)
        if format == :scientific
            x = Float64(value)
            raw = @sprintf("%.15e", x)
            mantissa, exp_str = split(raw, 'e')
            exponent = parse(Int, exp_str)
            if occursin('.', mantissa)
                mantissa = replace(mantissa, r"0+$" => "")
                mantissa = replace(mantissa, r"\.$" => "")
            end
            val_str = mantissa * "e" * exp_str
            err_str = "0e$(exponent)"
            return val_str, err_str
        elseif format == :decimal
            return string(value), "0"
        end
    end

    # Step 1: Determine the digits of the error based on its significant digits and error
    # Round the error first, then derive its order; this keeps value precision in sync
    rounded_error = round(error, RoundUp, sigdigits=error_sig_digits)
    # For non-zero errors, calculate order of magnitude (err = x.x × 10^error_order)
    # For example
    # - If error = 23456, error_order = 4;
    # - If error = 2.3456, error_order = 0;
    # - If error = 0.00943, error_sig_digits = 2, rounded_error = 0.0095, error_order = -3;
    # - If error = 0.00943, error_sig_digits = 1, rounded_error = 0.01, error_order = -2;
    error_order = floor(log10(abs(rounded_error)))
    # Check if error_order is -Inf (can happen with very small numbers due to floating point precision)
    if isinf(error_order)
        # Use the smallest representable order for very small numbers
        error_order = -324  # Approximately the smallest exponent for Float64
    end
    # For example
    # - If error = 23456, error_sig_digits = 1, error_digits = - 4;
    # - If error = 23456, error_sig_digits = 2, error_digits = - 4 + 2 - 1 = - 3;
    # - If error = 2.3456, error_sig_digits = 1, error_digits = 0;
    # - If error = 2.3456, error_sig_digits = 3, error_digits = 0 + 3 - 1 = 2;
    # - If error = 0.00943, error_sig_digits = 2, error_digits = 3 + 2 - 1 = 4;
    # - If error = 0.00943, error_sig_digits = 1, error_digits = 2 + 1 - 1 = 2;
    error_digits = - Int(error_order - error_sig_digits + 1)

    # Step 2: Round the value to match the precision of the error
    value_digits = error_digits
    rounded_val = round(value, digits=value_digits)

    # Step 3: Determine the number of significant digits in the value
    # For example
    # - If value = 2.367, value_digits = 3, value_order = 0, value_sig_digits = 0 + 1 + 3 = 4;
    # - If value = 2367, value_digits = 0, value_order = 3, value_sig_digits = 3 + 1 + 0 = 4;
    # - If value = 2300, value_digits = -2, value_order = 3, value_sig_digits = - 2 + 1 + 3 = 2;
    # - If value = 0, handle as a special case
    if rounded_val == 0.0  # Exactly zero
        value_order = 0
        value_sig_digits = 1 + value_digits
    else
        value_order = floor(log10(abs(rounded_val)))
        # Check if value_order is -Inf (can happen with very small numbers due to floating point precision)
        if isinf(value_order)
            # Use the smallest representable order for very small numbers
            value_order = -324  # Approximately the smallest exponent for Float64
        end
        value_sig_digits = Int(value_order) + 1 + value_digits
    end

    # Step 4: Format the value and error as strings
    err_exponent = Int(value_order)  # Safe now because we've handled the zero case above

    if format == :scientific
        # Format as scientific notation
        # For clarity, the exponents are printed to be the same for both values
        # Value is automatically formatted, while error is formatted manually
        # For example
        # - If value = 2.36738, error = 0.0023, error_sig_digits = 1, then we want 2.367e+00 and 0.003e+00
        # - If value = 2.36738, error = 0.0023, error_sig_digits = 2, then we want 2.3674e+00 and 0.0023e+00
        # - If value = 2367.38, error = 23, error_sig_digits = 1, then we want 2.37e+03 and 0.03e+03
        # - If value = 2367.38, error = 23, error_sig_digits = 2, then we want 2.367e+03 and 0.023e+03
        val_str = @sprintf("%.*e", value_sig_digits - 1, rounded_val)

        if rounded_error == 0.0
            err_str = "0e$(err_exponent)"
        else
            err_magnitude = abs(rounded_error) / 10.0^err_exponent
            rounded_error_format = round(err_magnitude, sigdigits=error_sig_digits)
            err_str = "$(rounded_error_format)e$(err_exponent)"
        end
    elseif format == :decimal
        # Format as plain decimal
        # Both value and error are rounded to match precision, then printed with fixed decimal places
        # For example
        # - If value = 2.36738, error = 0.0023, error_sig_digits = 1, then we want "2.367" and "0.003"
        # - If value = 2367.38, error = 23, error_sig_digits = 1, then we want "2370" and "30"
        value_decimal_digits = max(value_digits, 0)
        error_decimal_digits = max(error_digits, 0)

        val_str = @sprintf("%.*f", value_decimal_digits, rounded_val)
        err_str = @sprintf("%.*f", error_decimal_digits, rounded_error)
    end

    return val_str, err_str
end


## -------------------------------------------------------------------------- ##
##                            Statistical routines                            ##
## -------------------------------------------------------------------------- ##

"""
    compute_stats(real_values, imag_values; auto_digits=true)

计算复数数据的统计量（平均值、误差和格式化字符串）。

参数:
- `real_values`: 实部值数组
- `imag_values`: 虚部值数组
- `auto_digits`: 是否自动确定有效数字 (默认: true)

返回:
- 包含以下字段的命名元组:
  - `mean_real`, `mean_imag`: 平均值（实部和虚部）
  - `err_real`, `err_imag`: 误差（实部和虚部）
  - `formatted_real`, `formatted_imag`: 格式化后的结果字符串
"""
function compute_stats(real_values, imag_values; auto_digits=true)
    # 计算统计量
    mean_real = mean(real_values)
    mean_imag = mean(imag_values)
    
    # 计算误差（完整精度和格式化精度）
    err_real = error(real_values, sigma=1, bessel=true, auto_digits=false)
    err_imag = error(imag_values, sigma=1, bessel=true, auto_digits=false)
    err_real_fmt = error(real_values, sigma=1, bessel=true, auto_digits=auto_digits)
    err_imag_fmt = error(imag_values, sigma=1, bessel=true, auto_digits=auto_digits)
    
    # 格式化显示
    formatted_real, formatted_real_err = format_value_error(mean_real, err_real_fmt)
    formatted_imag, formatted_imag_err = format_value_error(mean_imag, err_imag_fmt)
    
    return (;
        mean_real = mean_real,
        mean_imag = mean_imag,
        err_real = err_real,
        err_imag = err_imag,
        formatted_real = "$(formatted_real) ± $(formatted_real_err)",
        formatted_imag = "$(formatted_imag) ± $(formatted_imag_err)"
    )
end

"""
    compute_cr_m2_covariance(m2_value, m2_error, cr_value, scaling_factor; return_components=false)

计算相关比（Correlation Ratio）和缩放m²之间的协方差。

# 数学推导

给定独立测量量：
- A = C(0) = m² （Q点的结构因子）
- B = C(δk) （Q+δq点的结构因子）
- Cov(A, B) = 0 （假设独立测量）

定义：
- K = L^(-(1+η_φ)) （缩放常数）
- y = K·A （缩放后的m²）
- x = 1 - B/A （相关比）

使用误差传播公式：
```
σ_xy ≈ (∂x/∂A)(∂y/∂A)σ_A² + (∂x/∂B)(∂y/∂B)σ_B²
```

由于 ∂y/∂B = 0，第二项为零。计算偏导数：
- ∂y/∂A = K
- ∂x/∂A = B/A²

因此：
```
σ_xy = (B/A²)·K·σ_A² = y·(1-x)·(σ_A/A)²
```

最终公式：
```
σ_xy = y·(1-x)·(σ_C(0)/C(0))²
```

# 参数

- `m2_value::Real`: m² = C(0) 的值（Q点的结构因子）
- `m2_error::Real`: σ_C(0)，m²的误差
- `cr_value::Real`: x，相关比的值（通常为 1 - C(δk)/C(0)）
- `scaling_factor::Real`: 例如K = L^(1+η_φ)，缩放因子
- `return_components::Bool`: 是否返回中间计算值（默认：false）

# 返回值

如果 `return_components=false`（默认）：
- `(; covariance = σ_xy)` - 仅返回协方差值

如果 `return_components=true`：
- `(; covariance, scaled_m2, relative_error_m2)` - 返回协方差及中间值

# 示例

```julia
# 基本用法
m2 = 2.0
σ_m2 = 0.1
cr = 0.5
K = 0.8
result = compute_cr_m2_covariance(m2, σ_m2, cr, K)
println("协方差: ", result.covariance)

# 返回中间值用于调试
result = compute_cr_m2_covariance(m2, σ_m2, cr, K; return_components=true)
println("缩放m²: ", result.scaled_m2)
println("相对误差: ", result.relative_error_m2)
```

# 参考

- AFMCorrelationRatio: 计算反铁磁相关比
- CDWCorrelationRatio: 计算CDW相关比
- AFMStructureFactor: 计算反铁磁结构因子（m²）
"""
function compute_cr_m2_covariance(
    m2_value::Real,
    m2_error::Real,
    cr_value::Real,
    scaling_factor::Real;
    return_components::Bool=false
)
    # 输入验证
    @assert m2_value > 0 "m2_value must be positive, got $(m2_value)"
    @assert m2_error >= 0 "m2_error must be non-negative, got $(m2_error)"
    @assert scaling_factor > 0 "scaling_factor must be positive, got $(scaling_factor)"

    # 计算缩放后的m²: y = K·A
    scaled_m2 = scaling_factor * m2_value

    # 计算相对误差: σ_A/A
    relative_error = m2_error / m2_value

    # 计算协方差: σ_xy = y·(1-x)·(σ_A/A)²
    covariance = scaled_m2 * (1.0 - cr_value) * relative_error^2

    # 根据标志返回结果
    if return_components
        return (;
            covariance = covariance,
            scaled_m2 = scaled_m2,
            relative_error_m2 = relative_error
        )
    else
        return (; covariance = covariance)
    end
end

"""
    compute_cr_m2_covariance(scaled_m2, scaled_error, cr_value)

计算相关比和缩放m²之间的协方差（简化版，直接使用 scaled 值）。

由于 y = K·A 是线性缩放，相对误差保持不变：σ_y/y = σ_A/A
因此可以直接使用 scaled 值计算协方差。

# 参数
- `scaled_m2::Real`: Y = K·m²，缩放后的 m² 值
- `scaled_error::Real`: E = K·σ_m²，缩放后的误差
- `cr_value::Real`: x，相关比的值

# 返回值
- `(; covariance)` - 协方差值 σ_xy

# 示例
```julia
Y = 0.5    # scaled m²
E = 0.02   # scaled error
x = 0.3    # correlation ratio
result = compute_cr_m2_covariance(Y, E, x)
println("协方差: ", result.covariance)
```

# 参考
- `compute_cr_m2_covariance(m2_value, m2_error, cr_value, scaling_factor)`: 4参数版本
"""
function compute_cr_m2_covariance(
    scaled_m2::Real,
    scaled_error::Real,
    cr_value::Real
)
    @assert scaled_m2 > 0 "scaled_m2 must be positive, got $(scaled_m2)"
    @assert scaled_error >= 0 "scaled_error must be non-negative, got $(scaled_error)"

    relative_error = scaled_error / scaled_m2
    covariance = scaled_m2 * (1.0 - cr_value) * relative_error^2

    return (; covariance = covariance)
end
