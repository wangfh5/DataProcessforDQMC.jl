using Statistics
using DataFrames
using Printf

export error, round_error, format_value_error
export iqr_fence_filter

## -------------------------------------------------------------------------- ##
##                                   Filters                                  ##
## -------------------------------------------------------------------------- ##

"""
    bin_filter(df::DataFrame)
Remove the first `ntherm` bins from `df`. `df` should have a column named `bin`.
"""
@inline function bin_filter(df::DataFrame;ntherm=1::Int64)
    # filter data
    df_filter = df[df.bin .> ntherm,:]
    # remove the bin column
    select!(df_filter, Not(:bin))
    df_filter
end

"""
    k_filter(df::DataFrame, k::Tuple{Float64,Float64})
Filter the data at k = (kx,ky).
"""
@inline function k_filter(df::DataFrame, k::Tuple{Float64,Float64})
    # filter data
    df_filter = df[(df.kx .== k[1]) .& (df.ky .== k[2]),:]
    # remove the kx, ky column
    select!(df_filter, Not(:kx))
    select!(df_filter, Not(:ky))
    df_filter
end

"""
Remove the max and min of the data
"""
function filter(onecolumndata::Array{Float64,1}, dropnum::Int64=1)
    onecolumndata_sort = sort(onecolumndata)
    return onecolumndata_sort[1+dropnum:end-dropnum]
end
function filter(multicolumnsdata::Array{Float64,2}, dropnum::Int64=1)
    return hcat([filter(multicolumnsdata[:,i], dropnum) for i in 1:size(multicolumnsdata,2)]...)
end
function filter(data::Vector{Complex{T}}, dropnum::Int64=1) where T <: AbstractFloat
    if length(data) <= 2*dropnum
        return data
    end

    # Sort by magnitude
    sorted_indices = sortperm(abs.(data))

    # Remove the dropnum smallest and largest values
    return data[sorted_indices[(dropnum+1):(end-dropnum)]]
end

"""
    iqr_fence_filter(values; k=10.0, min_n=8)

Apply Tukey's IQR fence to remove extreme outliers.

Given samples `x₁,…,xₙ`, define the quartiles:

    Q₁ = quantile(x, 0.25),   Q₃ = quantile(x, 0.75),   IQR = Q₃ − Q₁

and the Tukey fences:

    L = Q₁ − k·IQR,    U = Q₃ + k·IQR

We keep `S = { xᵢ | L ≤ xᵢ ≤ U }` and treat the rest as outliers.

This function is designed for situations where **a tiny number of extreme algorithmic outliers**
(e.g. pathological bootstrap fits) can destroy plain mean/std.

# Behaviour / edge cases
- Non-finite values (`NaN`, `Inf`) and `missing` are ignored before computing quantiles.
- If the number of finite samples is `< min_n`, no fence is applied (returns input finite values).
- If `IQR == 0` or non-finite, no fence is applied.
- If the fence would remove all samples, we fall back to the original finite values.

# Returns
`(filtered_values::Vector{Float64}, n_removed::Int)`

where `n_removed` counts only the number of values removed by the fence (not NaN/Inf/missing).

# Parameters
- `k`: fence multiplier. Typical boxplot values: 1.5 (mild) or 3 (extreme).
       Use large `k` (e.g. 10) if you only want to remove **very extreme** outliers.
- `min_n`: minimum number of samples required to apply the fence.
"""
function iqr_fence_filter(values; k::Real=10.0, min_n::Int=8)
    # Collect finite numeric values; ignore missing/NaN/Inf.
    vals = Float64[]
    for v in values
        if ismissing(v)
            continue
        end
        x = Float64(v)
        isfinite(x) || continue
        push!(vals, x)
    end

    n = length(vals)
    if n < min_n
        return vals, 0
    end

    q1 = quantile(vals, 0.25)
    q3 = quantile(vals, 0.75)
    iqr = q3 - q1
    if !(isfinite(iqr) && iqr > 0)
        return vals, 0
    end

    lo = q1 - Float64(k) * iqr
    hi = q3 + Float64(k) * iqr

    # NOTE: This module defines its own `filter(x::Vector, dropnum)` for trimming
    # extrema, so use a comprehension here (or `Base.filter`) to avoid name clash.
    filtered = [x for x in vals if lo <= x <= hi]
    if isempty(filtered)
        return vals, 0
    end
    return filtered, n - length(filtered)
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
    # Step 1: Determine the digits of the error based on its significant digits and error
    if error == 0.0
        # Handle the case where error is exactly zero
        rounded_error = 0.0
        error_order = 0
        error_digits = error_sig_digits - 1
    else
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
    end

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
    else
        throw(ArgumentError("format must be :scientific or :decimal, got $format"))
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
