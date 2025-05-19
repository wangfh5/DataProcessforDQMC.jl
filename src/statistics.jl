using Statistics
using DataFrames
using Printf

export error, round_error, format_value_error

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

"""
    mean_ns(data)
mean without spikes: Calculate the mean of the data after removing the largest and smallest value.
#!deprecated
"""
function mean_ns(data)
    sort_data = sort(data)
    mean(sort_data[2:end-1])
end

"""
    error_ns(data)
error without spikes: Calculate the error of the data after removing the largest and smallest value.
#!deprecated
"""
function error_ns(data)
    n = length(data)
    sort_data = sort(data)
    3*std(sort_data[2:end-1])/sqrt(n-2)
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
    format_value_error(value, error, error_sig_digits=1)
Format a value and its error with appropriate precision using scientific notation.

The error is formatted to the specified number of significant digits.
The value is rounded to match the precision of the error.

Arguments:
- `value`: The main value to format
- `error`: The error/uncertainty of the value
- `error_sig_digits`: Number of significant digits to use for the error (default: 1)

Returns:
- Tuple of (formatted_value, formatted_error) as strings in scientific notation

Examples:
```julia
# Default 1 significant digit for error
format_value_error(2.36738, 0.0023)     # Returns ("2.367e+00", "0.003e0")
format_value_error(2.36738, 0.00023)    # Returns ("2.3674e+00", "0.0003e0")

# With 2 significant digits for error
format_value_error(2.36738, 0.0023, 2)  # Returns ("2.3674e+00", "0.0023e0")

# Large numbers
format_value_error(2367.38, 23)         # Returns ("2.37e+03", "0.03e3")
format_value_error(2367.38, 23, 2)      # Returns ("2.367e+03", "0.023e3")
```
"""
function format_value_error(value::Number, error::Number, error_sig_digits::Int=1)

    # Step 1: Determine the digits of the error based on its significant digits and error
    # For example
    # - If error = 23456, error_sig_digits = 1, error_digits = - 4;
    # - If error = 23456, error_sig_digits = 2, error_digits = - 4 + 2 - 1 = - 3;
    # - If error = 2.3456, error_sig_digits = 1, error_digits = 0;
    # - If error = 2.3456, error_sig_digits = 3, error_digits = 0 + 3 - 1 = 2;
    # Handle the case where error is exactly zero
    if error == 0.0
        error_order = 0
        error_digits = error_sig_digits - 1
        rounded_error = 0.0
    else
        # For non-zero errors, calculate order of magnitude
        error_order = floor(log10(abs(error)))
        # Check if error_order is -Inf (can happen with very small numbers due to floating point precision)
        if isinf(error_order)
            # Use the smallest representable order for very small numbers
            error_order = -324  # Approximately the smallest exponent for Float64
        end
        error_digits = - Int(error_order - error_sig_digits + 1)
        rounded_error = round(error, RoundUp, sigdigits=error_sig_digits)
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

    # Step 4: Format the value and error as strings in scientific notation
    # For clarity, the exponents are printed to be the same for both values
    # Value is automatically formatted, while error is formatted manually
    # For example
    # - If value = 2.36738, error = 0.0023, error_sig_digits = 1, then we want 2.367e+00 and 0.003e+00
    # - If value = 2.36738, error = 0.0023, error_sig_digits = 2, then we want 2.3674e+00 and 0.0023e+00
    # - If value = 2367.38, error = 23, error_sig_digits = 1, then we want 2.37e+03 and 0.03e+03
    # - If value = 2367.38, error = 23, error_sig_digits = 2, then we want 2.367e+03 and 0.023e+03
    val_str = @sprintf("%.*e", value_sig_digits - 1, rounded_val)
    err_exponent = Int(value_order)  # Safe now because we've handled the zero case above

    # Handle the case where error is effectively zero
    if rounded_error == 0.0
        err_str = "0e$(err_exponent)"
    else
        err_magnitude = abs(rounded_error) / 10.0^err_exponent
        rounded_error_format = round(err_magnitude, sigdigits=error_sig_digits)
        err_str = "$(rounded_error_format)e$(err_exponent)"
    end

    return val_str, err_str
end


## -------------------------------------------------------------------------- ##
##                            Statistical routines                            ##
## -------------------------------------------------------------------------- ##

"""
    statistics_columns(df_tmp::DataFrame)
average and error of each column of a dataframe with only `bin` and data columns.
"""
function statistics_columns(df_tmp::DataFrame)
    df_tmp = bin_filter(df_tmp)
    # calculate the mean and error of each column
    mean_values = mean.(eachcol(df_tmp))
    error_values = error.(eachcol(df_tmp))
    # 创建一个存储统计信息的空 DataFrame
    df_stats = DataFrame()
    # 将每一列的平均值和标准差添加到新的 DataFrame 中
    for (colname, mean_val, error_val) in zip(names(df_tmp), mean_values, error_values)
        df_stats[!, Symbol(colname, "_mean")] = [mean_val]
        df_stats[!, Symbol(colname, "_error")] = [error_val]
    end
    return df_stats
end


"""
    statistics_columns_withparas(df_tmp::DataFrame, paras)
average and error of a dataframe with `bin` column, parameter columns specified by `paras` (symbols) and data columns.
"""
function statistics_columns_withparas(df_tmp; paras=nothing)
    df_tmp = bin_filter(df_tmp)
    # find the data names
    datavec = setdiff(Symbol.(names(df_tmp)), paras)
    datamat = reshape(datavec, 1, :)
    # group and combine
    gdf = groupby(df_tmp, paras)
    df_stat = combine(gdf, datamat .=> [mean, error])
    return df_stat
end
