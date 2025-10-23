#= 
单参数分析模块 (Single Parameter Analysis) 之关联比率分析

此模块提供了一组函数，用于计算不同结构因子的关联比率 （Correlation ratio）。

主要功能:
- 反铁磁相关比率分析 (AFMCorrelationRatio)
- CDW相关比率分析 (CDWCorrelationRatio)

=#

export AFMCorrelationRatio, CDWCorrelationRatio


"""    
    AFMCorrelationRatio(shift_point::Tuple{<:Real,<:Real}, Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                     filename::String="afm_sf_k.bin", source_file::String="spsm_k.bin", filedir::String=pwd();
                     force_rebuild::Bool=false, startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                     orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                     orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                     auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)

Calculate the correlation ratio R_{m^2} for antiferromagnetic structure factor, used to quantify disorder-order transitions.

## Formula
R_{m^2} = 1 - S_AFM(Q+δq) / S_AFM(Q)
where δq is a small shift in reciprocal space

## Parameters
- `shift_point`: Momentum space shift (δq_x, δq_y)
- `Q_point`: Antiferromagnetic vector Q (default: (0.0, 0.0))
- `filename`: Structure factor file name (default: "afm_sf_k.bin")
- `source_file`: Source file to generate structure factor if not exists (default: "spsm_k.bin")
- `filedir`: File directory (default: current directory)
- `force_rebuild`: Force rebuild structure factor file even if exists (default: false)
- `startbin`: Starting bin (default: 2)
- `endbin`: Ending bin (default: all bins)
- `dropmaxmin`: Number of max/min values to drop (default: 0)
- `orbital_columns`: Orbital column indices (default: [(3,4), (5,6), (7,8), (9,10)])
- `orbital_labels`: Orbital labels (default: ["AA", "AB", "BA", "BB"])
- `auto_digits`: Whether to automatically determine significant digits (default: true)
- `tolerance`: Tolerance for matching k-points (default: 1e-6)
- `verbose`: Whether to output detailed information (default: true)

## Returns
A named tuple with the following fields:
- `Q_point`: Antiferromagnetic vector Q
- `shift_point`: Momentum space shift (δq_x, δq_y)
- `Q_shifted`: Shifted k-point
- `S_AFM_Q`: Antiferromagnetic structure factor at Q
- `err_S_AFM_Q`: Error of the antiferromagnetic structure factor at Q
- `S_AFM_Q_shifted`: Antiferromagnetic structure factor at the shifted point
- `err_S_AFM_Q_shifted`: Error of the antiferromagnetic structure factor at the shifted point
- `correlation_ratio`: Calculated correlation ratio
- `err_correlation_ratio`: Error of the correlation ratio
- `formatted_correlation_ratio`: Formatted correlation ratio result

## Example
```julia
# Calculate correlation ratio with (0.25,0) shift
result = AFMCorrelationRatio((0.25, 0.0))
println("Correlation Ratio: \$(result.formatted_correlation_ratio)")
```

"""
function AFMCorrelationRatio(shift_point::Tuple{<:Real,<:Real}, Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                         filename::String="afm_sf_k.bin", filedir::String=pwd();
                         source_file::String="spsm_k.bin", force_rebuild::Bool=false, startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                         orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                         orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                         auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    
    # 计算Q点处的反铁磁结构因子
    if verbose
        println("计算Q点 $Q_point 处的反铁磁结构因子")
    end
    
    S_AFM_Q = AFMStructureFactor(Q_point, filename, filedir;
                               force_rebuild=force_rebuild, source_file=source_file, startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
                               auto_digits=auto_digits, tolerance=tolerance, verbose=verbose)
    
    # 计算偏移后的k点坐标
    shift_x, shift_y = shift_point
    Q_shifted = (Q_point[1] + shift_x, Q_point[2] + shift_y)
    
    if verbose
        println("\n计算偏移点 $Q_shifted 处的反铁磁结构因子")
        println("偏移量: ($(round(shift_x, digits=6)), $(round(shift_y, digits=6)))")
    end
    
    # 计算偏移点处的反铁磁结构因子（文件已在上一步生成，无需再次 rebuild）
    S_AFM_Q_shifted = AFMStructureFactor(Q_shifted, filename, filedir;
                                       force_rebuild=false, source_file=source_file, startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
                                       auto_digits=auto_digits, tolerance=tolerance, verbose=verbose)
    
    # 使用实部计算相关比（通常反铁磁结构因子主要是实部）
    ratio = S_AFM_Q_shifted.mean_real / S_AFM_Q.mean_real
    correlation_ratio = 1.0 - ratio
    
    # 误差传播 (使用一阶近似)：
    # 如果 R = 1 - S2/S1，则 δR ≈ sqrt((S2/S1^2 * δS1)^2 + (1/S1 * δS2)^2)
    S1 = S_AFM_Q.mean_real
    S2 = S_AFM_Q_shifted.mean_real
    dS1 = S_AFM_Q.err_real
    dS2 = S_AFM_Q_shifted.err_real
    
    err_ratio = sqrt((S2/(S1^2) * dS1)^2 + (1/S1 * dS2)^2)
    
    # 使用round_error函数保留一位有效数字
    rounded_err, _ = round_error(err_ratio, err_ratio/10)
    
    # 使用format_value_error函数格式化结果
    formatted_val, formatted_err = format_value_error(correlation_ratio, rounded_err, 1)
    formatted_correlation_ratio = "$(formatted_val) ± $(formatted_err)"
    
    # 打印结果
    if verbose
        println("\n自旋结构因子相关比（Correlation Ratio）分析:")
        println("---------------------------------------------------")
        println("Q点: $Q_point")
        println("偏移点: $Q_shifted (偏移量: $shift_point)")
        println("S_AFM(Q) = $(S_AFM_Q.formatted_real)")
        println("S_AFM(Q+δq) = $(S_AFM_Q_shifted.formatted_real)")
        println("相关比 R = 1 - S_AFM(Q+δq)/S_AFM(Q) = $formatted_correlation_ratio")
        println("---------------------------------------------------")
    end
    
    # 返回结果
    return (
        Q_point = Q_point,
        shift_point = shift_point,
        Q_shifted = Q_shifted,
        S_AFM_Q = S_AFM_Q.mean_real,
        err_S_AFM_Q = S_AFM_Q.err_real,
        S_AFM_Q_shifted = S_AFM_Q_shifted.mean_real,
        err_S_AFM_Q_shifted = S_AFM_Q_shifted.err_real,
        correlation_ratio = correlation_ratio,
        err_correlation_ratio = err_ratio,
        formatted_correlation_ratio = formatted_correlation_ratio
    )
end

"""    
    CDWCorrelationRatio(shift_point::Tuple{<:Real,<:Real}, Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                     filename::String="cdwpair_sf_k.bin", source_file::String="cdwpair_k.bin", filedir::String=pwd();
                     force_rebuild::Bool=false, startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                     orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                     orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                     auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)

Calculate the correlation ratio R_{m^2} for charge density wave structure factor, used to quantify disorder-order transitions.

## Formula
R_{m^2} = 1 - S_CDW(Q+δq) / S_CDW(Q)
where δq is a small shift in reciprocal space

## Parameters
- `shift_point`: Momentum space shift (δq_x, δq_y)
- `Q_point`: CDW ordering vector Q (default: (0.0, 0.0))
- `filename`: Structure factor file name (default: "cdwpair_sf_k.bin")
- `source_file`: Source file to generate structure factor if not exists (default: "cdwpair_k.bin")
- `filedir`: File directory (default: current directory)
- `force_rebuild`: Force rebuild structure factor file even if exists (default: false)
- `startbin`: Starting bin (default: 2)
- `endbin`: Ending bin (default: all bins)
- `dropmaxmin`: Number of max/min values to drop (default: 0)
- `orbital_columns`: Orbital column indices (default: [(3,4), (5,6), (7,8), (9,10)])
- `orbital_labels`: Orbital labels (default: ["AA", "AB", "BA", "BB"])
- `auto_digits`: Whether to automatically determine significant digits (default: true)
- `tolerance`: Tolerance for matching k-points (default: 1e-6)
- `verbose`: Whether to output detailed information (default: true)

## Returns
A named tuple with the following fields:
- `Q_point`: CDW ordering vector Q
- `shift_point`: Momentum space shift (δq_x, δq_y)
- `Q_shifted`: Shifted k-point
- `S_CDW_Q`: CDW structure factor at Q
- `err_S_CDW_Q`: Error of the CDW structure factor at Q
- `S_CDW_Q_shifted`: CDW structure factor at the shifted point
- `err_S_CDW_Q_shifted`: Error of the CDW structure factor at the shifted point
- `correlation_ratio`: Calculated correlation ratio
- `err_correlation_ratio`: Error of the correlation ratio
- `formatted_correlation_ratio`: Formatted correlation ratio result

## Example
```julia
# Calculate correlation ratio with (0.25,0) shift
result = CDWCorrelationRatio((0.25, 0.0))
println("Correlation Ratio: \$(result.formatted_correlation_ratio)")
```
"""
function CDWCorrelationRatio(shift_point::Tuple{<:Real,<:Real}, Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                         filename::String="cdwpair_sf_k.bin", filedir::String=pwd();
                         source_file::String="cdwpair_k.bin", force_rebuild::Bool=false, startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                         orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                         orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                         auto_digits::Bool=true, tolerance::Float64=1e-6, verbose::Bool=true)
    
    # Calculate CDW structure factor at Q point
    if verbose
        println("Calculating CDW structure factor at Q point $Q_point")
    end
    
    S_CDW_Q = CDWStructureFactor(Q_point, filename, filedir;
                               force_rebuild=force_rebuild, source_file=source_file, startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
                               auto_digits=auto_digits, tolerance=tolerance, verbose=verbose)
    
    # Calculate shifted k-point coordinates
    shift_x, shift_y = shift_point
    Q_shifted = (Q_point[1] + shift_x, Q_point[2] + shift_y)
    
    if verbose
        println("\nCalculating CDW structure factor at shifted point $Q_shifted")
        println("Shift: ($(round(shift_x, digits=6)), $(round(shift_y, digits=6)))")
    end
    
    # Calculate CDW structure factor at shifted point (file already generated, no need to rebuild)
    S_CDW_Q_shifted = CDWStructureFactor(Q_shifted, filename, filedir;
                                       force_rebuild=false, source_file=source_file, startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
                                       auto_digits=auto_digits, tolerance=tolerance, verbose=verbose)
    
    # Calculate correlation ratio using real part (typically CDW structure factor is mainly real)
    ratio = S_CDW_Q_shifted.mean_real / S_CDW_Q.mean_real
    correlation_ratio = 1.0 - ratio
    
    # Error propagation (using first-order approximation):
    # If R = 1 - S2/S1, then δR ≈ sqrt((S2/S1^2 * δS1)^2 + (1/S1 * δS2)^2)
    S1 = S_CDW_Q.mean_real
    S2 = S_CDW_Q_shifted.mean_real
    dS1 = S_CDW_Q.err_real
    dS2 = S_CDW_Q_shifted.err_real
    
    err_ratio = sqrt((S2/(S1^2) * dS1)^2 + (1/S1 * dS2)^2)
    
    # Round error to one significant digit
    rounded_err, _ = round_error(err_ratio, err_ratio/10)
    
    # Format result with value and error
    formatted_val, formatted_err = format_value_error(correlation_ratio, rounded_err, 1)
    formatted_correlation_ratio = "$(formatted_val) ± $(formatted_err)"
    
    # Print results
    if verbose
        println("\nCDW Structure Factor Correlation Ratio Analysis:")
        println("---------------------------------------------------")
        println("Q point: $Q_point")
        println("Shifted point: $Q_shifted (Shift: $shift_point)")
        println("S_CDW(Q) = $(S_CDW_Q.formatted_real)")
        println("S_CDW(Q+δq) = $(S_CDW_Q_shifted.formatted_real)")
        println("Correlation Ratio R = 1 - S_CDW(Q+δq)/S_CDW(Q) = $formatted_correlation_ratio")
        println("---------------------------------------------------")
    end
    
    # Return results
    return (
        Q_point = Q_point,
        shift_point = shift_point,
        Q_shifted = Q_shifted,
        S_CDW_Q = S_CDW_Q.mean_real,
        err_S_CDW_Q = S_CDW_Q.err_real,
        S_CDW_Q_shifted = S_CDW_Q_shifted.mean_real,
        err_S_CDW_Q_shifted = S_CDW_Q_shifted.err_real,
        correlation_ratio = correlation_ratio,
        err_correlation_ratio = err_ratio,
        formatted_correlation_ratio = formatted_correlation_ratio
    )
end
