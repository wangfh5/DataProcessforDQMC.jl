#=
多参数分析功能 (Multiple Parameter Analysis) 之关联比率分析
=#

# 导出多参数分析函数
export analyze_correlation_ratio_multi_parameter,
       analyze_AFM_correlation_ratio_multi_parameter,
       analyze_CDW_correlation_ratio_multi_parameter


"""
    analyze_correlation_ratio_multi_parameter(correlation_ratio_function::Function,
                                           base_dir::AbstractString=pwd(); 
                                           shift_point::Tuple{<:Real,<:Real}=(0.25, 0.0),
                                           Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                           filename::String="afm_sf_k.bin",
                                           source_file::String="spsm_k.bin",
                                           result_columns::Vector{Symbol}=[:correlation_ratio, :err_correlation_ratio],
                                           result_prefix::String="R",
                                           startbin::Int=2,
                                           endbin::Union{Int,Nothing}=nothing,
                                           dropmaxmin::Int=0,
                                           auto_digits::Bool=true,
                                           tolerance::Float64=1e-6,
                                           verbose::Bool=false,
                                           filter_options::Union{Dict, NamedTuple}=Dict(),
                                           pattern::Regex=r"^proj_fft_honeycomb") -> DataFrame

Generic multi-parameter correlation ratio analysis function, suitable for analyzing different types of correlation ratios.

# Parameters
- `correlation_ratio_function::Function`: Function for analyzing correlation ratio, such as AFMCorrelationRatio or CDWCorrelationRatio
- `base_dir::AbstractString`: Base directory containing parameter directories (default: current directory)
- `shift_point::Tuple{<:Real,<:Real}`: Momentum space shift (δq_x, δq_y) (default: (0.25, 0.0))
- `Q_point::Tuple{<:Real,<:Real}`: Ordering vector Q (default: (0.0, 0.0))
- `filename::String`: Structure factor file name (default: "afm_sf_k.bin")
- `source_file::String`: Source file to generate structure factor if not exists (default: "spsm_k.bin")
- `result_columns::Vector{Symbol}`: Columns to include in results (default: [:correlation_ratio, :err_correlation_ratio])
- `result_prefix::String`: Prefix for result column names (default: "R")
- `startbin::Int`: Starting bin for analysis (default: 2)
- `endbin::Union{Int,Nothing}`: Ending bin for analysis (default: all bins)
- `dropmaxmin::Int`: Number of max/min values to drop (default: 0)
- `auto_digits::Bool`: Whether to automatically determine significant digits (default: true)
- `tolerance::Float64`: Tolerance for matching k-points (default: 1e-6)
- `verbose::Bool`: Whether to output detailed information (default: false)
- `filter_options::Union{Dict, NamedTuple}`: Options for filtering parameter directories (default: empty Dict)
- `pattern::Regex`: Pattern for matching parameter directories (default: r"^proj_fft_honeycomb")

# Returns
- `DataFrame`: DataFrame containing parameters and correlation ratio results
"""
function analyze_correlation_ratio_multi_parameter(correlation_ratio_function::Function,
                                           base_dir::AbstractString=pwd(); 
                                           shift_point::Tuple{<:Real,<:Real}=(0.25, 0.0),
                                           Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                           filename::String="afm_sf_k.bin",
                                           source_file::String="spsm_k.bin",
                                           result_columns::Vector{Symbol}=[:correlation_ratio, :err_correlation_ratio],
                                           result_prefix::String="R",
                                           startbin::Int=2,
                                           endbin::Union{Int,Nothing}=nothing,
                                           dropmaxmin::Int=0,
                                           auto_digits::Bool=true,
                                           tolerance::Float64=1e-6,
                                           verbose::Bool=false,
                                           filter_options::Union{Dict, NamedTuple}=Dict(),
                                           pattern::Regex=r"^proj_fft_honeycomb")
    # Scan parameter directories
    param_dirs = scan_parameter_directories(base_dir; filter_options=filter_options, pattern=pattern)
    
    if isempty(param_dirs)
        @warn "No parameter directories found in $(base_dir)"
        return DataFrame()
    end
    
    # Create a DataFrame with basic columns
    df = DataFrame()
    
    # Add parameter columns based on the first directory
    params = extract_parameters_from_dirname(basename(first(param_dirs)))
    for param_name in keys(params)
        df[!, param_name] = typeof(params[param_name])[]
    end
    
    # Add result columns
    for col in result_columns
        col_name = Symbol("$(result_prefix)_$(col)")
        df[!, col_name] = Float64[]
    end
    
    # Add formatted result column
    df[!, Symbol("$(result_prefix)_formatted")] = String[]
    
    # Process each parameter directory
    for dir_path in param_dirs
        dir_name = basename(dir_path)
        # Extract parameters
        params = extract_parameters_from_dirname(dir_name)
        
        # Analyze correlation ratio
        try
            result = correlation_ratio_function(
                shift_point, Q_point, 
                filename, dir_path;
                source_file=source_file,
                startbin=startbin, 
                endbin=endbin, 
                dropmaxmin=dropmaxmin,
                auto_digits=auto_digits, 
                tolerance=tolerance, 
                verbose=verbose
            )
            
            # Create a new row
            row = Dict{Symbol, Any}()
            
            # Add parameters
            for (param_name, param_value) in params
                row[param_name] = param_value
            end
            
            # Add results
            for col in result_columns
                col_name = Symbol("$(result_prefix)_$(col)")
                row[col_name] = getproperty(result, col)
            end
            
            # Add formatted result
            row[Symbol("$(result_prefix)_formatted")] = result.formatted_correlation_ratio
            
            # Add row to DataFrame
            push!(df, row)
            
        catch e
            @warn "Error analyzing directory $(dir_name): $(e)"
        end
    end
    
    # Sort DataFrame by parameters
    sort!(df, [:b, :U, :gw])
    
    return df
end

"""
    analyze_AFM_correlation_ratio_multi_parameter(base_dir::AbstractString=pwd();
                                               shift_point::Tuple{<:Real,<:Real}=(0.25, 0.0),
                                               Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="afm_sf_k.bin",
                                               source_file::String="spsm_k.bin",
                                               startbin::Int=2,
                                               endbin::Union{Int,Nothing}=nothing,
                                               dropmaxmin::Int=0,
                                               auto_digits::Bool=true,
                                               tolerance::Float64=1e-6,
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict(),
                                               pattern::Regex=r"^proj_fft_honeycomb") -> DataFrame

Analyze AFM correlation ratio across multiple parameter directories.

# Parameters
- `base_dir::AbstractString`: Base directory containing parameter directories (default: current directory)
- `shift_point::Tuple{<:Real,<:Real}`: Momentum space shift (δq_x, δq_y) (default: (0.25, 0.0))
- `Q_point::Tuple{<:Real,<:Real}`: AFM ordering vector Q (default: (0.0, 0.0))
- `filename::String`: Structure factor file name (default: "afm_sf_k.bin")
- `source_file::String`: Source file to generate structure factor if not exists (default: "spsm_k.bin")
- `startbin::Int`: Starting bin for analysis (default: 2)
- `endbin::Union{Int,Nothing}`: Ending bin for analysis (default: all bins)
- `dropmaxmin::Int`: Number of max/min values to drop (default: 0)
- `auto_digits::Bool`: Whether to automatically determine significant digits (default: true)
- `tolerance::Float64`: Tolerance for matching k-points (default: 1e-6)
- `verbose::Bool`: Whether to output detailed information (default: false)
- `filter_options::Union{Dict, NamedTuple}`: Options for filtering parameter directories (default: empty Dict)
- `pattern::Regex`: Pattern for matching parameter directories (default: r"^proj_fft_honeycomb")

# Returns
- `DataFrame`: DataFrame containing parameters and AFM correlation ratio results
"""
function analyze_AFM_correlation_ratio_multi_parameter(base_dir::AbstractString=pwd();
                                               shift_point::Tuple{<:Real,<:Real}=(0.25, 0.0),
                                               Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="afm_sf_k.bin",
                                               source_file::String="spsm_k.bin",
                                               startbin::Int=2,
                                               endbin::Union{Int,Nothing}=nothing,
                                               dropmaxmin::Int=0,
                                               auto_digits::Bool=true,
                                               tolerance::Float64=1e-6,
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict(),
                                               pattern::Regex=r"^proj_fft_honeycomb")
    return analyze_correlation_ratio_multi_parameter(
        AFMCorrelationRatio,
        base_dir;
        shift_point=shift_point,
        Q_point=Q_point,
        filename=filename,
        source_file=source_file,
        result_columns=[:correlation_ratio, :err_correlation_ratio],
        result_prefix="R_AFM",
        startbin=startbin,
        endbin=endbin,
        dropmaxmin=dropmaxmin,
        auto_digits=auto_digits,
        tolerance=tolerance,
        verbose=verbose,
        filter_options=filter_options,
        pattern=pattern
    )
end

"""
    analyze_CDW_correlation_ratio_multi_parameter(base_dir::AbstractString=pwd();
                                               shift_point::Tuple{<:Real,<:Real}=(0.25, 0.0),
                                               Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="cdwpair_sf_k.bin",
                                               source_file::String="cdwpair_k.bin",
                                               startbin::Int=2,
                                               endbin::Union{Int,Nothing}=nothing,
                                               dropmaxmin::Int=0,
                                               auto_digits::Bool=true,
                                               tolerance::Float64=1e-6,
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict(),
                                               pattern::Regex=r"^proj_fft_honeycomb") -> DataFrame

Analyze CDW correlation ratio across multiple parameter directories.

# Parameters
- `base_dir::AbstractString`: Base directory containing parameter directories (default: current directory)
- `shift_point::Tuple{<:Real,<:Real}`: Momentum space shift (δq_x, δq_y) (default: (0.25, 0.0))
- `Q_point::Tuple{<:Real,<:Real}`: CDW ordering vector Q (default: (0.0, 0.0))
- `filename::String`: Structure factor file name (default: "cdwpair_sf_k.bin")
- `source_file::String`: Source file to generate structure factor if not exists (default: "cdwpair_k.bin")
- `startbin::Int`: Starting bin for analysis (default: 2)
- `endbin::Union{Int,Nothing}`: Ending bin for analysis (default: all bins)
- `dropmaxmin::Int`: Number of max/min values to drop (default: 0)
- `auto_digits::Bool`: Whether to automatically determine significant digits (default: true)
- `tolerance::Float64`: Tolerance for matching k-points (default: 1e-6)
- `verbose::Bool`: Whether to output detailed information (default: false)
- `filter_options::Union{Dict, NamedTuple}`: Options for filtering parameter directories (default: empty Dict)
- `pattern::Regex`: Pattern for matching parameter directories (default: r"^proj_fft_honeycomb")

# Returns
- `DataFrame`: DataFrame containing parameters and CDW correlation ratio results
"""
function analyze_CDW_correlation_ratio_multi_parameter(base_dir::AbstractString=pwd();
                                               shift_point::Tuple{<:Real,<:Real}=(0.25, 0.0),
                                               Q_point::Tuple{<:Real,<:Real}=(0.0, 0.0),
                                               filename::String="cdwpair_sf_k.bin",
                                               source_file::String="cdwpair_k.bin",
                                               startbin::Int=2,
                                               endbin::Union{Int,Nothing}=nothing,
                                               dropmaxmin::Int=0,
                                               auto_digits::Bool=true,
                                               tolerance::Float64=1e-6,
                                               verbose::Bool=false,
                                               filter_options::Union{Dict, NamedTuple}=Dict(),
                                               pattern::Regex=r"^proj_fft_honeycomb")
    return analyze_correlation_ratio_multi_parameter(
        CDWCorrelationRatio,
        base_dir;
        shift_point=shift_point,
        Q_point=Q_point,
        filename=filename,
        source_file=source_file,
        result_columns=[:correlation_ratio, :err_correlation_ratio],
        result_prefix="R_CDW",
        startbin=startbin,
        endbin=endbin,
        dropmaxmin=dropmaxmin,
        auto_digits=auto_digits,
        tolerance=tolerance,
        verbose=verbose,
        filter_options=filter_options,
        pattern=pattern
    )
end
