#= 
单参数分析模块 (Single Parameter Analysis)

此模块提供了一组函数，用于分析单个参数文件夹中的DQMC模拟数据。
通常在模拟输出目录中直接运行，处理该目录下的 .bin 文件。

主要功能:
- 能量分析 (EnergyAnalysis)
- 单轨道关联函数分析 (CorrelationAnalysis)
- 多轨道关联函数分析 (MultiOrbitalCorrelationAnalysis)
- Renyi负值计算 (RenyiNegativity, RenyiNegativity_all)

使用示例:
1. 在Julia REPL中:
   using DataProcessforDQMC
   EnergyAnalysis("energy.bin", ".")

2. 作为脚本运行:
   julia -e 'using DataProcessforDQMC; EnergyAnalysis("energy.bin", ".")'
=#

using Printf
using Statistics
using DelimitedFiles

export EnergyAnalysis, CorrelationAnalysis, MultiOrbitalCorrelationAnalysis
export RenyiNegativity, RenyiNegativity_all

## -------------------------------------------------------------------------- ##
##                              Renyi Negativity                              ##
## -------------------------------------------------------------------------- ##

function RenyiNegativity(filename::String, filedir::String=pwd();
    printLA=nothing, startbin::Int=3, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=1)
    # Add .bin extension if not present
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end

    # 打开文件
    filepath = joinpath(filedir, filename)

    # 使用 readdlm 读取整个文件
    data = readdlm(filepath, Float64)
    # Apply start and end bin selection
    data = data[startbin:end,:]
    if !isnothing(endbin)
        data = data[1:endbin,:]
    end
    # remove the max and min
    data = filter(data, dropmaxmin)

    # deduce L and rank from the data
    L = size(data,2) ÷ 2 - 1
    # example filename: expRenyiN3.bin, expRenyiN4_TW.bin
    rank = parse(Int, split(filename, "expRenyiN")[2][1])
    quantity_name = split(filename, ".bin")[1]

    # 输出结果
    expRenyiN = zeros(L+1,2)
    RenyiN = zeros(L+1,2)
    for i in 1:L+1
        vectmp = data[:,2i-1] # take real part
        if mean(vectmp) > 0.0
            mn = mean(vectmp)

            # Calculate error with auto_digits
            err = error(vectmp, sigma=1, bessel=true, auto_digits=true)

            expRenyiN[i,:] = [mn, err]
            RenyiN[i,:] = [log(mn)/(1-rank), abs(err/mn/(1-rank))]

            if (printLA == nothing) || ((i-1) ∈ printLA)
                # Format values for display
                formatted_mn, formatted_err = format_value_error(mn, err)
                renyi_val = RenyiN[i,1]
                renyi_err = RenyiN[i,2]
                formatted_renyi, formatted_renyi_err = format_value_error(renyi_val, renyi_err)

                println("$(quantity_name) $(i-1) $formatted_mn ± $formatted_err $formatted_renyi ± $formatted_renyi_err")
            end
        end
    end
    # @show expRenyiN
    # @show RenyiN
    return expRenyiN, RenyiN
end

function RenyiNegativity_all(filedir::String=pwd();maxrank::Int=4, kwargs...)
    # find all the files with name expRenyiN*.bin or expRenyiN*_TW.bin, where * is an integer
    filenames = Base.filter(x->occursin(r"expRenyiN\d+\.bin", x) || occursin(r"expRenyiN\d+_TW\.bin", x), readdir(filedir))
    # sort the filenames
    filenames = sort(filenames)
    # get the number of files
    nfiles = length(filenames)
    # for each file, create an array to store the average and error of Renyi negativity
    # and then append the result to an opened JLD2 file
    file = jldopen("RenyiNall.jld2", "w") do file
        for i in 1:nfiles
            filename = filenames[i]
            suffix = split(filename, "expRenyiN")[2][1:end-4]
            rank = parse(Int, suffix[1])
            if rank <= maxrank
                expRenyiN, RenyiN = RenyiNegativity(filename, filedir;kwargs...)
                # save the result to a JLD2 file
                file["expRenyiN$suffix"] = expRenyiN
                file["RenyiN$suffix"] = RenyiN
            end
        end
    end
end

## -------------------------------------------------------------------------- ##
##                                   Energy                                   ##
## -------------------------------------------------------------------------- ##

"""
    EnergyAnalysis(filename="energy.bin", filedir=pwd();
                   startbin=2, endbin=nothing, dropmaxmin=0,
                   columns=[5,7,9,11], labels=["E_kin", "E_pot", "E_tot", "E_tot^2"],
                   verbose=true)

Analyze energy.bin file and calculate Monte Carlo averages and errors for specified columns.

Arguments:
- `filename`: Name of the energy file (default: "energy.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 2)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 0)
- `columns`: Array of column indices to analyze (default: [5,7,9,11])
- `labels`: Array of labels for the columns (default: ["E_kin", "E_pot", "E_tot", "E_tot^2"])
- `verbose`: Whether to print results to console (default: true)

Returns:
- Dictionary with column labels as keys and [mean, error] arrays as values

Example:
```julia
# Analyze energy.bin in current directory
results = EnergyAnalysis()

# Analyze specific columns with custom labels
results = EnergyAnalysis(columns=[5,7], labels=["Kinetic", "Potential"])

# Access results
kinetic_energy = results["Kinetic"][1]  # mean value
kinetic_error = results["Kinetic"][2]   # error
```
"""
function EnergyAnalysis(filename::String="energy.bin", filedir::String=pwd();
                       startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                       columns::Vector{Int}=[5,7,9],
                       labels::Vector{String}=["E_kin", "d_occ", "E_tot"],
                       verbose::Bool=true)
    # Validate inputs
    if length(columns) != length(labels)
        error("Number of columns ($(length(columns))) must match number of labels ($(length(labels)))")
    end

    # Add .bin extension if not present
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end

    # Open and read the file
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Read data
    data = readdlm(filepath, Float64)

    # Check if file has enough columns
    if size(data, 2) < maximum(columns)
        error("File has only $(size(data, 2)) columns, but analysis requested column $(maximum(columns))")
    end

    # Apply start and end bin selection
    data = data[startbin:end,:]
    if !isnothing(endbin)
        data = data[1:endbin,:]
    end

    # Remove the max and min values
    data = filter(data, dropmaxmin)

    # Calculate statistics for each column
    results = Dict{String, Vector{Float64}}()

    if verbose
        println("\nMonte Carlo analysis of $filename:")
        println("----------------------------------------")
        println("Using $(size(data, 1)) measurements after filtering")
        println("----------------------------------------")
    end

    for (i, col) in enumerate(columns)
        label = labels[i]
        values = data[:, col]

        # Calculate mean and error
        mean_val = mean(values)

        # Calculate error with appropriate precision
        err_val = error(values, sigma=1, bessel=true, auto_digits=true)

        # Store results
        results[label] = [mean_val, err_val]

        if verbose
            # Format value and error for display
            formatted_val, formatted_err = format_value_error(mean_val, err_val)
            println("$label = $formatted_val ± $formatted_err")
        end
    end

    if verbose
        println("----------------------------------------\n")
    end

    return results
end

## -------------------------------------------------------------------------- ##
##                      Real-space Correlation functions                      ##
## -------------------------------------------------------------------------- ##

## ---------------- Helper functions for correlation analysis --------------- ##

"""Calculate Euclidean distance between two points"""
function euclidean_distance(x, y)
    return sqrt(x^2 + y^2)
end

"""Calculate i coordinate from imj when j is fixed at (1,1)"""
function calculate_i_coord(imj_x, imj_y, Lx, Ly)
    i_x = mod1(1 + imj_x, Lx)  # mod1 ensures result is in range 1:Lx
    i_y = mod1(1 + imj_y, Ly)  # mod1 ensures result is in range 1:Ly
    return (i_x, i_y)
end

"""Format a row for correlation output table"""
function format_correlation_row(coord, distance, i_coord, real_data, imag_data; col_width=25)
    coord_str = @sprintf("(%d,%d)", coord[1], coord[2])
    dist_str = @sprintf("%.3f", distance)
    i_coord_str = @sprintf("(%d,%d)", i_coord[1], i_coord[2])

    return @sprintf("%-15s %-15s %-15s", coord_str, dist_str, i_coord_str) *
           " " * lpad(real_data, col_width) *
           " " * lpad(imag_data, col_width)
end

"""Format a header for correlation output table"""
function format_correlation_header(col_width=25)
    return @sprintf("%-15s %-15s %-15s", "imj", "distance", "i [to j=(1,1)]") *
           " " * lpad("real part", col_width) *
           " " * lpad("imag part", col_width)
end

"""Process correlation data for a single coordinate"""
function process_correlation_values(values::Vector{Complex{T}}) where T <: AbstractFloat
    # Extract real and imaginary parts
    real_values = real.(values)
    imag_values = imag.(values)

    # Calculate mean for real and imaginary parts
    mean_real = mean(real_values)
    mean_imag = mean(imag_values)

    # Calculate error for real and imaginary parts (with full precision)
    # 不使用 auto_digits，保留全精度
    err_real = error(real_values, sigma=1, bessel=true, auto_digits=false)
    err_imag = error(imag_values, sigma=1, bessel=true, auto_digits=false)

    # 为了显示格式化，单独计算格式化后的误差值
    err_real_formatted = error(real_values, sigma=1, bessel=true, auto_digits=true)
    err_imag_formatted = error(imag_values, sigma=1, bessel=true, auto_digits=true)

    # Format values with appropriate precision for display only
    formatted_real, formatted_real_err = format_value_error(mean_real, err_real_formatted)
    formatted_imag, formatted_imag_err = format_value_error(mean_imag, err_imag_formatted)

    # Return results: 
    # 1. 全精度的复数值用于后续计算
    # 2. 格式化的字符串用于显示
    return [Complex(mean_real, mean_imag), Complex(err_real, err_imag)],
           "$(formatted_real) ± $(formatted_real_err)",
           "$(formatted_imag) ± $(formatted_imag_err)"
end

"""
    CorrelationAnalysis(filename="spsm_r.bin", filedir=pwd();
                        startbin=2, endbin=nothing, dropmaxmin=0,
                        real_column=3, imag_column=4,
                        auto_digits=true,
                        verbose=true)

Analyze correlation function data files like `spsm_r.bin` or `nn_r.bin`, where the first two columns
represent imj coordinates (distance vector between points i and j), and calculate Monte Carlo averages.
#TODO for Honeycomb: Add support for honeycomb lattice with multiple sublattices

Arguments:
- `filename`: Name of the correlation file (default: "spsm_r.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 2)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 0)
- `real_column`: Column index containing the real part of correlation values (default: 3)
- `imag_column`: Column index containing the imaginary part of correlation values (default: 4)
- `auto_digits`: Whether to automatically determine significant digits based on error of error (default: true)
- `verbose`: Whether to print results to console (default: true)

Returns:
- Dictionary with imj coordinate pairs as keys and [mean, error] arrays as values,
  where both mean and error are complex numbers

Example:
```julia
# Analyze spsm_r.bin in current directory
results = CorrelationAnalysis()

# Analyze nn_r.bin with custom columns
results = CorrelationAnalysis("nn_r.bin", real_column=3, imag_column=4)

# Use automatic determination of significant digits
results = CorrelationAnalysis(auto_digits=true)

# Access results
corr_1_1 = results[(1,1)][1]  # mean value at imj=(1,1)
corr_1_1_real = real(corr_1_1)  # real part
corr_1_1_imag = imag(corr_1_1)  # imaginary part
```
"""
function CorrelationAnalysis(filename::String="spsm_r.bin", filedir::String=pwd();
                            startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                            real_column::Int=3, imag_column::Int=4,
                            auto_digits::Bool=true,
                            verbose::Bool=true)
    # Add .bin extension if not present
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end

    # Open and read the file
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Read data
    data = readdlm(filepath, Float64)

    # Infer lattice size from data - maximum values of first and second columns
    Lx = Int(maximum(data[:, 1]))
    Ly = Int(maximum(data[:, 2]))

    # Analyze the structure of the data to determine the bin structure
    nrows = size(data, 1)
    unique_coords = Set([(Int(data[i,1]), Int(data[i,2])) for i in 1:nrows])
    n_coords = length(unique_coords)
    n_bins = nrows ÷ n_coords

    if verbose
        println("File structure analysis:")
        println("  Total rows: $nrows")
        println("  Unique coordinates: $n_coords")
        println("  Estimated bins: $n_bins")
    end

    # Group data by imj coordinates
    coord_groups = Dict{Tuple{Int,Int}, Vector{Complex{Float64}}}()

    # Initialize the coordinate groups
    for coord in unique_coords
        coord_groups[coord] = Complex{Float64}[]
    end

    # Process data based on bin structure
    if n_bins > 1
        # Process multi-bin data
        selected_bins, filtered_bins = process_bins(data, n_bins, n_coords, startbin, endbin, dropmaxmin, real_column, verbose)

        # Extract values for each coordinate from filtered bins
        for coord in unique_coords
            values = extract_coord_values(data, coord, filtered_bins, n_coords, real_column, imag_column)
            coord_groups[coord] = values
        end
    else
        # We only have one bin - can't do bin selection or filtering
        if verbose
            println("Warning: Only one bin found. Cannot perform bin selection or filtering.")
        end

        # Process each row as a single measurement
        for i in 1:nrows
            imj_x = Int(data[i, 1])
            imj_y = Int(data[i, 2])
            real_val = data[i, real_column]
            imag_val = data[i, imag_column]
            value = Complex(real_val, imag_val)

            key = (imj_x, imj_y)
            coord_groups[key] = [value]
        end
    end

    # Calculate statistics for each group
    results = Dict{Tuple{Int,Int}, Vector{Complex{Float64}}}()

    # Prepare output table
    if verbose
        print_correlation_header(filename, Lx, Ly, n_bins, n_bins > 1 ? length(selected_bins) : 0,
                               n_bins > 1 ? length(filtered_bins) : 0, dropmaxmin, coord_groups,
                               startbin, endbin)
    end

    # Sort coordinates for ordered output
    sorted_keys = sort(collect(keys(coord_groups)),
                      by = k -> (euclidean_distance(k[1], k[2]), k[1], k[2]))

    # Process each coordinate
    for key in sorted_keys
        imj_x, imj_y = key
        values = coord_groups[key]

        # Process values and get statistics
        result, real_data, imag_data = process_correlation_values(values)
        results[key] = result

        # Print results if verbose
        if verbose
            # Calculate distance and i coordinate
            distance = euclidean_distance(imj_x, imj_y)
            i_coord = calculate_i_coord(imj_x, imj_y, Lx, Ly)

            # Format and print row
            row = format_correlation_row((imj_x, imj_y), distance, i_coord, real_data, imag_data)
            println(row)
        end
    end

    if verbose
        println("----------------------------------------\n")
    end

    return results
end

"""Process bins and apply filtering"""
function process_bins(data, n_bins, n_coords, startbin, endbin, dropmaxmin, real_column, verbose)
    # First, identify which rows belong to which bin
    bin_rows = Dict{Int, Vector{Int}}()

    # Assuming coordinates repeat in the same order for each bin
    for bin_idx in 1:n_bins
        start_row = (bin_idx - 1) * n_coords + 1
        end_row = bin_idx * n_coords
        bin_rows[bin_idx] = collect(start_row:end_row)
    end

    # Apply bin selection
    selected_bins = collect(startbin:min(isnothing(endbin) ? n_bins : endbin, n_bins))

    # Apply filtering if needed
    if dropmaxmin > 0 && length(selected_bins) > 2*dropmaxmin
        # Sort bins by some metric (e.g., average value across all coordinates)
        bin_avg_values = Float64[]
        for bin_idx in selected_bins
            bin_sum = 0.0
            for row_idx in bin_rows[bin_idx]
                bin_sum += data[row_idx, real_column]
            end
            push!(bin_avg_values, bin_sum / n_coords)
        end

        # Sort bins by their average values
        sorted_indices = sortperm(bin_avg_values)

        # Remove the lowest and highest bins
        filtered_bins = selected_bins[sorted_indices[(dropmaxmin+1):(end-dropmaxmin)]]
    else
        filtered_bins = selected_bins
    end

    if verbose
        println("  Selected bins: $(length(selected_bins)) out of $n_bins")
        println("  Filtered bins: $(length(filtered_bins)) after removing $dropmaxmin max/min")
    end

    return selected_bins, filtered_bins
end

"""Extract values for a specific coordinate from filtered bins"""
function extract_coord_values(data, coord, filtered_bins, n_coords, real_column, imag_column)
    values = Complex{Float64}[]

    # For each filtered bin, find the row with this coordinate
    for bin_idx in filtered_bins
        start_row = (bin_idx - 1) * n_coords + 1
        end_row = bin_idx * n_coords

        for row_idx in start_row:end_row
            if (Int(data[row_idx, 1]), Int(data[row_idx, 2])) == coord
                real_val = data[row_idx, real_column]
                imag_val = data[row_idx, imag_column]
                value = Complex(real_val, imag_val)
                push!(values, value)
                break  # Found the coordinate in this bin, move to next bin
            end
        end
    end

    return values
end

"""Print header for correlation analysis results"""
function print_correlation_header(filename, Lx, Ly, n_bins, selected_bins, filtered_bins, dropmaxmin, coord_groups, startbin=1, endbin=nothing)
    println("\nCorrelation analysis of $filename:")
    println("----------------------------------------")
    println("Lattice size: $(Lx)×$(Ly)")

    # Print information about bins
    if n_bins > 1
        # Calculate how many bins were selected before max/min filtering
        selected_count = min(isnothing(endbin) ? n_bins : endbin, n_bins) - startbin + 1

        # Describe all filtering steps
        if startbin > 1 || !isnothing(endbin)
            bin_range = "$startbin-$(isnothing(endbin) ? n_bins : min(endbin, n_bins))"
            println("Bins: $filtered_bins (selected range $bin_range from $n_bins total, then removed $dropmaxmin max/min)")
        else
            println("Bins: $filtered_bins (removed $dropmaxmin max/min from $n_bins total)")
        end
    else
        println("Only one bin found - no filtering applied")
    end

    # Print table header
    println("\n----------------------------------------")
    println(format_correlation_header())
    println("----------------------------------------")
end

"""
    MultiOrbitalCorrelationAnalysis(filename="nn_r.bin", filedir=pwd();
                                   startbin=2, endbin=nothing, dropmaxmin=0,
                                   orbital_columns=[(3,4), (5,6), (7,8), (9,10)],
                                   orbital_labels=["AA", "AB", "BA", "BB"],
                                   auto_digits=true,
                                   verbose=true)

Analyze multi-orbital correlation function data files, where each orbital pair has its own
real and imaginary columns. This is particularly useful for Honeycomb lattice models with
multiple sublattices (e.g., AA, AB, BA, BB orbitals).

Arguments:
- `filename`: Name of the multi-orbital correlation file (default: "nn_r.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 2)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 0)
- `orbital_columns`: Array of (real_column, imag_column) tuples for each orbital (default: [(3,4), (5,6), (7,8), (9,10)])
- `orbital_labels`: Array of labels for each orbital (default: ["AA", "AB", "BA", "BB"])
- `auto_digits`: Whether to automatically determine significant digits based on error of error (default: true)
- `verbose`: Whether to print results to console (default: true)

Returns:
- Dictionary with imj coordinate pairs as keys and a dictionary of orbital results as values,
  where each orbital result contains [mean, error] arrays with complex numbers

Example:
```julia
# Analyze nn_r.bin in current directory with default orbital columns
results = MultiOrbitalCorrelationAnalysis()

# Access results
corr_1_1_AA = results[(1,1)]["AA"][1]  # mean value at imj=(1,1) for AA orbital
corr_1_1_AB = results[(1,1)]["AB"][1]  # mean value at imj=(1,1) for AB orbital
```
"""
function MultiOrbitalCorrelationAnalysis(filename::String="nn_r.bin", filedir::String=pwd();
                                       startbin::Int=2, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=0,
                                       orbital_columns::Vector{Tuple{Int,Int}}=[(3,4), (5,6), (7,8), (9,10)],
                                       orbital_labels::Vector{String}=["AA", "AB", "BA", "BB"],
                                       auto_digits::Bool=true,
                                       verbose::Bool=true)
    # Validate inputs
    if length(orbital_columns) != length(orbital_labels)
        error("Number of orbital columns ($(length(orbital_columns))) must match number of labels ($(length(orbital_labels)))")
    end

    # Add .bin extension if not present
    if !endswith(filename, ".bin")
        filename = filename * ".bin"
    end

    # Open and read the file to check if it exists
    filepath = joinpath(filedir, filename)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    # Read data to infer lattice size and provide diagnostic information
    data = readdlm(filepath, Float64)
    Lx = Int(maximum(data[:, 1]))
    Ly = Int(maximum(data[:, 2]))

    # Print diagnostic information if verbose
    if verbose
        print_diagnostic_info(data, filename)
    end

    # Analyze each orbital separately
    orbital_results = analyze_orbitals(filename, filedir, orbital_columns, orbital_labels,
                                      startbin, endbin, dropmaxmin, auto_digits)

    # Combine results into a single dictionary
    combined_results, sorted_coords = combine_orbital_results(orbital_results, orbital_labels)

    # Print combined results if verbose
    if verbose
        print_multi_orbital_results(filename, Lx, Ly, data, dropmaxmin, orbital_labels,
                                   combined_results, sorted_coords, startbin, endbin)
    end

    return combined_results
end

# Helper functions for MultiOrbitalCorrelationAnalysis

"""Print diagnostic information about the data file"""
function print_diagnostic_info(data, filename)
    println("\nData structure for $filename:")
    println("----------------------------------------")

    # Determine if the file has a bin structure
    unique_coords = Set([(data[i,1], data[i,2]) for i in 1:size(data,1)])
    n_coords = length(unique_coords)
    n_bins = size(data,1) ÷ n_coords

    println("Lattice size: $(Int(maximum(data[:, 1])))×$(Int(maximum(data[:, 2])))")
    println("Total rows: $(size(data, 1))")
    println("Unique coordinates: $n_coords")
    println("Total bins: $n_bins")

    println("----------------------------------------")
end

"""Analyze each orbital separately"""
function analyze_orbitals(filename, filedir, orbital_columns, orbital_labels, startbin, endbin, dropmaxmin, auto_digits)
    orbital_results = Dict{String, Dict{Tuple{Int,Int}, Vector{Complex{Float64}}}}()

    for (i, (real_col, imag_col)) in enumerate(orbital_columns)
        orbital_label = orbital_labels[i]

        # Call CorrelationAnalysis for this orbital
        results = CorrelationAnalysis(filename, filedir;
                                     startbin=startbin, endbin=endbin, dropmaxmin=dropmaxmin,
                                     real_column=real_col, imag_column=imag_col,
                                     auto_digits=auto_digits,
                                     verbose=false)  # Don't print individual orbital results

        orbital_results[orbital_label] = results
    end

    return orbital_results
end

"""Combine results from all orbitals"""
function combine_orbital_results(orbital_results, orbital_labels)
    combined_results = Dict{Tuple{Int,Int}, Dict{String, Vector{Complex{Float64}}}}()

    # Get all unique imj coordinates across all orbitals
    all_coords = Set{Tuple{Int,Int}}()
    for orbital_result in values(orbital_results)
        union!(all_coords, keys(orbital_result))
    end

    # Sort coordinates for ordered output
    sorted_coords = sort(collect(all_coords),
                         by = k -> (euclidean_distance(k[1], k[2]), k[1], k[2]))

    # Combine results for each coordinate
    for coord in sorted_coords
        combined_results[coord] = Dict{String, Vector{Complex{Float64}}}()

        for orbital_label in orbital_labels
            if haskey(orbital_results[orbital_label], coord)
                combined_results[coord][orbital_label] = orbital_results[orbital_label][coord]
            end
        end
    end

    return combined_results, sorted_coords
end

"""Print results for all orbitals"""
function print_multi_orbital_results(filename, Lx, Ly, data, dropmaxmin, orbital_labels, combined_results, sorted_coords, startbin=1, endbin=nothing)
    # Calculate bin information
    nrows = size(data, 1)
    unique_coords = Set([(Int(data[i,1]), Int(data[i,2])) for i in 1:size(data,1)])
    n_coords = length(unique_coords)
    n_bins = nrows ÷ n_coords

    # Calculate how many bins were selected before max/min filtering
    selected_count = min(isnothing(endbin) ? n_bins : endbin, n_bins) - startbin + 1

    # Calculate how many bins remain after all filtering
    filtered_bins = max(1, selected_count - 2 * dropmaxmin)

    # Print header
    println("\nMulti-orbital correlation analysis of $filename:")
    println("----------------------------------------")
    println("Lattice size: $(Lx)×$(Ly)")
    println("Orbitals analyzed: $(join(orbital_labels, ", "))")
    # Describe all filtering steps
    if startbin > 1 || !isnothing(endbin)
        bin_range = "$startbin-$(isnothing(endbin) ? n_bins : min(endbin, n_bins))"
        println("Bins: $filtered_bins (selected range $bin_range from $n_bins total, then removed $dropmaxmin max/min)")
    else
        println("Bins: $filtered_bins (removed $dropmaxmin max/min from $n_bins total)")
    end

    # Print detailed results for each orbital
    print_orbital_details(orbital_labels, combined_results, sorted_coords, Lx, Ly)

    # Print compact summary
    print_compact_summary(orbital_labels, combined_results, sorted_coords)
end

"""Print detailed results for each orbital"""
function print_orbital_details(orbital_labels, combined_results, sorted_coords, Lx, Ly)
    col_width = 25  # Width for value columns

    for orbital_label in orbital_labels
        println("\n----------------------------------------")
        println("Orbital: $(orbital_label)")
        println("----------------------------------------")
        println(format_correlation_header())
        println("----------------------------------------")

        for coord in sorted_coords
            imj_x, imj_y = coord
            distance = euclidean_distance(imj_x, imj_y)
            i_coord = calculate_i_coord(imj_x, imj_y, Lx, Ly)

            if haskey(combined_results[coord], orbital_label)
                mean_val, err_val = combined_results[coord][orbital_label]

                # Format values with appropriate precision
                formatted_real, formatted_real_err = format_value_error(real(mean_val), real(err_val))
                formatted_imag, formatted_imag_err = format_value_error(imag(mean_val), imag(err_val))

                real_data = "$(formatted_real) ± $(formatted_real_err)"
                imag_data = "$(formatted_imag) ± $(formatted_imag_err)"

                # Format and print row
                row = format_correlation_row(coord, distance, i_coord, real_data, imag_data)
                println(row)
            else
                # Format and print row with N/A values
                row = format_correlation_row(coord, distance, i_coord, "N/A", "N/A")
                println(row)
            end
        end
    end
end

"""Print compact summary with real parts only"""
function print_compact_summary(orbital_labels, combined_results, sorted_coords)
    col_width = 25  # Width for each orbital column

    println("\n----------------------------------------")
    println("\nCompact summary (real parts only):")
    println("----------------------------------------")

    # Print header
    header = @sprintf("%-15s", "imj")
    for label in orbital_labels
        header *= " " * lpad(label, col_width)
    end
    println(header)
    println("----------------------------------------")

    # Print each coordinate
    for coord in sorted_coords
        imj_x, imj_y = coord

        # Format coordinate
        coord_str = "($imj_x,$imj_y)"
        row = @sprintf("%-15s", coord_str)

        # Add each orbital's data
        for orbital_label in orbital_labels
            if haskey(combined_results[coord], orbital_label)
                mean_val, err_val = combined_results[coord][orbital_label]
                formatted_real, formatted_real_err = format_value_error(real(mean_val), real(err_val))
                orbital_data = "$(formatted_real) ± $(formatted_real_err)"
                row *= " " * lpad(orbital_data, col_width)
            else
                row *= " " * lpad("N/A", col_width)
            end
        end
        println(row)
    end

    println("----------------------------------------\n")
end
