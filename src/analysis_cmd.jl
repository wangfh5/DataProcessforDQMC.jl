## These function accept the filename of the data file and then do analysis
export RenyiNegativity, RenyiNegativity_all, EnergyAnalysis, CorrelationAnalysis

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

"""
    EnergyAnalysis(filename="energy.bin", filedir=pwd();
                   startbin=1, endbin=nothing, dropmaxmin=1,
                   columns=[5,7,9,11], labels=["E_kin", "E_pot", "E_tot", "E_tot^2"],
                   verbose=true)

Analyze energy.bin file and calculate Monte Carlo averages and errors for specified columns.

Arguments:
- `filename`: Name of the energy file (default: "energy.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 1)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 1)
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
                       startbin::Int=1, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=1,
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

# Distance functions for correlation analysis

"""
    CorrelationAnalysis(filename="spsm_r.bin", filedir=pwd();
                        startbin=1, endbin=nothing, dropmaxmin=1,
                        real_column=3, imag_column=4,
                        auto_digits=true,
                        verbose=true)

Analyze correlation function data files like `spsm_r.bin` or `nn_r.bin`, where the first two columns
represent imj coordinates (distance vector between points i and j), and calculate Monte Carlo averages.
#TODO for Honeycomb: Add support for honeycomb lattice with multiple sublattices

Arguments:
- `filename`: Name of the correlation file (default: "spsm_r.bin")
- `filedir`: Directory containing the file (default: current directory)
- `startbin`: First bin to include in analysis (default: 1)
- `endbin`: Last bin to include in analysis (default: all bins)
- `dropmaxmin`: Number of maximum and minimum values to drop (default: 1)
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
                            startbin::Int=1, endbin::Union{Int,Nothing}=nothing, dropmaxmin::Int=1,
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

    # We will apply bin selection for each imj group separately

    # Infer lattice size from data - maximum values of first and second columns
    Lx = Int(maximum(data[:, 1]))
    Ly = Int(maximum(data[:, 2]))

    # Group data by imj coordinates
    coord_groups = Dict{Tuple{Int,Int}, Vector{Complex{Float64}}}()

    # Process each row
    for i in 1:size(data, 1)
        imj_x = Int(data[i, 1])
        imj_y = Int(data[i, 2])
        real_val = data[i, real_column]
        imag_val = data[i, imag_column]
        value = Complex(real_val, imag_val)

        # Use imj coordinate pair as key
        key = (imj_x, imj_y)

        # Add value to the appropriate group
        if haskey(coord_groups, key)
            push!(coord_groups[key], value)
        else
            coord_groups[key] = [value]
        end
    end

    # Apply bin selection and filtering to each group
    for key in keys(coord_groups)
        values = coord_groups[key]

        # Apply startbin and endbin selection
        start_idx = startbin
        end_idx = isnothing(endbin) ? length(values) : min(endbin, length(values))

        if start_idx <= end_idx && start_idx <= length(values)
            values = values[start_idx:end_idx]

            # Apply filtering if needed
            if dropmaxmin > 0
                values = filter(values, dropmaxmin)
            end

            coord_groups[key] = values
        else
            # Empty group if no valid data after bin selection
            coord_groups[key] = Complex{Float64}[]
        end
    end

    # Calculate statistics for each group
    results = Dict{Tuple{Int,Int}, Vector{Complex{Float64}}}()

    # Prepare output table
    if verbose
        println("\nCorrelation analysis of $filename:")
        println("----------------------------------------")
        println("Lattice size: $(Lx)×$(Ly)")
        println("Total measurements after filtering: $(length(coord_groups[(1,1)]))")
        println("----------------------------------------")
        println("imj\t\tdistance\t\ti [to j=(1,1)]\t\treal part\t\timag part")
        println("----------------------------------------")
    end

    # Sort keys for ordered output (first by distance, then by coordinates)
    function euclidean_distance(x, y)
        return sqrt(x^2 + y^2)
    end

    sorted_keys = sort(collect(keys(coord_groups)),
                      by = k -> (euclidean_distance(k[1], k[2]), k[1], k[2]))

    for key in sorted_keys
        imj_x, imj_y = key
        values = coord_groups[key]

        # Extract real and imaginary parts
        real_values = real.(values)
        imag_values = imag.(values)

        # Calculate mean for real and imaginary parts
        mean_real = mean(real_values)
        mean_imag = mean(imag_values)

        # Calculate error for real and imaginary parts
        err_real = error(real_values, sigma=1, bessel=true, auto_digits=true)
        err_imag = error(imag_values, sigma=1, bessel=true, auto_digits=true)

        # Store results as complex numbers
        results[key] = [Complex(mean_real, mean_imag), Complex(err_real, err_imag)]

        # Calculate distance
        distance = euclidean_distance(imj_x, imj_y)

        # Calculate example i coordinate (when j is fixed at (1,1))
        # i = j + imj with periodic boundary conditions
        i_x = mod1(1 + imj_x, Lx)  # mod1 ensures result is in range 1:Lx
        i_y = mod1(1 + imj_y, Ly)  # mod1 ensures result is in range 1:Ly

        # Format values with appropriate precision
        formatted_real, formatted_real_err = format_value_error(mean_real, err_real)
        formatted_imag, formatted_imag_err = format_value_error(mean_imag, err_imag)

        if verbose
            println("($imj_x,$imj_y)\t\t$(round(distance, digits=3))\t\t\t($i_x,$i_y)\t\t\t$formatted_real ± $formatted_real_err\t\t$formatted_imag ± $formatted_imag_err")
        end
    end

    if verbose
        println("----------------------------------------\n")
    end

    return results
end

