#=
Bin File Export Functions (bin-file-export.jl)

This file provides simplified functions to export bin files to tabular formats
(CSV and JLD2) without metadata or statistical analysis. It's designed for
quick data export and visualization preparation.

Main functions:
- export_bin_to_dataframe: Convert bin file to DataFrame
- export_bin_to_csv: Export bin file to CSV
- export_bin_to_jld2: Export bin file to JLD2
- export_directory_bins: Batch export all bin files in a directory
=#

export export_bin_to_dataframe, export_bin_to_csv, export_bin_to_jld2
export export_directory_bins

const DEFAULT_ORBITAL_LABELS = ["AA", "AB", "BA", "BB"]

"""
    export_bin_to_dataframe(
        input_file::String,
        dir::String=pwd();
        iscorrelation::Bool=true,
        orbital_columns::Union{Vector{Tuple{Int,Int}}, Nothing}=nothing,
        orbital_labels::Union{Vector{String}, Nothing}=nothing
    )

Convert a bin file to a DataFrame with bin index column.

This function reads a bin file and converts it to a DataFrame format suitable
for analysis and visualization. It automatically detects the file structure
(k-space vs r-space) and handles both multi-orbital and single-column data.

# Arguments
- `input_file::String`: Input .bin file name
- `dir::String=pwd()`: Directory containing the bin file
- `iscorrelation::Bool=true`: Whether to interpret the file as a correlation function (coordinates + real/imag pairs)
- `orbital_columns::Union{Vector{Tuple{Int,Int}}, Nothing}=nothing`: Column indices for each orbital pair (real, imag). Auto-generated when `nothing`
- `orbital_labels::Union{Vector{String}, Nothing}=nothing`: Labels for each orbital pair (defaults to AA/AB/... or pair1/pair2/...)

# Returns
- `DataFrame`: Converted data with columns: bin, coord1, coord2, data columns...

# Data Structure
The function handles multiple formats:
1. **Multi-orbital** (default 10 columns): kx, ky, AA_real, AA_imag, AB_real, AB_imag, BA_real, BA_imag, BB_real, BB_imag
   - Note: AA, AB, BA, BB represent orbital pairs in two-point correlation functions <c†_{i,α} c_{j,β}>
2. **Custom orbital pairs**: Specify any number of orbital pairs via `orbital_columns` and `orbital_labels`
3. **Single column** (4 columns): kx, ky, value_real, value_imag

# Examples
```julia
# Export k-space correlation function (default 4 orbital pairs: AA, AB, BA, BB)
df = export_bin_to_dataframe("nn_k.bin", dir="data/")

# Export only diagonal orbital pairs (AA and BB)
df = export_bin_to_dataframe("nn_k.bin", 
    orbital_columns=[(3,4), (9,10)],
    orbital_labels=["AA", "BB"])

# Export r-space correlation function
df = export_bin_to_dataframe("nn_r.bin", dir="data/")

# Export AFM structure factor (single column)
df = export_bin_to_dataframe("afm_sf_k.bin", dir="data/", iscorrelation=true, orbital_columns=[(3,4)], orbital_labels=["value"])
```
"""
function export_bin_to_dataframe(
    input_file::String,
    dir::String=pwd();
    iscorrelation::Bool=true,
    orbital_columns::Union{Vector{Tuple{Int,Int}}, Nothing}=nothing,
    orbital_labels::Union{Vector{String}, Nothing}=nothing
)
    # Read the bin file
    input_path = joinpath(dir, input_file)
    if !isfile(input_path)
        @error "Input file not found: $input_path"
        return DataFrame()
    end
    
    data = readdlm(input_path)
    
    n_rows, n_cols = size(data, 1), size(data, 2)

    # Handle scalar (bin-indexed) data: two columns => value_real/value_imag
    if n_cols == 2
        return DataFrame(
            :bin => 1:n_rows,
            :value_real => data[:, 1],
            :value_imag => data[:, 2]
        )
    end

    if iscorrelation && n_cols >= 4
        coords = unique(data[:, 1:2], dims=1)
        n_coords = size(coords, 1)
        n_bins = n_rows ÷ n_coords
        bin_indices = repeat(1:n_bins, inner=n_coords)

        # Auto-detect coordinate type (k-space vs r-space)
        sample_coords = coords[1:min(10, n_coords), :]
        is_k_space = !all(x -> isapprox(x, round(x), atol=1e-6), sample_coords)
        coord_names = is_k_space ? [:kx, :ky] : [:rx, :ry]

        col_count = n_cols - 2

        if col_count == 2
            # Single orbital correlation (one real/imag pair)
            return DataFrame(
                :bin => bin_indices,
                coord_names[1] => data[:, 1],
                coord_names[2] => data[:, 2],
                :value_real => data[:, 3],
                :value_imag => data[:, 4]
            )
        elseif col_count % 2 == 0
            n_pairs = col_count ÷ 2
            cols = isnothing(orbital_columns) ? Tuple{Int,Int}[(2 + (i-1)*2 + 1, 2 + (i-1)*2 + 2) for i in 1:n_pairs] : orbital_columns
            labels = if isnothing(orbital_labels)
                n_pairs <= length(DEFAULT_ORBITAL_LABELS) ? DEFAULT_ORBITAL_LABELS[1:n_pairs] : ["pair$(i)" for i in 1:n_pairs]
            else
                orbital_labels
            end

            if length(cols) != length(labels)
                @error "Number of orbital_columns ($(length(cols))) must match orbital_labels ($(length(labels)))"
                return DataFrame()
            end

            df = DataFrame(
                :bin => bin_indices,
                coord_names[1] => data[:, 1],
                coord_names[2] => data[:, 2]
            )

            valid_pairs = false
            for (pair_idx, (real_col, imag_col)) in enumerate(cols)
                orbital = labels[pair_idx]
                if real_col <= n_cols && imag_col <= n_cols
                    df[!, Symbol("$(orbital)_real")] = data[:, real_col]
                    df[!, Symbol("$(orbital)_imag")] = data[:, imag_col]
                    valid_pairs = true
                else
                    @warn "Orbital $orbital: columns ($real_col, $imag_col) exceed file columns ($n_cols), skipping"
                end
            end

            if valid_pairs
                return df
            else
                @warn "No valid orbital pairs detected in $input_file; falling back to generic DataFrame" input_file
            end
        else
            @warn "Unable to interpret correlation columns for $input_file (column count = $n_cols); falling back to generic DataFrame" input_file n_cols
        end
    elseif iscorrelation && n_cols < 4
        @warn "File $input_file has fewer than four columns; falling back to generic DataFrame" input_file n_cols
    end

    # Generic fallback: treat each column as col_i and provide sequential bin index
    col_names = [Symbol("col_$i") for i in 1:n_cols]
    df = DataFrame(data, col_names)
    insertcols!(df, 1, :bin => 1:n_rows)
    return df
end

"""
    export_bin_to_csv(
        input_file::String,
        dir::String=pwd();
        output_dir::Union{String, Nothing}=nothing,
        output_file::Union{String, Nothing}=nothing,
        kwargs...
    )

Export a bin file to CSV format.

# Arguments
- `input_file::String`: Input .bin file name
- `dir::String=pwd()`: Directory containing the bin file
- `output_dir::Union{String, Nothing}=nothing`: Output directory (defaults to dir/exported_csv)
- `output_file::Union{String, Nothing}=nothing`: Output filename (defaults to input_file with .csv extension)
- `kwargs...`: Additional arguments passed to `export_bin_to_dataframe`

# Returns
- `String`: Path to the exported CSV file

# Examples
```julia
# Export to default location (data/exported_csv/)
csv_file = export_bin_to_csv("nn_k.bin", dir="data/")

# Export to custom location
csv_file = export_bin_to_csv("nn_k.bin", dir="data/", output_dir="results/")

# Export with custom filename
csv_file = export_bin_to_csv("nn_k.bin", output_file="correlation_nn.csv")
```
"""
function export_bin_to_csv(
    input_file::String,
    dir::String=pwd();
    output_dir::Union{String, Nothing}=nothing,
    output_file::Union{String, Nothing}=nothing,
    kwargs...
)
    # Convert to DataFrame
    df = export_bin_to_dataframe(input_file, dir; kwargs...)
    
    if isempty(df)
        @error "Failed to convert bin file to DataFrame"
        return ""
    end
    
    # Determine output directory
    if isnothing(output_dir)
        output_dir = joinpath(dir, "exported_csv")
    end
    mkpath(output_dir)
    
    # Determine output filename
    if isnothing(output_file)
        output_file = replace(input_file, ".bin" => ".csv")
    end
    
    # Ensure .csv extension
    if !endswith(output_file, ".csv")
        output_file = output_file * ".csv"
    end
    
    output_path = joinpath(output_dir, output_file)
    
    # Write CSV
    CSV.write(output_path, df)
    
    println("✓ Exported to CSV: $output_path")
    return output_path
end

"""
    export_bin_to_jld2(
        input_file::String,
        dir::String=pwd();
        output_dir::Union{String, Nothing}=nothing,
        output_file::Union{String, Nothing}=nothing,
        dataset_name::Union{String, Nothing}=nothing,
        kwargs...
    )

Export a bin file to JLD2 format.

# Arguments
- `input_file::String`: Input .bin file name
- `dir::String=pwd()`: Directory containing the bin file
- `output_dir::Union{String, Nothing}=nothing`: Output directory (defaults to dir/exported_jld2)
- `output_file::Union{String, Nothing}=nothing`: Output filename (defaults to input filename with .jld2 extension, e.g., nn_k.jld2)
- `dataset_name::Union{String, Nothing}=nothing`: Dataset name in JLD2 file (defaults to base filename, e.g., "nn_k")
- `kwargs...`: Additional arguments passed to `export_bin_to_dataframe`

# Returns
- `String`: Path to the exported JLD2 file

# Examples
```julia
# Export to default JLD2 file (nn_k.jld2 with dataset "nn_k")
jld2_file = export_bin_to_jld2("nn_k.bin", dir="data/")
# Creates: data/exported_jld2/nn_k.jld2 (dataset: "nn_k")
# Read: load("nn_k.jld2", "nn_k")

# Export to custom JLD2 file
jld2_file = export_bin_to_jld2("nn_k.bin", output_file="correlations.jld2")
# Creates: correlations.jld2 (dataset: "nn_k")

# Multiple datasets in one file (if needed)
export_bin_to_jld2("nn_k.bin", output_file="all_data.jld2")  # dataset: "nn_k"
export_bin_to_jld2("nn_r.bin", output_file="all_data.jld2")  # dataset: "nn_r"
# Read: load("all_data.jld2", "nn_k") or load("all_data.jld2", "nn_r")
```
"""
function export_bin_to_jld2(
    input_file::String,
    dir::String=pwd();
    output_dir::Union{String, Nothing}=nothing,
    output_file::Union{String, Nothing}=nothing,
    dataset_name::Union{String, Nothing}=nothing,
    kwargs...
)
    # Convert to DataFrame
    df = export_bin_to_dataframe(input_file, dir; kwargs...)
    
    if isempty(df)
        @error "Failed to convert bin file to DataFrame"
        return ""
    end
    
    # Determine output directory
    if isnothing(output_dir)
        output_dir = joinpath(dir, "exported_jld2")
    end
    mkpath(output_dir)
    
    # Determine base name from input file
    base_name = replace(input_file, ".bin" => "")
    
    # Determine output filename (default: same as input file with .jld2 extension)
    if isnothing(output_file)
        output_file = base_name * ".jld2"
    end
    
    # Ensure .jld2 extension
    if !endswith(output_file, ".jld2")
        output_file = output_file * ".jld2"
    end
    
    output_path = joinpath(output_dir, output_file)
    
    # Determine dataset name (default: same as base name)
    if isnothing(dataset_name)
        dataset_name = base_name
    end
    
    # Write to JLD2 file
    if isfile(output_path)
        # File exists: append or overwrite
        jldopen(output_path, "r+") do file
            if haskey(file, dataset_name)
                @warn "Dataset '$dataset_name' already exists in $output_path, overwriting..."
                delete!(file, dataset_name)
            end
            file[dataset_name] = df
        end
        println("✓ Updated JLD2 file '$output_file' with dataset '$dataset_name'")
    else
        # Create new file
        jldopen(output_path, "w") do file
            file[dataset_name] = df
        end
        println("✓ Exported to JLD2: $output_path")
    end
    
    return output_path
end

"""
    export_directory_bins(
        dir::String=pwd();
        file_patterns::Vector{String}=["*_k.bin", "*_r.bin"],
        output_format::Symbol=:both,
        output_dir::Union{String, Nothing}=nothing,
        verbose::Bool=true,
        iscorrelation::Union{Bool,Nothing}=nothing,
        orbital_columns::Union{Vector{Tuple{Int,Int}}, Nothing}=nothing,
        orbital_labels::Union{Vector{String}, Nothing}=nothing
    )

Batch export all matching bin files in a directory.

# Arguments
- `dir::String=pwd()`: Directory containing bin files
- `file_patterns::Vector{String}=["*_k.bin", "*_r.bin"]`: Glob patterns for files to export
- `output_format::Symbol=:both`: Output format (:csv, :jld2, or :both) - default is :both
- `output_dir::Union{String, Nothing}=nothing`: Output directory
- `verbose::Bool=true`: Print progress information
- `iscorrelation::Union{Bool,Nothing}=nothing`: Whether files should be treated as correlation data. `nothing` enables auto-detection based on filename (contains `_k` or `_r`).
- `orbital_columns::Union{Vector{Tuple{Int,Int}}, Nothing}=nothing`: Optional orbital column mapping passed to export functions (used when `iscorrelation` is `true`).
- `orbital_labels::Union{Vector{String}, Nothing}=nothing`: Optional orbital labels passed to export functions (used when `iscorrelation` is `true`).

# Returns
- `Vector{String}`: Paths to exported files

# Examples
```julia
# Export all k-space and r-space correlation files to both CSV and JLD2 (default)
export_directory_bins("data/")

# Export only to CSV
export_directory_bins("data/", output_format=:csv)

# Export only k-space files
export_directory_bins("data/", file_patterns=["*_k.bin"])

# Export with custom output directory
export_directory_bins("data/", output_dir="results/")
```
"""
function export_directory_bins(
    dir::String=pwd();
    file_patterns::Vector{String}=["*_k.bin", "*_r.bin"],
    output_format::Symbol=:both,
    output_dir::Union{String, Nothing}=nothing,
    verbose::Bool=true,
    iscorrelation::Union{Bool, Nothing}=nothing,
    orbital_columns::Union{Vector{Tuple{Int,Int}}, Nothing}=nothing,
    orbital_labels::Union{Vector{String}, Nothing}=nothing
)
    @assert output_format in [:csv, :jld2, :both] "output_format must be :csv, :jld2, or :both"
    
    exported_files = String[]
    
    verbose && println("=" ^ 70)
    verbose && println("Batch Bin File Export")
    verbose && println("=" ^ 70)
    verbose && println("Directory: $dir")
    verbose && println("Patterns: $(join(file_patterns, ", "))")
    verbose && println("Format: $output_format")
    verbose && println()
    
    # Find all matching files
    all_files = readdir(dir)
    matched_files = String[]
    
    for pattern in file_patterns
        # Simple glob matching (convert * to regex)
        regex_pattern = replace(pattern, "*" => ".*")
        regex_pattern = "^" * regex_pattern * "\$"
        
        for file in all_files
            if occursin(Regex(regex_pattern), file) && file ∉ matched_files
                push!(matched_files, file)
            end
        end
    end
    
    if isempty(matched_files)
        verbose && println("⚠  No files found matching patterns")
        return exported_files
    end
    
    sort!(matched_files)
    verbose && println("Found $(length(matched_files)) files to export")
    verbose && println()
    
    # Export each file
    for (i, file) in enumerate(matched_files)
        verbose && println("[$i/$(length(matched_files))] Processing: $file")
        
        try
            auto_correlation = occursin("_k", file) || occursin("_r", file)
            use_correlation = isnothing(iscorrelation) ? auto_correlation : iscorrelation
            base_kwargs = (; iscorrelation=use_correlation)
            csv_kwargs = base_kwargs
            jld2_kwargs = csv_kwargs
            if use_correlation
                if orbital_columns !== nothing
                    csv_kwargs = merge(csv_kwargs, (; orbital_columns=orbital_columns))
                    jld2_kwargs = merge(jld2_kwargs, (; orbital_columns=orbital_columns))
                end
                if orbital_labels !== nothing
                    csv_kwargs = merge(csv_kwargs, (; orbital_labels=orbital_labels))
                    jld2_kwargs = merge(jld2_kwargs, (; orbital_labels=orbital_labels))
                end
            end
            if output_format in [:csv, :both]
                csv_path = export_bin_to_csv(
                    file,
                    dir,
                    output_dir=output_dir;
                    csv_kwargs...
                )
                push!(exported_files, csv_path)
            end
            
            if output_format in [:jld2, :both]
                jld2_path = export_bin_to_jld2(
                    file,
                    dir,
                    output_dir=output_dir;
                    jld2_kwargs...
                )
                push!(exported_files, jld2_path)
            end
            
            verbose && println()
        catch e
            @warn "Failed to export $file: $e"
            verbose && println()
        end
    end
    
    verbose && println("=" ^ 70)
    verbose && println("✨ Batch export completed!")
    verbose && println("Exported $(length(exported_files)) files")
    verbose && println("=" ^ 70)
    
    return exported_files
end

