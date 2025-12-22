#=
Example: Multi-k-point structure factor analysis

This demonstrates the new multi-k-point analysis functionality,
which analyzes multiple k-points in a single scan of the data file.
=#

using DataProcessforDQMC
using DataFrames

println("="^80)
println("Multi-k-point Structure Factor Analysis Test")
println("="^80)

# Path to example data directory
base_dir = joinpath(@__DIR__, "proj_bt_honeycomb_exact_b=36.000_U=4.00_L=24_dtau=0.10_xmag=0.0")

# ============================================================================ #
#                  Part 1: Single Directory Multi-k Analysis                   #
# ============================================================================ #

println("\n--- Part 1: Single Directory Multi-k Analysis ---")

# Define multiple k-points to analyze
k_points = [
    (0.0, 0.0),    # Γ point
    (0.5, 0.5),    # (π,π) point
    (0.5, 0.0),    # (π,0) point
    (0.0, 0.5),    # (0,π) point
]

println("Analyzing cpcm_k.bin at $(length(k_points)) k-points:")
for (i, k) in enumerate(k_points)
    println("  $i. k = $k")
end

# Call the multi-k version of StructureFactorAnalysis
results = StructureFactorAnalysis(
    k_points,             # Vector of k-points
    "cpcm_k.bin",
    base_dir;
    real_column=5,        # AB orbital real part
    imag_column=6,        # AB orbital imag part
    startbin=2,
    outlier_mode=:dropmaxmin,
    outlier_param=0,
    k_point_tolerance=1e-6,
    verbose=true
)

println("\n--- Multi-k Analysis Results ---")
println("Returned $(length(results)) results:")
for (i, result) in enumerate(results)
    println("  $i. k=$(result.k_point): $(result.formatted_real)")
end

# ============================================================================ #
#              Part 2: Multi-Directory Multi-k Batch Analysis                  #
# ============================================================================ #

println("\n" * "="^80)
println("--- Part 2: Multi-Directory Multi-k Batch Analysis ---")
println("="^80)

# Scan parent directory for multiple parameter sets
scan_base_dir = dirname(base_dir)

# Use the multi-k version of analyze_structure_factor_multi_parameter
df_multi_k = analyze_structure_factor_multi_parameter(
    (k_points, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(
            k_points,          # Pass multiple k-points
            filename, 
            dir;
            real_column=5,     # AB orbital real part
            imag_column=6,     # AB orbital imag part
            kwargs...
        ),
    k_points,              # k-points to analyze (now a positional argument)
    scan_base_dir;         # Base directory to scan
    filename="cpcm_k.bin",
    result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
    result_prefix="G_AB",
    startbin=2,
    outlier_mode=:dropmaxmin,
    outlier_param=0,
    verbose=true,
    filter_options=(prefix="proj_bt_honeycomb_exact",)
)

println("\n--- Multi-k Results DataFrame (Long Table Format) ---")
println("Total rows: $(nrow(df_multi_k)) (parameters × k-points)")
println("\nDataFrame preview:")
println(df_multi_k)

println("\n--- Verifying Long Table Structure ---")
println("Columns: ", names(df_multi_k))
println("Has :kx column: ", :kx in names(df_multi_k))
println("Has :ky column: ", :ky in names(df_multi_k))
println("Unique k-points: ", unique([(row.kx, row.ky) for row in eachrow(df_multi_k)]))

println("\n✅ Multi-k-point analysis test completed successfully!")

