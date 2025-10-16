#=
Example: Using generic functions to analyze custom k-space correlators
without writing specialized wrapper functions.

This demonstrates:
1. Single-directory analysis with StructureFactorAnalysis
2. Multi-directory batch analysis with analyze_structure_factor_multi_parameter

Key idea: By specifying real_column and imag_column, you can extract any 
orbital pair from multi-orbital k-space files like cpcm_k.bin, nn_k.bin, etc.
=#

using DataProcessforDQMC

# ============================================================================ #
#                       Part 1: Single Directory Analysis                      #
# ============================================================================ #

println("="^80)
println("Part 1: Single Directory Analysis - Direct use of StructureFactorAnalysis")
println("="^80)

# Path to example data directory
example_dir = joinpath(@__DIR__, "proj_bt_honeycomb_exact_b=36.000_U=4.00_L=24_dtau=0.10_xmag=0.0")

# Example 1.1: Analyze cpcm_k.bin, AA orbital (columns 3, 4)
# File format: kx ky [AA_real AA_imag] [AB_real AB_imag] [BA_real BA_imag] [BB_real BB_imag]
println("\n--- Example 1.1: cpcm_k.bin at Γ point, AA orbital ---")
result_AA = StructureFactorAnalysis(
    (0.0, 0.0),           # k-point: Γ point
    "cpcm_k.bin",         # filename
    example_dir;          # directory
    real_column=3,        # AA orbital real part (column 3)
    imag_column=4,        # AA orbital imag part (column 4)
    startbin=2,           # skip first bin
    dropmaxmin=0,         # don't drop outliers
    verbose=true          # print detailed results
)

println("Result: G_AA(Γ) = $(result_AA.formatted_real)")

# Example 1.2: Analyze cpcm_k.bin, AB orbital (columns 5, 6) at different k-point
println("\n--- Example 1.2: cpcm_k.bin at (π,π) point, AB orbital ---")
result_AB = StructureFactorAnalysis(
    (0.5, 0.5),               # k-point: (π,π) point
    "cpcm_k.bin",
    example_dir;
    real_column=5,        # AB orbital real part (column 5)
    imag_column=6,        # AB orbital imag part (column 6)
    startbin=2,
    dropmaxmin=0,
    tolerance=1e-6,       # tolerance for k-point matching
    verbose=true
)

println("Result: G_AB(π,π) = $(result_AB.formatted_real)")

# Example 1.3: Analyze nn_k.bin (density-density correlator), BB orbital
println("\n--- Example 1.3: nn_k.bin at Γ point, BB orbital ---")
result_nn_BB = StructureFactorAnalysis(
    (0.0, 0.0),
    "nn_k.bin",
    example_dir;
    real_column=9,        # BB orbital real part (column 9)
    imag_column=10,       # BB orbital imag part (column 10)
    startbin=2,
    dropmaxmin=0,
    verbose=true
)

println("Result: N_BB(Γ) = $(result_nn_BB.formatted_real)")

# ============================================================================ #
#                    Part 2: Multi-Directory Batch Analysis                    #
# ============================================================================ #

println("\n" * "="^80)
println("Part 2: Multi-Directory Batch Analysis")
println("="^80)

# For demonstration, we'll use the parent directory
# In practice, you would have multiple parameter directories
base_dir = dirname(example_dir)

# Example 2.1: Batch analysis of cpcm_k.bin, AB orbital at Γ point
println("\n--- Example 2.1: Batch analysis of cpcm_k.bin, AB orbital ---")

df_AB = analyze_structure_factor_multi_parameter(
    # Anonymous function that wraps StructureFactorAnalysis
    # This function receives (k_point, filename, dir; kwargs...) from the driver
    # Note: filter out force_rebuild and source_file as they're not used by StructureFactorAnalysis
    (k_point, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(
            k_point, 
            filename, 
            dir;
            real_column=5,     # Fix: AB orbital real part
            imag_column=6,     # Fix: AB orbital imag part
            kwargs...          # Forward relevant parameters (startbin, endbin, etc.)
        ),
    base_dir;              # Base directory to scan
    k_point=(0.0, 0.0),    # k-point to analyze
    filename="cpcm_k.bin", # Target file
    result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],  # Match StructureFactorAnalysis output
    result_prefix="G_AB",  # Prefix for result columns in DataFrame
    startbin=2,
    dropmaxmin=0,
    verbose=true,
    filter_options=(prefix="proj_bt_honeycomb_exact",)  # Use complete prefix for exact match
)

println("\n--- Results DataFrame ---")
println(df_AB)

# Example 2.2: Batch analysis of nn_k.bin, AA orbital at (π,π)
println("\n--- Example 2.2: Batch analysis of nn_k.bin, AA orbital at (π,π) ---")

df_nn_AA = analyze_structure_factor_multi_parameter(
    (k_point, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(
            k_point, 
            filename, 
            dir;
            real_column=3,     # Fix: AA orbital real part
            imag_column=4,     # Fix: AA orbital imag part
            kwargs...
        ),
    base_dir;
    k_point=(0.5, 0.5),
    filename="nn_k.bin",
    result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
    result_prefix="N_AA",
    startbin=2,
    dropmaxmin=0,
    verbose=true,
    filter_options=(prefix="proj_bt_honeycomb_exact", U=4.0)  # Can filter by multiple parameters
)

println("\n--- Results DataFrame ---")
println(df_nn_AA)