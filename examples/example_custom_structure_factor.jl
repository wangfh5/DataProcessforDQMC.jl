#=
Example: Using generic functions to analyze custom correlators in both 
k-space and real-space without writing specialized wrapper functions.

This demonstrates:
1. Single-directory analysis with StructureFactorAnalysis (k-space and r-space)
2. Multi-directory batch analysis with analyze_structure_factor_multi_parameter

Key idea: By specifying real_column and imag_column, you can extract any 
orbital pair from multi-orbital files:
- K-space files: cpcm_k.bin, nn_k.bin, etc. (columns: kx, ky, ...)
- R-space files: cpcm_r.bin, nn_r.bin, etc. (columns: imj_x, imj_y, ...)

The framework treats coordinates generically - whether they are k-points or 
r-points doesn't matter to the underlying analysis machinery.
=#

using DataProcessforDQMC

# ============================================================================ #
#                       Part 1: Single Directory Analysis                      #
# ============================================================================ #

println("="^80)
println("Part 1: Single Directory Analysis (K-space and R-space)")
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
    k_point_tolerance=1e-6,  # tolerance for k-point matching
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

# Example 1.4: Analyze cpcm_r.bin (real-space Green's function), AB orbital
# Note: The framework treats (imj_x, imj_y) coordinates just like (kx, ky)
println("\n--- Example 1.4: cpcm_r.bin at r=(1,1) point, AB orbital ---")
result_r_AB = StructureFactorAnalysis(
    (1, 1),               # r-point: (imj_x, imj_y) = (1, 1) in unit cell coordinates
    "cpcm_r.bin",         # real-space file
    example_dir;
    real_column=7,        # AB orbital real part (column 7 in cpcm_r.bin)
    imag_column=8,        # AB orbital imag part (column 8)
    startbin=2,
    dropmaxmin=0,
    k_point_tolerance=1e-6,  # For r-space integer coordinates, this ensures exact match
    verbose=true
)

println("Result: G_AB(r=(1,1)) = $(result_r_AB.formatted_real)")

# ============================================================================ #
#                    Part 2: Multi-Directory Batch Analysis                    #
# ============================================================================ #

println("\n" * "="^80)
println("Part 2: Multi-Directory Batch Analysis (K-space and R-space)")
println("="^80)

# For demonstration, we'll use the parent directory
# In practice, you would have multiple parameter directories
base_dir = dirname(example_dir)

# Example 2.1: Batch analysis of cpcm_k.bin, AB orbital at Γ point
println("\n--- Example 2.1: Batch analysis of cpcm_k.bin, AB orbital ---")

df_AB = analyze_structure_factor_multi_parameter(
    # Anonymous function that wraps StructureFactorAnalysis
    # Note: filter out force_rebuild and source_file (not used by StructureFactorAnalysis)
    (k_point, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(
            k_point, 
            filename, 
            dir;
            real_column=7,     # Fix: AB orbital real part
            imag_column=8,     # Fix: AB orbital imag part
            kwargs...          # Auto-inherit: startbin, dropmaxmin, etc.
        ),
    (0.0, 0.0),            # Γ point
    base_dir;              # Base directory to scan
    filename="cpcm_k.bin",
    # Note: DO NOT set force_rebuild here to allow file existence check
    result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
    result_prefix="G_AB",
    startbin=2,
    dropmaxmin=0,
    verbose=true,
    filter_options=(prefix="proj_bt_honeycomb_exact", U=4.0)
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
            kwargs...          # Auto-inherit: startbin, dropmaxmin, etc.
        ),
    (0.5, 0.5),            # (π,π) point
    base_dir;
    filename="nn_k.bin",
    # Note: DO NOT set force_rebuild here to allow file existence check
    result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
    result_prefix="N_AA",
    startbin=2,
    dropmaxmin=0,
    verbose=true,
    filter_options=(prefix="proj_bt_honeycomb_exact", U=4.0)
)

println("\n--- Results DataFrame ---")
println(df_nn_AA)

# Example 2.3: Batch analysis of cpcm_r.bin (real-space), AB orbital at specific r-point
# This demonstrates that the same multi-parameter framework works seamlessly for r-space
println("\n--- Example 2.3: Batch analysis of cpcm_r.bin, AB orbital at r=(12,12) ---")

df_r_AB = analyze_structure_factor_multi_parameter(
    # Anonymous function wrapper - identical pattern as k-space examples
    # The framework doesn't care whether coordinates are k-points or r-points
    (r_point, filename, dir; force_rebuild=false, source_file="", kwargs...) -> 
        StructureFactorAnalysis(
            r_point, 
            filename, 
            dir;
            real_column=7,     # Fix: AB orbital real part in cpcm_r.bin
            imag_column=8,     # Fix: AB orbital imag part in cpcm_r.bin
            kwargs...          # Auto-inherit: startbin, dropmaxmin, etc.
        ),
    (12, 12),              # r-point to analyze (unit cell coordinates)
    base_dir;              # Base directory to scan
    filename="cpcm_r.bin", # Real-space file
    # Note: DO NOT set force_rebuild here to allow file existence check
    result_columns=[:mean_real, :err_real, :mean_imag, :err_imag],
    result_prefix="G_AB_r",
    startbin=2,
    dropmaxmin=0,
    verbose=true,
    filter_options=(prefix="proj_bt_honeycomb_exact", L=24)
)

println("\n--- Results DataFrame ---")
println(df_r_AB)