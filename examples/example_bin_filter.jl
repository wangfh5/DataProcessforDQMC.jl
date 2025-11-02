#!/usr/bin/env julia

"""
Example script demonstrating the bin file filtering functionality.

This script shows how to extract bin data for specific coordinates from
correlation function files, which is useful for bin convergence analysis.
"""

using DataProcessforDQMC

# Set the data directory
data_dir = joinpath(@__DIR__, "proj_bt_honeycomb_exact_b=36.000_U=4.00_L=24_dtau=0.10_xmag=0.0")


# Example 1: Filter k-space data for a specific momentum point
begin
    println("\n📍 Example 1: K-space filtering")
    println("-" ^ 70)
    println("Extracting all bins for k = (-0.458, -0.458) from nn_k.bin...")

    output_k = filter_bin_file(
        "nn_k.bin",
        (-0.45833333, -0.45833333),
        dir=data_dir,
        verbose=true
    )

    println("\n✓ Output saved to: $(basename(output_k))")
end

# Example 2: Filter r-space data for a specific position
begin
    println("\n\n📍 Example 2: R-space filtering")
    println("-" ^ 70)
    println("Extracting all bins for r = (1, 1) from nn_r.bin...")

    output_r = filter_bin_file(
        "nn_r.bin",
        (1, 1),
        dir=data_dir,
        verbose=true
    )

    println("\n✓ Output saved to: $(basename(output_r))")
end

# Example 3: Generate derived quantity and then filter
begin
    println("\n\n📍 Example 3: Generate AFM structure factor and filter")
    println("-" ^ 70)
    println("Demonstrating the workflow: derive → filter → analyze")
    println()
    
    println("Step 1: Generating afm_sf_k.bin from spsm_k.bin...")
    println("  (This computes S_AFM = AA + BB - AB - BA from XY spin component)")
    
    # Generate AFM structure factor directly from spsm_k.bin (XY component only)
    afm_sf_path = merge_afm_sf(
        "spsm_k.bin",  # Use spsm_k.bin as source instead of ss_k.bin
        "afm_sf_k.bin",
        data_dir,
        data_dir,
        verbose=true
    )
    
    println("\nStep 2: Filtering k=(π,π) from afm_sf_k.bin for bin analysis...")
    println("  (Extract all bin data at the AFM ordering wave vector)")
    
    # Filter the generated afm_sf_k.bin for a specific k-point
    output_afm = filter_bin_file(
        "afm_sf_k.bin",
        (0.0, 0.0),  # k=(0,0)
        dir=data_dir,
        verbose=true
    )
    
    println("\n✓ Generated and filtered AFM structure factor")
    println("  Generated: $(basename(afm_sf_path)) (5760 rows: 576 k-points × 10 bins)")
    println("  Filtered: $(basename(output_afm)) (10 rows: 1 k-point × 10 bins)")
    println("\n💡 Use case: Check if AFM order at (0,0) has converged across bins")
end
