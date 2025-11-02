#!/usr/bin/env julia

"""
Example script demonstrating bin file export functionality.

This script shows how to export bin files to CSV and JLD2 formats
for visualization and further analysis in Excel, Python, or Julia.
"""

using DataProcessforDQMC

# Set the data directory
data_dir = joinpath(@__DIR__, "proj_bt_honeycomb_exact_b=36.000_U=4.00_L=24_dtau=0.10_xmag=0.0")

# Example 1: Export single file to CSV
begin
    println("\n📊 Example 1: Export to CSV")
    println("-" ^ 70)
    
    csv_path = export_bin_to_csv("nn_k.bin", data_dir)
    
    println("\n✓ Exported nn_k.bin to CSV")
    println("  Location: exported_csv/nn_k.csv")
    println("  Use case: Open in Excel or load in Python/Julia for plotting")
end

# Example 2: Export single file to JLD2
begin
    println("\n\n💾 Example 2: Export to JLD2")
    println("-" ^ 70)
    
    jld2_path = export_bin_to_jld2("nn_k.bin", data_dir)
    
    println("\n✓ Exported nn_k.bin to JLD2")
    println("  Location: exported_jld2/nn_k.jld2")
    println("  Dataset name: nn_k")
    println("  Use case: Fast loading in Julia for analysis")
end

# Example 3: Export derived quantity (AFM structure factor)
begin
    println("\n\n🔬 Example 3: Export derived quantity")
    println("-" ^ 70)
    println("Step 1: Generate AFM structure factor...")
    
    # Generate AFM structure factor
    afm_sf_path = merge_afm_sf("spsm_k.bin", "afm_sf_k.bin", data_dir, data_dir, verbose=false)
    
    println("Step 2: Export to CSV...")
    csv_path = export_bin_to_csv("afm_sf_k.bin", data_dir)
    
    println("\n✓ Generated and exported AFM structure factor")
    println("  CSV file: exported_csv/afm_sf_k.csv")
    println("  Contains: 5760 rows (576 k-points × 10 bins)")
end

# Example 4: Custom orbital configuration
begin
    println("\n\n🔧 Example 4: Custom orbital pair selection")
    println("-" ^ 70)
    println("Exporting only specific orbital pairs...")
    println()
    
    # Export only diagonal terms (AA and BB)
    csv_path = export_bin_to_csv(
        "nn_k.bin", 
        data_dir,
        orbital_columns=[(3,4), (9,10)],
        orbital_labels=["AA", "BB"]
    )
    
    println("\n✓ Custom export completed")
    println("  Only 2 orbital pairs: AA, BB (diagonal terms)")
    println("  Useful for: analyzing same-orbital correlations")
end

# Example 5: Batch export multiple files
begin
    println("\n\n📦 Example 5: Batch export")
    println("-" ^ 70)
    println("Exporting multiple files to both CSV and JLD2 (default)...")
    println()
    
    exported_files = export_directory_bins(
        data_dir,
        file_patterns=["nn_k.bin", "nn_r.bin", "spsm_k.bin", "b.bin"],
        verbose=false  # output_format=:both is now the default
    )
    
    println("\n✓ Batch export completed")
    println("  Exported $(length(exported_files)) files")
    println("  Default format: CSV + JLD2")
end

# Example 6: Workflow - filter then export
begin
    println("\n\n🔄 Example 6: Filter → Export workflow")
    println("-" ^ 70)
    println("Demonstrating: filter specific k-point → export for plotting")
    println()
    
    println("Step 1: Filter k=(0,0) from afm_sf_k.bin...")
    filtered_file = filter_bin_file(
        "afm_sf_k.bin",
        (0.0, 0.0),
        dir=data_dir,
        verbose=false
    )
    
    println("Step 2: Export filtered data to CSV...")
    csv_path = export_bin_to_csv(
        basename(filtered_file),
        data_dir
    )
    
    println("\n✓ Workflow completed")
    println("  Filtered file: $(basename(filtered_file))")
    println("  CSV file: exported_csv/$(basename(filtered_file) |> x -> replace(x, ".bin" => ".csv"))")
    println("\n💡 Use case: Plot bin convergence at Γ-point")
    println("   - CSV has columns: bin, kx, ky, value_real, value_imag")
    println("   - Ready for Excel or Python matplotlib")
end

println("\n" * "=" ^ 70)
println("✨ All export examples completed!")
println("=" ^ 70)
println("\n📁 Output directories:")
println("  • exported_csv/  - CSV files for Excel/Python")
println("  • exported_jld2/ - JLD2 files for Julia")
println("\n💡 Next steps:")
println("  • Load CSV in Excel: File → Open → Select .csv")
println("  • Load in Python: pandas.read_csv('nn_k.csv')")
println("  • Load in Julia: using JLD2, DataFrames; df = load('nn_k.jld2', 'nn_k')")
println("\n✨ New features:")
println("  • Custom orbital configuration (orbital_columns & orbital_labels)")
println("  • Default batch export format: :both (CSV + JLD2)")
println("  • JLD2 files: one file per bin (nn_k.jld2, nn_r.jld2, etc.)")
println("  • Bin indexing: 1-based (1, 2, 3, ...)")

