#!/usr/bin/env julia

# Script to generate SSH deploy key for Documenter.jl
# Run this once to set up automatic documentation deployment

using DocumenterTools

println("="^60)
println("Generating SSH deploy key for documentation deployment")
println("="^60)
println()

try
    DocumenterTools.genkeys(; user="wangfh5", repo="DataProcessforDQMC.jl")
    
    println()
    println("="^60)
    println("✅ Key generation successful!")
    println("="^60)
    println()
    println("Next steps:")
    println("1. Go to: https://github.com/wangfh5/DataProcessforDQMC.jl/settings/secrets/actions")
    println("2. Click 'New repository secret'")
    println("3. Name: DOCUMENTER_KEY")
    println("4. Value: [Copy the private key shown above]")
    println("5. Click 'Add secret'")
    println()
    println("The public key has been automatically added to your GitHub deploy keys.")
    println("="^60)
catch e
    println("Error: $e")
    println()
    println("Make sure you have DocumenterTools installed:")
    println("  using Pkg; Pkg.add(\"DocumenterTools\")")
end

