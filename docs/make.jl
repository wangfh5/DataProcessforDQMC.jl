import Pkg

# Activate the docs environment (isolated from global)
Pkg.activate(@__DIR__)

try
    # Ensure docs environment uses the local, unregistered package
    Pkg.develop(path=joinpath(@__DIR__, ".."))
    Pkg.instantiate()
catch e
    @warn "Docs environment setup (Pkg.develop/instantiate) failed" e
end

using Documenter
using DataProcessforDQMC

DocMeta.setdocmeta!(DataProcessforDQMC, :DocTestSetup, :(using DataProcessforDQMC); recursive=true)

makedocs(;
    modules=[DataProcessforDQMC],
    authors="wangfh5 <wangfohong@126.com> and contributors",
    sitename="DataProcessforDQMC.jl",
    format=Documenter.HTML(;
        canonical="https://wangfh5.github.io/DataProcessforDQMC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "预编译指南" => "precompilation.md",
        "Bin Analysis与数据处理指南" => "bin_analysis.md",
        "多参数分析指南" => "multi_parameter_analysis.md",
        "API Reference" => "api.md",
    ],
    checkdocs=:none,  # Don't check for missing docstrings (internal modules are not documented)
)

deploydocs(;
    repo="github.com/wangfh5/DataProcessforDQMC.jl",
    devbranch="main",
)

# Optional: Start local preview server
# Usage: DOCS_PREVIEW=true julia docs/make.jl
if get(ENV, "DOCS_PREVIEW", "false") == "true"
    build_dir = joinpath(@__DIR__, "build")
    port = parse(Int, get(ENV, "DOCS_PORT", "8275"))
    
    println("\n" * "="^60)
    println("📚 Starting local documentation preview server...")
    println("🌐 URL: http://localhost:$(port)")
    println("⏹️  Stop: Press Ctrl+C")
    println("="^60 * "\n")
    
    cd(build_dir) do
        run(`python3 -m http.server $(port)`)
    end
end
