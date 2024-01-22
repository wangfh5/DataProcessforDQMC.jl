using DataProcessforDQMC
using Documenter

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
    ],
)

deploydocs(;
    repo="github.com/wangfh5/DataProcessforDQMC.jl",
    devbranch="main",
)
