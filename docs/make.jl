using Documenter
using MicrostructureNoise
DocMeta.setdocmeta!(MicrostructureNoise, :DocTestSetup, :(using MicrostructureNoise); recursive=true)

makedocs(;
    modules=[MicrostructureNoise],
    authors="mschauer <moritzschauer@web.de> and contributors",
    sitename="MicrostructureNoise.jl",
    format=Documenter.HTML(;
        canonical="https://mschauer.github.io/MicrostructureNoise.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mschauer/MicrostructureNoise.jl",
    devbranch="main",
)
