using Documenter
using MicrostructureNoise

makedocs(
    modules = [MicrostructureNoise],
    format = :html,
    sitename = "MicrostructureNoise.jl",
    authors = "Moritz Schauer and contributors",
    pages = Any[ 
        "Home" => "index.md",
#        "Manual" => "manual.md",
#        "Library" => "library.md",
        ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/mschauer/MicrostructureNoise.jl.git",
    julia  = "0.7",
    target = "build",
    deps = nothing,
    make = nothing,
)
