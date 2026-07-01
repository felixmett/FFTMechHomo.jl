using Documenter
using DocumenterCitations
using FFTMechHomo

DocMeta.setdocmeta!(
    FFTMechHomo,
    :DocTestSetup,
    :(using FFTMechHomo)
)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(
    modules = [FFTMechHomo],
    plugins = [bib],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://felixmett.github.io/FFTMechHomo.jl"
    ),
    sitename = "FFTMechHomo.jl",
    authors = "Felix Mett",
    pages = [
        "Home" => "index.md",
        "Background" => "background.md",
        "Examples" => [
			"Custom Material Definition" => "examples/custommaterial.md",
			"Effective Elasticity Tensor" => "examples/effectivetensor.md",
		],
        "API" => [
            "Overview" => "api/index.md",
            "Discretization" => "api/discretization.md",
            "Microstructure" => "api/microstructure.md",
            "Materials" => "api/material.md",
            "Solver" => "api/solver.md",
        ],
        "References" => "references.md",
    ],
    warnonly = false,
    draft = false,
    source = "src",
    build = "build",
)

deploydocs(
    repo = "github.com/felixmett/FFTMechHomo.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "master",
    push_preview = true,
)