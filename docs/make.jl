using FFTMechHomo
using Documenter

DocMeta.setdocmeta!(FFTMechHomo, :DocTestSetup, :(using FFTMechHomo); recursive=true)

makedocs(;
    modules=[FFTMechHomo],
    authors="Felix Mett",
    sitename="FFTMechHomo.jl",
    format=Documenter.HTML(;
        canonical="https://Cr0gan.github.io/FFTMechHomo.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Cr0gan/FFTMechHomo.jl",
    devbranch="master",
)
