using SolidState
using Documenter

DocMeta.setdocmeta!(SolidState, :DocTestSetup, :(using SolidState); recursive=true)

makedocs(;
    modules=[SolidState],
    authors="Andrew Smith",
    repo="https://github.com/smith-and/SolidState.jl/blob/{commit}{path}#{line}",
    sitename="SolidState.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://smith-and.github.io/SolidState.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)


deploydocs(;
    repo="github.com/smith-and/SolidState.jl",
    # devbranch="main",
)
