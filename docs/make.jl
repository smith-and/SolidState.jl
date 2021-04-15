using Documenter, SolidState
# include("$(ENV["HOME"])/.julia/config/startup.jl")
# makedocs(sitename="SolidState")

DocMeta.setdocmeta!(SolidState, :DocTestSetup, :(using SolidState); recursive=true)

println("starting with docs!")

makedocs(;
    modules=[SolidState],
    authors="Andrew Smith <asmith.nic@gmail.com> and contributors",
    # repo="https://gitlab.com/smith-and/SolidState",
    sitename="SolidState",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        # canonical="https://smith-and.gitlab.io/SolidState",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

println("done with docs!")
