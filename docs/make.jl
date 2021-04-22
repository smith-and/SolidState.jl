using Documenter, SolidState, Git
# include("$(ENV["HOME"])/.julia/config/startup.jl")
# makedocs(sitename="SolidState")

DocMeta.setdocmeta!(SolidState, :DocTestSetup, :(using SolidState); recursive=true)

println("starting with docs!")

makedocs(;
    modules=[SolidState],
    authors="Andrew Smith <asmith.nic@gmail.com> and contributors",
    repo="https://gitlab.com/solidstateapps/SolidState",
    sitename="SolidState",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://solidstateapps.gitlab.io/SolidState",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

println("done with docs!")
