using Documenter, Literate, SolidState, Hop

# Make Documentation Source
Literate.markdown("$(@__DIR__)/src/dev.jl","$(@__DIR__)/src")

# Make the Documentation Pages
println("starting with docs!")
DocMeta.setdocmeta!(SolidState, :DocTestSetup, :(using Documenter); recursive=true)
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
    pages = [
        "Homepage" => "index.md",
        "ReferenceCheck" => "dev.md"
    ],
)
println("done with docs!")
