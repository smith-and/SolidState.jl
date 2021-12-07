using Distributed
    addprocs(2)
    @everywhere using SolidState, OrderedCollections, LinearAlgebra, BSON, Measures

function precompile_main()

    # if ENV["HOME"]=="/jet/home/asmithc"
    #     ENV["scriptdir"] = "/ocean/projects/phy190028p/asmithc/scripts"
    #     ENV["cachedir"]  = "/ocean/projects/phy190028p/asmithc/scripts/.cache"
    # else
    #     ENV["scriptdir"] = "$(ENV["HOME"])/Dropbox/Graduate/scripts"
    #     ENV["cachedir"] = "$(ENV["HOME"])/Dropbox/Graduate/scripts/.cache"
    # end

    SolidState.Main.integral(
        "compilation",
        ASD2,
        (1,2),
        (SHG,
        [(2,2,2)],
        [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.0,0.0,1)],
        [(:ω,0.0,0.0,1)],
        10),
        default_worker_pool(),
    )

    # SolidState.Main.section(
    #     "compilation",
    #     ASD2,
    #     (1,1),
    #     SHG,
    #     [(2,2,2)],
    #     [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.0,0.0,1)],
    #     [(:ω,0.0,0.0,1)],
    #     10,
    #     default_worker_pool(),
    # )
end

precompile_main()
precompile_main()
