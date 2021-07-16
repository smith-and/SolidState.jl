#Global Module Loading
using Distributed

#Augment LOAD_PATH
@everywhere using SolidState

function precompile_main()
    #Set Paths for running
    args = (
        "compilation",
        ASD2,
        (1,1),
        SHG,
        [(2,2,2)],
        [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.0,0.0,1)],
        [(:ω,0.0,0.0,1)],
        10,
        default_worker_pool(),
    )

    SolidState.Main.section(args...)

    args = (
        "compilation",
        ASD2,
        (1,2),
        SHG,
        [(2,2,2)],
        [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.0,0.0,1)],
        [(:ω,0.0,0.0,1)],
        10,
        default_worker_pool(),
    )

    SolidState.Main.integral(args...)

    args = (
        "compilation",
        ASD2,
        (1,1),
    )

    SolidState.Main.bands(args...)


end

precompile_main()
precompile_main()
