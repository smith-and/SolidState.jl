module SolidState
    using Revise
    using Distributed, Dates, OrderedCollections, BSON, Plots
    using LinearAlgebra, SharedArrays, StaticArrays
    using Mux, WebIO, Interact, InteractiveUtils
    using CubicSplines, Roots, SpecialFunctions, HCubature
    using Base: Threads

    ### Core Methods

    export KinematicOperators, KinematicDensity
    export HamiltonianOperators, HamiltonianDensity
    include("ops/lattice.jl")
    include("ops/structures.jl")
    include("ops/kinematic_ops.jl")

    export TightBindingInfo, TightBindingDensity
    include("ops/tightbinding.jl")

    export TensorChart, DataMap, DataIntegral, DataSection, data_export, data_import
    include("ops/datastructures.jl")
    include("ops/spectral_charts.jl")
    include("ops/optical_charts.jl")

    # Execute scaling analysis
    include("scaling.jl")

    # Method Bridge to Bashing
    include("main.jl")

    # Post Processing
    include("extraction.jl")
    include("plots.jl")
    include("figures.jl")

    # Report Processing
    include("report.jl")

    # Submodule for different interactive apps
    include("apps.jl")

end
