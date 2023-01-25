module SolidState
    using Revise
    using Distributed, Dates, OrderedCollections, BSON, Plots, PackageCompiler
    using LinearAlgebra, SharedArrays, StaticArrays
    using CubicSplines, Roots, SpecialFunctions, HCubature
    using Base: Threads

    ## Core Methods
    include("group_theory.jl")

    export HamiltonianOperators, HamiltonianDensity
    export KinematicOperators, KinematicDensity
    include("lattice.jl")
    include("structures.jl")
    include("kinematic_ops.jl")
    include("projectors.jl")

    ## Model Methods
    ### Tight Binding Methods
    export TightBindingInfo, TightBindingDensity
    include("tightbinding.jl")

    # ### Central Equation Methods
    # include("planewave_expansion.jl")
    #
    # ### Spiral Methods
    # include("spiral_models.jl")


    ## Response Methods
    export TensorChart, DataMap, DataIntegral, DataSection, data_export, data_import
    include("datastructures.jl")

    ### Spectral Response Function 
    include("spectral_charts.jl")

    ### Length Gauge Perturbation Theory Results
    include("optical_charts.jl")

    ## Program Methods
    include("main.jl")



end
