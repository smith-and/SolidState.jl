module SolidState
    using Revise
    using Distributed, Dates, OrderedCollections, BSON, Plots, PackageCompiler
    using LinearAlgebra, SharedArrays, StaticArrays
    using CubicSplines, Roots, SpecialFunctions, HCubature
    using Base: Threads

    ### Core Methods
    include("group_theory.jl")

    export HamiltonianOperators, HamiltonianDensity
    export KinematicOperators, KinematicDensity
    include("lattice.jl")
    include("structures.jl")
    include("kinematic_ops.jl")
    include("projectors.jl")

    export TightBindingInfo, TightBindingDensity
    include("tightbinding.jl")

    export TensorChart, DataMap, DataIntegral, DataSection, data_export, data_import
    include("datastructures.jl")
    include("spectral_charts.jl")
    include("optical_charts.jl")

    include("main.jl")


end
