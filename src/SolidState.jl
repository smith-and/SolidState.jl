module SolidState
    using Revise
    using Distributed, Dates, OrderedCollections, BSON
    using LinearAlgebra, SharedArrays, StaticArrays
    using CubicSplines, Roots, SpecialFunctions, HCubature
    using Base: Threads
    #Generalized Slater-Koster Functions
    include("opt/GSK.jl")

    # Atomic Site Descriptions
    include.(readdir(string(@__DIR__)*"/opt/structures",join=true))
    include("opt/lattice.jl")

    # Hamiltonian Information
    export HamiltonianOperators, HamiltonianDensity
    #include("opt/model0/hamiltonian_ops.jl")
    include("opt/model/hamiltonian_ops.jl")
    include("opt/model/projection_ops.jl")
    include("opt/model/cartan.jl")

    # Tight Binding Models
    export TightBindingInfo, TightBindingDensity
    include("opt/model/TightBinding/tightbinding.jl")
    include("opt/model_making.jl")

    # Position and Velocity Operator Information
    export KinematicOperators, KinematicDensity
    #include("opt/kinematic_operators.jl")
    include("opt/model/kinematic_ops.jl")

    #Domain Information
    include("opt/domains.jl")

    #Data Charts & Maps & Sections
    export TensorChart, DataMap, DataIntegral, DataSection
    export data_export, data_import
    # include("opt/datastructures.jl")
    export logstep

    """
        function logstep(p,n)::AbstractRange
    """
    function logstep(p,n)::AbstractRange
        p^(n-1):p^(n-1):p^n
    end

    export log2step
    """
        function log2step(p::Int,n::Int)::AbstractRange
    """
    function log2step(p::Int,n::Int)::AbstractRange
        p^(n-1):9*p^(n-1):p^n
    end

    export log3step
    function log3step(p::Int,n::Int)::AbstractRange
        p^(n-1):3*p^(n-1):p^n
    end

end
