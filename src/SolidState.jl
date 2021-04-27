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
    include("opt/model/kinematic_ops.jl")

    #Domain Information
    include("opt/domains.jl")

    #Data Charts & Maps & Sections
    export TensorChart, DataMap, DataIntegral, DataSection
    export data_export, data_import
    include("opt/datastructures.jl")
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


    module Scaling
        using OrderedCollections, LinearAlgebra, Distributed, PackageCompiler, BSON, SharedArrays
        using Plots, LsqFit, Dates, BenchmarkTools

        using ..SolidState
        using ..SolidState: make_models, cθ, CommensurateASD,  integrate, cointegrate, data_import

        #Include Test Functions
        include("scaling/dim_scaling.jl")
        include("scaling/loop_test.jl")
        include("scaling/blas_threads.jl")
        include("scaling/convergence_test.jl")
        include("scaling/julia_workers.jl")
        include("scaling/blas_julia_thread_tradeoff.jl")
    end

    module Apps
        using Distributed
        using OrderedCollections, LinearAlgebra, Distributed, PackageCompiler, BSON, SharedArrays
        using Plots, Mux, WebIO, Interact, InteractiveUtils, BenchmarkTools, LsqFit, Dates

        using ..SolidState
        using ..SolidState: make_models, cθ, CommensurateASD,  integrate, cointegrate, data_import

        #Running Methods
        include("apps/compiler.jl")
        include("apps/header.jl")
        include("apps/serving.jl")

        #Runtime Methods
        include("apps/metric.jl")
        include("apps/integral_plot.jl")
        include("apps/hsp_spectra.jl")

        #Band Methods
        include("apps/bands/bands.jl")
        include("apps/bands/state_projector.jl")

        #Extraction Methods
        include("apps/extraction/extraction.jl")
        include("apps/extraction/extraction_plots.jl")
        include("apps/extraction/extraction_plotter.jl")

        # DataSection Selector
        include("apps/ds_selector/abs_angle_selector.jl")

        #Component Deck
        include("apps/component_deck/widget_deck.jl")
        include("apps/component_deck/slide_wrapper.jl")

        # Snipe Spray Interface
        include("apps/snipe_spray/snipe_inputs.jl")

    end

end
