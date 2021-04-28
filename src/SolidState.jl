module SolidState
    using Revise
    using Distributed, Dates, OrderedCollections, BSON
    using LinearAlgebra, SharedArrays, StaticArrays
    using Mux, WebIO, Interact, InteractiveUtils
    using CubicSplines, Roots, SpecialFunctions, HCubature
    using Base: Threads

    # Compiling Methods

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

    include("utility/compiler.jl")
    include("utility/serving.jl")
    include("utility/header.jl")

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
        include("scaling/metric.jl")
    end

    # Snipe Spray Scripting Methods
    # include("scripts/snipe_inputs.jl")
    include("scripts/chart_series.jl")
    include("scripts/bands.jl")
    include("scripts/hsp_spectra.jl")

    module Apps
        using Distributed
        using OrderedCollections, LinearAlgebra, Distributed, PackageCompiler, BSON, SharedArrays
        using Plots, BenchmarkTools, LsqFit, Dates
        using Mux, WebIO, Interact, InteractiveUtils

        using ..SolidState
        using ..SolidState: make_models, cθ, CommensurateASD,  integrate, cointegrate, data_import

        # Extraction Methods
        include("apps/completion/extraction.jl")
        include("apps/completion/extraction_plots.jl")
        include("apps/completion/extraction_plotter.jl")

        # Interactive Apps
        include("apps/state_projector.jl")
        include("apps/abs_angle_selector.jl")

        # Component Deck
        include("apps/component_deck/widget_deck.jl")
        include("apps/component_deck/slide_wrapper.jl")
    end



end
