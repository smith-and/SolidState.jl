#Global Module Loading
using Distributed

#Augment LOAD_PATH
@everywhere include("$(ENV["HOME"])/.julia/config/startup.jl")

@everywhere using SolidState
@everywhere using OrderedCollections, LinearAlgebra

@everywhere using SolidState: DOS, LP, SHG
@everywhere using SolidState.Apps: chart_study, chart_plot

function precompile_main()
    #Set Paths for running
    pargs = SolidStateApps.program_args(:Test, pool=:none, headdir = "$(@__DIR__)");

    #Set Model Set Information
    margs = SolidStateApps.model_args(:BNABBA, (:hull,1,1), pargs)

    #Evaluations for Integration
    Nevals  = nworkers()*(1:1:2);
    priors  = [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)];

    #Integration Chart Arguments
    chart_args = (
        (DOS, [(1,)], priors,(:ω,:bandwidth,400),Nevals, imag),
        (LP,  [(1,1)], priors, [(:ω,0.0,10.0,400)], Nevals,real),
        (SHG, [(2,2,2)], priors, [(:ω,0.0,10.0,400)], Nevals,abs)
    )


    for arg ∈ chart_args
        chart_study(arg,margs)
        chart_plot(arg,margs)
    end

    #Delete the results
    rm(pargs.cachedirr,recursive=true)
    rm("$(pargs.jobname)",recursive=true)
end

precompile_main()
precompile_main()
