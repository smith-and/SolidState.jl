#This is a series integration method based around comargs = [mn1,mn2......]


#Wrapper function for
function chart_study(int_args, model_args)
    #Print Messages
    println("")
    println("Doing Chart Integral");

    #asd Comargs
    dict_print(int_args)

    #Do integration
    println("Doing Integration")
    @time integrate_series(int_args...; model_args...)

    nothing
end

#Open Plots, has same args...;kargs...
function plot_integrals(datatype::BC, indices, priors, base, Neval, f; asd::Symbol, comargs, cachedirr, scriptdir, plotdir, RN=length(readdir(scriptdir)),  pool, kargs...)
    cachedir = "$cachedirr/$asd";
    rootdir  = mkpath("$scriptdir/E$(RN)")
    plotdirr = mkpath("$plotdir/E$(RN)")
    for (i,outdir) ∈ enumerate(readdir("$rootdir"))
        println("Plotting Charts: $outdir");flush(stdout)
        #Extract the Commensurate Indices
        idxs = findall(isequal('-'),outdir)
        m = parse(Int,outdir[(idxs[1]+1):(idxs[2]-1)])
        n = parse(Int,outdir[(idxs[2]+1):(end)])
        #Plot the Chart Information
        mkpath(plotdir)
        args = (
            name     = "$(fieldname(datatype,1))-$Neval",
            rootdir  ="$rootdir/$outdir",
            plotdir  = "$plotdirr/$outdir",
            plotname = "$(fieldname(datatype,1))",
            title    = "$asd $datatype - $(round(cθ(m,n)*180/π,digits=2))",
            f        = f
        )
        di = data_import("$rootdir/$name.bson");

        plt = plot(frame=:box)
        for i ∈ 1:size(di.data[1],1)
            plot!( f.(di.data[1][i, :]) ; linetype=:steppre, label = "")
        end

        println(f.(di.data[1][i, :]))

        Plots.pdf(plt,"$(args.plotdir)/$(args.plotname)-$(di.evals[i])")
    end

    nothing
end

function plot_integrals(datatype::Type{T} where T <: DataChart, indices, priors, base, Neval, f; asd::Symbol, comargs, cachedirr, scriptdir, plotdir, RN=length(readdir(scriptdir)),  pool, kargs...)
    cachedir = "$cachedirr/$asd";
    rootdir  = mkpath("$scriptdir/E$(RN)")
    plotdirr = mkpath("$plotdir/E$(RN)")
    for (i,outdir) ∈ enumerate(readdir("$rootdir"))
        println("Plotting Charts: $outdir");flush(stdout)
        #Extract the Commensurate Indices
        idxs = findall(isequal('-'),outdir)
        m = parse(Int,outdir[(idxs[1]+1):(idxs[2]-1)])
        n = parse(Int,outdir[(idxs[2]+1):(end)])
        #Plot the Chart Information
        di_plot("$(fieldname(datatype,1))-$Neval";
            rootdir="$rootdir/$outdir",
            plotdir = "$plotdirr/$outdir",
            plotname = "$(fieldname(datatype,1))",
            title = "$asd $datatype - $(round(cθ(m,n)*180/π,digits=2))",
            f=f
        )
    end

    nothing
end

#Convenience Wrapper for plot_integrals
function chart_plot(int_args, model_args)
    #Plot integration
    println("Plotting Integration");flush(stdout)
    @time plot_integrals(int_args...; model_args...)
    nothing
end



#Print Args
function dict_print(args)
    for  (i,key) ∈ enumerate(keys(args))
        println("$key: $(args[i])")
    end
    flush(stdout)
end

#Program Args
function program_args(jobname::Symbol; headdir = pwd(), cachedir="$headdir/.cache", pool::Union{AbstractWorkerPool,Symbol}=default_worker_pool(), key="R", kargs...)
     # asd header dir

    if typeof(pool)==Symbol
        poolz=pool
    elseif length(pool)==0
        addprocs(1)
        poolz = default_worker_pool()
    else
        poolz = pool
    end

    keyN = length(filter(x->occursin(key,x),readdir(mkpath("$headdir/$jobname/out"))))+1

    kargs = (
        jobname          = "$headdir/$jobname",
        cachedirr       = mkpath(cachedir),
        scriptdir       = mkpath("$headdir/$jobname/out/$key$keyN"),
        plotdir         = mkpath("$headdir/$jobname/plot/$key$keyN"),
        pool            = poolz, #:none,
    )
end

#Model Args
function model_args(asd::Symbol,(series,minn,maxn)::Tuple{Symbol,Int,Int},pargs)
    margs = (
        asd             = asd,
        comargs         = make_models(asd,series,minn,maxn,cachedir=pargs.cachedirr),
    )
    args = merge(pargs, margs)
    dict_print(args)

    args
end

function main_script(asd::Symbol, comargs::Tuple{Symbol,Int,Int}, datatype, indices, base, f, evals::AbstractRange; name=asd, saving = true, pool::Union{AbstractWorkerPool,Symbol}=default_worker_pool())
    #Data Chart Argument
    priors = [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)]
    arg    = (datatype, indices, priors, base, evals, f)
    #Set Paths for running
    pargs   = program_args(name,
        headdir  = "$(@__DIR__)/..",
        cachedir = "$(@__DIR__)/../.cache",
        pool     = pool,
        key      = "$datatype"
    );
    margs = model_args(asd, comargs, pargs) #Set Model Set Information

    #Compute & Plot Chart Integrals
    chart_study(arg,margs)
    chart_plot(arg,margs)

    #Delete the results
    !saving && rm.(("$(pargs.scriptdir)","$(pargs.plotdir)"),recursive=true)
end

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
