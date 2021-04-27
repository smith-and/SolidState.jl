#This is a series integration method based around comargs = [mn1,mn2......]
function integrate_series(datatype::Type{T} where T <: DataChart, indices, priors, base, Neval, f; asd::Symbol, comargs, cachedirr, scriptdir, pool, kargs...)

    cachedir = "$cachedirr/$asd";
    rootdir  = mkpath("$scriptdir/E$(length(readdir(scriptdir))+1)")

    for (k,mn) ∈ enumerate(comargs)
        println("Integrating $asd $(mn[1])-$(mn[2])");flush(stdout)
        #Make a Path for the model results
        outdir = mkpath("$rootdir/$asd-$(mn[1])-$(mn[2])")
        #Do some chart integration
        di = DataIntegral(datatype, indices, priors, base, cachedir, mn)

        did = di(Neval, pool)

        if outdir!="none"
            data_export("$outdir/$(fieldname(datatype,1))-$Neval.bson",did);
        end

    end
end

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

#Generic Plotting Method for Charts
function di_plot(di::DataIntegral; rootdir::String=pwd(),plotdir::String=rootdir, plotname::String=name, f::Function=identity, args...)
    mkpath(plotdir)
    plts = Vector{typeof(plot())}(undef,length(di.err))
    for i∈1:length(di.err)
        plts[i] = plot(
            getindex.(getfield(di.dm.chart,1).base,1),
            f.(di.data[i][:]),
            ribbon = di.err[i],
            xlims  = (getindex.(getfield(di.dm.chart,1).base,1)[1],getindex.(getfield(di.dm.chart,1).base,1)[end]),
            legend = :topleft,
            label  = "$(di.evals[i])",
            frame  = :box,
            margins = 8Plots.mm;
            ylims = (min(0.0,(f.(di.data[i][:]))...),max(f.(di.data[i][:])...,1e-15)),
            args...
        )

        Plots.pdf(plts[i],"$plotdir/$plotname-$(di.evals[i])")
    end

    plts
end
# anim = Animation()
# frame(anim,plts[i])
# gif(anim,"$plotdir/$plotname.gif",fps=4)

#neval_error_di(di, rootdir, plotdir, plotname)
#neval_su_di(di, rootdir, plotdir, plotname)

function add_linecut!(plt_spectra,di,i,f;args...)
    base = getindex.(getfield(di.dm.chart,1).base,1)
    plot!(plt_spectra,[base[i],base[i]],[min(f.(di.data[di.evals.==min(di.evals...)][1])[:]...),max(f.(di.data[di.evals.==max(di.evals...)][1])[:]...)];
        label="",
        args...
        )

    plt_spectra
end

#Generic Plotting Method for Charts
function di_plot(name::String;rootdir::String=pwd(),plotdir::String=rootdir, plotname::String=name, f::Function=identity, args...)
    mkpath(plotdir)
    di = data_import("$rootdir/$name.bson");
    di_plot(di, rootdir=rootdir, plotdir=plotdir, plotname=plotname, f=f; args...)
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
