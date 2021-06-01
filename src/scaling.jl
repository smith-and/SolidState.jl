module Scaling
    using OrderedCollections, LinearAlgebra, Distributed, PackageCompiler, BSON, SharedArrays
    using Plots, LsqFit, Dates, BenchmarkTools

    using ..SolidState
    using ..SolidState: cθ, CommensurateASD,  integrate, cointegrate, data_import

    #######################################
    #### Metric Fitting
    #######################################

    ### Fitting Diagnositics
    function modelfit(xdata,ydata, model, p0)
        p = curve_fit(model, xdata, ydata, p0).param
        f = x->model(x,p)
        Dict(:p =>p, :f=>f, :m=>model,
                :xdata=>xdata, :ydata=>ydata,
                :fdata=>(range(xdata[1],xdata[end],length=100),
                         model(range(xdata[1],xdata[2],length=100),p)
                         )
        )
    end

    function modelfit(xdata,ydata, model, j, p0)
        p = curve_fit(model, j, xdata, ydata, p0).param
        f = x->model(x,p)
        Dict(:p =>p, :f=>f, :m=>model,
                :xdata=>xdata, :ydata=>ydata,
                :fdata=>(collect(range(xdata[1],xdata[end],length=100)),
                         model(collect(range(xdata[1],xdata[2],length=100)),p)
                         )
        )
    end

    @. shiftedpowerlaw(x,p) = (x-p[1])^p[2]+p[3]

    @. powerlaw(x,p) = p[1]*x^p[2]+p[3]

    function powerlawjacobian(x,p)
        @. [ x^p[2] p[1]*p[2]*x^(p[2]-1) ones(length(x)) ]
    end

    @. linearmodel(x,p) = p[1]*x+p[2]

    function eval_su_plot(di)
        xdata = di.evals/di.evals[1]
        fitdict = modelfit(xdata,di.times, linearmodel, ones(2))
        resource_plt = scatter(
            xdata,di.times,
            label="",
            frame=:box,
            yguide="SU(cores⋅h)",
            yrot=60,
            xguide="Function Evaluations",
            title="Evaluations vs. Time",
            margins = 8Plots.mm
        )
        evaldata = range(xdata[1],xdata[end],length=100)|>collect

        plot!(resource_plt,
            evaldata,
            fitdict[:f](evaldata),
            label = "",
            title = "$(round.(fitdict[:p],sigdigits=3))"
        )

        fitdict, resource_plt
    end

    function neval_su_di(di, rootdir, plotdir, plotname)
        if di.evals|>length > 1
            fdict, rplt = eval_su_plot(di)
            Plots.pdf(rplt, "$plotdir/$plotname-su-eval")
            bson("$rootdir/$plotname-su-eval.bson",fdict)
        end
        nothing
    end

    function eval_error_plot(di)
        xdata = di.evals
        ydata = di.err./(di.err[end])

        resource_plt = scatter(
            xdata,ydata,
            label="",
            frame=:box,
            yguide="scaled error: E/E₀",
            yrot=60,
            xguide="scaled evaluations N/N₀",
            title="Evaluations vs. Error",
            margins = 8Plots.mm
        )
        fitdict = modelfit(xdata,ydata,linearmodel, ones(2))
        #evaldata = range(xdata[1],xdata[end],length=100)|>collect
        #=plot!(resource_plt,
            evaldata,
            fitdict[:f](evaldata),
            label = "",
            title = "$(round.(fitdict[:p],sigdigits=3))"
        )
        =#
        fitdict, resource_plt
    end

    function neval_error_di(di, rootdir, plotdir, plotname)
        if di.evals|>length > 1
            fdict, rplt = eval_error_plot(di)
            Plots.pdf(rplt, "$plotdir/$plotname-su-error")
            bson("$rootdir/$plotname-su-error.bson",fdict)
        end
        nothing
    end


    #######################################
    #### Hamiltonian Dimension Scaling
    #######################################
    """
        map_dim_scaling(; asd, datatype, indices, priors, base, comargs, cachedirr, datadir, plotdir, kargs...)

    This method evaluates scaling for the DataMaps associated with the inputs for each
        element in comargs. It reads the model information from `cachedirr`, saves the timing information
        in `datadir` and the plot in `plotdir`. A powerlaw fit is performed and the scaling exponent and principal
        are included in the plot.
    """
    function dm_scaling(; asd, datatype, indices, priors, base, comargs, cachedir, datadir, plotdir, kargs...)
        # Inputs
        modelcache = mkpath("$cachedir/$asd");
        mkpath(datadir);

        dims = zeros(Float64,length(comargs))
        avg  = zeros(Float64,length(comargs))
        std  = zeros(Float64,length(comargs))

        println("");flush(stdout)
        println("Dimension Scaling $asd $datatype");flush(stdout)
        for (i,mn) ∈ enumerate(comargs)
            k = rand(2)
            asd0 = BSON.load("$modelcache/asd-$(mn[1])-$(mn[2]).bson")
            hd  = data_import("$modelcache/hd-$(mn[1])-$(mn[2]).bson")
            dm  = DataMap(datatype,asd0,hd,indices,priors,base)

            dm(k)
            dm(k)
            stats = @benchmark $dm($k)

            dims[i] = size(hd.h_ops.h,1)
            avg[i] = mean(stats.times)
            std[i] = sqrt(sum((stats.times .- avg[i]).^2))/length(stats.times)/2

            println("---$mn----AVG $(round(avg[i],sigdigits=3))------STD $(round(std[i],sigdigits=3))----------------");flush(stdout)
        end

        fitd = modelfit(dims,avg,powerlaw,[1e3,2.0,0.0])
        dimx = range(dims[1],dims[end],length=100)|>collect
        avgx = powerlaw(dimx,fitd[:p])


        bson("$datadir/dim-scaling-$asd-$datatype-bh-$(comargs|>length)-fit.bson",
            Dict(
                :dims=>dims,
                :avg=>avg,
                :std => std,
                :asd => asd,
                :datatype => datatype,
                :comargs => comargs,
                :plotdir => "$plotdir",
                :handle => "dim-scaling-$asd-$datatype-bh-$(comargs|>length)-fit",
                :fitd => fitd,
                :dimx => dimx,
                :avgx => avgx
                )
            )
    end

    #######################################
    #### Worker Number Parallelization Quality/Scaling
    #######################################
    function core_scaling(; asd, datatype, indices, priors, base, comargs, cachedir, datadir, plotdir, pool = default_worker_pool(), Neval, kargs...)

        wrk_samples = vcat([1], collect(2:2:(length(pool))))

        avg  = Vector{Float64}(undef, length(wrk_samples))
        avgs = Matrix{Float64}(undef,(length(wrk_samples),length(comargs)))

        println("");flush(stdout)
        println("------Testing Parallel Scaling for $asd $datatype ------");flush(stdout)
        for (k,mn) ∈ enumerate(comargs)
            println("------ $mn model ------")
            avg = Vector{Float64}(undef, length(wrk_samples))
            for (i,nw) ∈ enumerate(wrk_samples)
                pool = WorkerPool(workers()[1:nw])
                dm = DataMap(asd, mn, datatype, indices, priors, base, cachedir=cachedir)

                ranges=SolidState.ChartInfo(datatype, indices, priors, base, mn)
                if pool|>length|>isequal(1)
                    stat = @timed(integrate(dm, [0.0,0.0],[1.0,1.0], :none, evals = Neval, ranges=ranges, cachedir=cachedir))
                else
                    stat = @timed(integrate(dm, [0.0,0.0],[1.0,1.0], pool, evals = Neval:Neval, ranges=ranges, cachedir="$cachedir/$asd"))
                end
                avg[i] = stat.time

                println("$nw worker(s): $(round(avg[i],sigdigits=5))s ");flush(stdout)
            end

        end

        bson("$datadir/core-scaling-$asd-$datatype-bh-$(comargs|>length)-fit.bson",
            Dict(
                :asd => asd,
                :datatype => datatype,
                :plotdir => plotdir,
                :cachedir => cachedir,
                :handle => "core-scaling-$asd-$datatype-bh-$(comargs|>length)-fit",
                :comargs => comargs,
                :avgs=>avgs,
                :wrk_samples=>wrk_samples,
                :Neval => Neval,
                )
            )
    end

    #######################################
    #### Convergence Testing
    #######################################
    function integral_convergence(; asd, datatype, indices, priors, base, f, comargs, cachedir, datadir, pool = default_worker_pool(), evals, kargs...)

        rootdir0 = mkpath("$datadir");
        rootdir  = mkpath("$(rootdir0)/R$(length(readdir(rootdir0))+1)")

        # outdir = mkpath("$rootdir/$asd-$(mn[1])-$(mn[2])-$(round(cθ(mn...)*180/π,digits=2))")
        # integrate(datatype, indices, priors, base, evals, mn, pool, cachedir, outdir)
        for (k,mn) ∈ enumerate(comargs)
            dm = DataMap(asd, datatype, indices, priors, base; cachedir=cachedir, mn=mn)
            ranges=SolidState.ChartInfo(datatype, indices, priors, base, mn)
            if pool|>length|>isequal(1)
                integrate(dm, [0.0,0.0],[1.0,1.0], :none, evals = evals, ranges=ranges, cachedir=cachedir)
            else
                integrate(dm, [0.0,0.0],[1.0,1.0], pool, evals = evals, ranges=ranges, cachedir="$cachedir/$asd")
            end

            for outdir ∈ readdir("$rootdir",join=true)
                di_plot("$(fieldname(datatype,1))"; title = "$asd $datatype", f=f, rootdir=outdir);
            end
        end

        nothing
    end

    #######################################
    #### Memory Use Test
    #######################################
    function mem_test(; asd, comargs::Vector{Tuple{Int64,Int64}}, datatype, indices, priors, base, cachedir, datadir, plotdir)
        priors  =[(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)]
        a       =[0.0,0.0]
        b       =[1.0,1.0]

        com_dict = OrderedDict{Int,OrderedDict{Symbol,Any}}()

        for (i,mn) ∈ enumerate(comargs)
            di = DataIntegral(
                    DataMap(
                        datatype,
                        BSON.load("$cachedir/$asd/asd-$(mn[1])-$(mn[2]).bson"),
                        data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson"),
                        indices,priors,base
                    ),
                    a,b;
                    ranges=SolidState.ChartInfo(datatype,indices,priors,base,mn),cachedir=cachedir*"/$asd"
            )

            di.dm(rand(2))

            top_out = read(`top -bn1 -p $(getpid())`, String)
            res_size = split(split(top_out,  "\n")[end-1])[6]

            sizes = OrderedDict{Symbol,Any}(
                :mn     => mn,
                :angle  => cθ(mn...)*180/π,
                :dim    => di.dm.d,
                #:asd    => Base.summarysize(asd)/1e9,
                #:hd     => Base.summarysize(hd)/1e9,
                :chart  => Base.summarysize(getfield(di,1).chart)/1e9,
                :di     => Base.summarysize(di)/1e9,
                :top    => res_size,
            )

            println("")
            println("sizes for $mn: ");flush(stdout)
            for p ∈ sizes
                println("$(p.first): $(p.second)")
            end
            push!(com_dict, i=>sizes)

        end


        handle = "mem-scaling-$asd-$datatype-bh-$(comargs|>length)"
        bson("$datadir/$handle.bson", Dict(
                :com_dict => com_dict,
                :handle => handle,
                :comargs => comargs,
            )
        )
    end

    function mem_scaling(; com_dict, handle, plotdir, comargs, args)
        plt = scatter(getindex.(com_dict|>values|>collect,:dim),getindex.(com_dict|>values|>collect,:top),label="")
        # SolidState.make_models(asd, comargs..., cachedir=cachedir)
        Plots.pdf(plt,"$plotdir/$handle")

        plt
    end

    function mem_test(asd, comargs::Tuple{Symbol,Int,Int})

        cachedir = "$(@__DIR__)/../.cache"

        mem_test(asd, SolidState.twist_series(comargs...), cachedir=cachedir)
    end

    #######################################
    #### BLAS Worker Tradeoff
    #######################################
    function blas_julia_thread_tradeoff(; asd, datatype, indices, priors, base, comargs, cachedir, datadir, plotdir, pool = default_worker_pool(), Neval, blas_thread_max=length(pool), kargs...)

        wrk_samples = vcat([1], collect(2:1:(max(length(pool),2))))

        avgs = Array{Float64,3}(undef,(blas_thread_max,length(wrk_samples),length(comargs)))

        rootdir0 = mkpath("$datadir/out");
        rootdir  = mkpath("$rootdir0/R$(length(readdir(rootdir0))+1)")
        cdims    = Vector{Int64}(undef,length(comargs))

        blas_plots = Dict{Int,typeof(plot())}()
        julia_plots = Dict{Int,typeof(plot())}()
        heatmap_plots = Dict{Int,typeof(plot())}()

        println("");flush(stdout)
        println("Determining the BLAS-Julia Thread Tradeoff");flush(stdout)
        for (k,mn) ∈ enumerate(comargs)
            println("------$mn------");flush(stdout)
            hd = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")
            cdims[k] = size(hd.h_ops.h,1)
            for (j,nw) ∈ enumerate(wrk_samples)
                kiddie_pool = WorkerPool(workers()[1:nw])
                for i ∈ 1:blas_thread_max
                    LinearAlgebra.BLAS.set_num_threads(i)

                    ranges=SolidState.ChartInfo(datatype, indices, priors, base, mn)
                    dm = DataMap(asd, datatype, indices, priors, base; cachedir=cachedir, mn=mn)
                    if kiddie_pool|>length|>isequal(1)
                        stat = @timed(integrate(dm, [0.0,0.0],[1.0,1.0], :none, evals = Neval, ranges=ranges, cachedir=cachedir))
                    else
                        stat = @timed(integrate(dm, [0.0,0.0],[1.0,1.0], kiddie_pool, evals = Neval:Neval, ranges=ranges, cachedir="$cachedir/$asd"))
                    end

                    avgs[i,j,k] = stat.time

                    println("$nw worker(s) w/ $i BLAS threads: $(round(avgs[i,j,k],sigdigits=5))s ");flush(stdout)
                end
            end

            #plot the Heat Map
            plt_hm = heatmap(1:blas_thread_max,wrk_samples,avgs[:,:,k]',
                aspectratio = length(wrk_samples)/blas_thread_max,
                frame       = :box,
                margins     = 10Plots.mm,
                title       = "BLAS vs. Julia Thread Tradeoff",
                xguide      = "BLAS Threads",
                yguide      = "Julia Workers",
            )
            Plots.png(plt_hm,"$plotdir/$asd-$datatype-$(mn[1])-$(mn[2])-BJT-tradeoff-heatmap")
            Plots.pdf(plt_hm,"$plotdir/$asd-$datatype-$(mn[1])-$(mn[2])-BJT-tradeoff-heatmap")

            push!(heatmap_plots, k => plt_hm)
        end

        for i ∈ 1:blas_thread_max
            julia_linecuts = plot(title = "Scaling w/ $i BLAS Threads", margin = 10Plots.mm, xguide = "1/Ncores", yguide = "Tₙ/T₀");
            for k ∈ 1:length(comargs)
                #plot the Worker line cut
                plot!(julia_linecuts,  1 ./ wrk_samples, avgs[i,:,k] ./ max(avgs[i,:,k]...), label  = "$(cdims[k])", legend = :topleft);
            end
            plot!(julia_linecuts, 1 ./ wrk_samples, 1 ./ wrk_samples, label  = "ref", color = :black, legend = :topleft);
            Plots.pdf(julia_linecuts, "$plotdir/$asd-$datatype-JT-cut-BT-$i")
            push!(julia_plots, i=>julia_linecuts)
        end

        for j ∈ 1:length(wrk_samples)
            blas_linecuts = plot(title = "BLAS Thread Scaling: $asd $datatype", xguide = "BLAS Threads", yguide = "Evaluation Timing")
            for k ∈ 1:length(comargs)
                #plot the BLAS line cut
                plot!(blas_linecuts, 1:blas_thread_max, avgs[:,j,k]/avgs[1,j,k], label  = "$(cdims[k])");
            end
            Plots.pdf(blas_linecuts, "$plotdir/$asd-$datatype-JT-$j-BT-cut")
            push!(blas_plots, j => blas_linecuts)
        end

        Dict(
            :blas    => blas_plots,
            :julia   => julia_plots,
            :heatmap => heatmap_plots
        )
    end

    #########################################
    #### BLAS Threads
    #########################################
    function blas_thread_map_test(; asd, datatype, indices, priors, base, comargs, cachedirr, scriptdir, blas_thread_max, kargs...)
        plt = plot(
            title = "BLAS Thread Scaling: $asd $datatype",
            xguide = "BLAS Threads",
            yguide = "Evaluation Timing",
        )

        avgs = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))
        stds = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))

        println("");flush(stdout)
        println("BLAS Thread Scaling with vendor $(LinearAlgebra.BLAS.vendor()) for $asd $datatype");flush(stdout)
        for (mn_id,mn) ∈ enumerate(comargs)
            k = [0.0,0.0];
            asd = BSON.load("$cachedirr/$asd/asd-$(mn[1])-$(mn[2]).bson")
            hd  = data_import("$cachedirr/$asd/hd-$(mn[1])-$(mn[2]).bson")
            dm  = DataMap(datatype,asd,hd,indices,priors,base)

            avg = Vector{Float64}(undef, blas_thread_max)
            std = Vector{Float64}(undef, blas_thread_max)

            # Thread Testing
            println("Testing Model $mn");flush(stdout)
            for i ∈ 1:blas_thread_max

                LinearAlgebra.BLAS.set_num_threads(i)
                dm(k);

                stat = @benchmark $dm($k);
                times = sort(sort(stat.times))[10:end-10]
                avg[i] = mean(times)
                std[i]  = sqrt(sum((times.-avg[i]).^2)/length(times))

                println("$i Thread(s)------AVG: $(round(avg[i],sigdigits=4))-----STD: $(round(std[i],sigdigits=4))-----------------");flush(stdout)
            end

            avgs[:,mn_id] .= avg
            stds[:,mn_id] .= std
            plot!(plt, 1:blas_thread_max, avg/avg[1], ribbon = std/avg[1]/2, label  = "$(size(hd.h_ops.h,1))");
        end

        rootdir  = mkpath("$scriptdir/R$(length(readdir("$scriptdir"))+1)")
        bson("$(rootdir)/blas-scaling-$asd-$datatype-$blas_thread_max.bson", Dict(:avg=>avgs,:std=>stds,:comargs=>comargs,:datatype=>datatype,:asd=>asd))
        Plots.pdf(plt, "$(rootdir)/blas-scaling-$asd-$datatype-$blas_thread_max")

        nothing
    end

    function blas_thread_integral_test(asd, datatype, indices, priors, base, args...; comargs, scriptdir, blas_thread_max, Neval, pool, kargs...)
        plt = plot(
            title  = "BLAS Thread Scaling: $(fieldname(datatype,1))",
            xguide = "BLAS Threads",
            yguide = "Evaluation Timing",
        )
        cachedir = mkpath("$cachedirr/$asd");
        rootdir  = mkpath("$scriptdir/R$(length(readdir(mkpath("$scriptdir")))+1)")

        avgs = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))
        stds = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))

        println("");flush(stdout)
        println("BLAS Thread Scaling with vendor $(LinearAlgebra.BLAS.vendor()) for $asd $datatype");flush(stdout)
        for mn ∈ comargs
            k       = [0.0,0.0];
            outdir  = mkpath("$rootdir/$asd-$(mn[1])-$(mn[2])-$(round(cθ(mn...)*180/π,digits=2))")
            avg     = Vector{Float64}(undef, blas_thread_max)
            std = Vector{Float64}(undef, blas_thread_max)
            # Thread Testing
            println("Testing Model $mn");flush(stdout)
            for i ∈ 1:blas_thread_max
                LinearAlgebra.BLAS.set_num_threads(i)

                ranges=SolidState.ChartInfo(datatype, indices, priors, base, mn)
                dm = DataMap(asd, datatype, indices, priors, base; cachedir=cachedir, mn=mn)
                if pool|>length|>isequal(1)
                    stat = @timed(integrate(dm, [0.0,0.0],[1.0,1.0], :none, evals = Neval, ranges=ranges, cachedir=cachedir))
                else
                    stat = @timed(integrate(dm, [0.0,0.0],[1.0,1.0], pool, evals = Neval:Neval, ranges=ranges, cachedir="$cachedir/$asd"))
                end

                # stat   = @benchmark integrate($datatype,$indices,$priors,$base,$Neval,$mn,$pool,$cachedir,$outdir)
                # times  = sort(sort(stat.times))[10:end-10]
                # avg[i] = mean(times)
                # std[i] = sqrt(sum((times.-avg[i]).^2)/length(times))
                avg[i] = stat.time
                std[i] = 0.0

                println("$i Thread(s)------AVG: $(round(avg[i],sigdigits=4))-----");flush(stdout)
            end

            avgs[:,mn_id] .= avg
            stds[:,mn_id] .= std
            hd = data_import("$cachedir/hd-$(mn[1])-$(mn[2]).bson")
            plot!(plt, 1:blas_thread_max, avg/avg[1], label  = "$(size(hd.h_ops.h,1))");
        end

        bson("$(rootdir)/blas-scaling-$asd-$datatype-$blas_thread_max.bson", Dict(:avg=>avgs,:std=>stds,:comargs=>comargs,:datatype=>datatype,:asd=>asd))
        Plots.pdf("$(rootdir)/scaling-test-integral-$blas_thread_max")
        nothing
    end

end
