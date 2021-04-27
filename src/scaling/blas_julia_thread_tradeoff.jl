function blas_julia_thread_tradeoff(model, datatype, indices, priors, base, args...; comargs, cachedirr, scriptdir, pool = default_worker_pool(), Neval, blas_thread_max, kargs...)

    wrk_samples = vcat([1], collect(2:2:(max(length(pool),2))))

    avgs = Array{Float64,3}(undef,(blas_thread_max,length(wrk_samples),length(comargs)))

    cachedir = "$cachedirr/$model";
    rootdir0 = mkpath("$scriptdir/out");
    rootdir  = mkpath("$rootdir0/R$(length(readdir(rootdir0))+1)")
    outdir   = mkpath("$rootdir/T0")
    cdims    = Vector{Int64}(undef,length(comargs))
    
    println("");flush(stdout)
    println("Determining the BLAS-Julia Thread Tradeoff");flush(stdout)
    for (k,mn) ∈ enumerate(comargs)
        println("------ $mn model ------");flush(stdout)
        hd = data_import("$cachedir/hd-$(mn[1])-$(mn[2]).bson")
        cdims[k] = size(hd.h_ops.h,1)
        for (j,nw) ∈ enumerate(wrk_samples)
            kiddie_pool = WorkerPool(workers()[1:nw])
            for i ∈ 1:blas_thread_max
                LinearAlgebra.BLAS.set_num_threads(i)

                stat = @timed(integrate(datatype, indices, priors, base, Neval, mn, kiddie_pool, cachedir, outdir))
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
        Plots.png(plt_hm,"$rootdir/$model-$datatype-$(mn[1])-$(mn[2])-BJT-tradeoff-heatmap")
        Plots.pdf(plt_hm,"$rootdir/$model-$datatype-$(mn[1])-$(mn[2])-BJT-tradeoff-heatmap")
    end

    for (k,mn) ∈ enumerate(comargs)
        julia_linecuts = plot(title = "Core Scaling: dim $(cdims[k])", margin = 10Plots.mm, xguide = "1/Ncores", yguide = "Tₙ/T₀");
        for i ∈ 1:blas_thread_max
            #plot the Worker line cut
            plot!(julia_linecuts,  1 ./ wrk_samples, avgs[i,:,k] ./ max(avgs[i,:,k]...), label  = "$i", legend = :topleft);
        end
        plot!(julia_linecuts, 1 ./ wrk_samples, 1 ./ wrk_samples, label  = "ref", color = :black, legend = :topleft);
        Plots.pdf(julia_linecuts, "$(rootdir)/$model-$datatype-JT-cut-$(mn[1])-$(mn[2])")
    end

    for i ∈ 1:blas_thread_max
        julia_linecuts = plot(title = "Scaling w/ $i BLAS Threads", margin = 10Plots.mm, xguide = "1/Ncores", yguide = "Tₙ/T₀");
        for k ∈ 1:length(comargs)
            #plot the Worker line cut
            plot!(julia_linecuts,  1 ./ wrk_samples, avgs[i,:,k] ./ max(avgs[i,:,k]...), label  = "$(cdims[k])", legend = :topleft);
        end
        plot!(julia_linecuts, 1 ./ wrk_samples, 1 ./ wrk_samples, label  = "ref", color = :black, legend = :topleft);
        Plots.pdf(julia_linecuts, "$(rootdir)/$model-$datatype-JT-cut-BT-$i")
    end

    for j ∈ 1:length(wrk_samples)
        blas_linecuts = plot(title = "BLAS Thread Scaling: $model $datatype", xguide = "BLAS Threads", yguide = "Evaluation Timing")
        for k ∈ 1:length(comargs)
            #plot the BLAS line cut
            plot!(blas_linecuts, 1:blas_thread_max, avgs[:,j,k]/avgs[1,j,k], label  = "$(cdims[k])");
        end
        Plots.pdf(blas_linecuts, "$(rootdir)/$model-$datatype-JT-$j-BT-cut")
    end

    nothing
end
