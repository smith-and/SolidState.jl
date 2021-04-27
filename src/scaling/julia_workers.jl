function julia_worker_test(model, datatype, indices, priors, base, args...; comargs, cachedirr, scriptdir, pool = default_worker_pool(), Neval, kargs...)

    wrk_samples = vcat([1], collect(2:2:(length(pool))))
    avg  = Vector{Float64}(undef, length(wrk_samples))
    avgs = Matrix{Float64}(undef,(length(wrk_samples),length(comargs)))

    cachedir = "$cachedirr/$model";
    rootdir0 = mkpath("$scriptdir");
    rootdir  = mkpath("$rootdir0/R$(length(readdir(rootdir0))+1)")
    outdir   = mkpath("$rootdir/T0")

    plt = plot(
            title = "Core Scaling",
            xguide = "1/Ncores",
            yguide = "Tₙ/T₀",
            margin = 10Plots.mm
        );
    println("");flush(stdout)
    println("------Testing Parallel Scaling for $model $datatype ------");flush(stdout)
    for (k,mn) ∈ enumerate(comargs)
        println("------ $mn model ------")
        avg = Vector{Float64}(undef, length(wrk_samples))
        for (i,nw) ∈ enumerate(wrk_samples)
            pool = WorkerPool(workers()[1:nw])

            stat = @timed(integrate(datatype, indices, priors, base, Neval, mn, pool, cachedir, outdir))
            avg[i] = stat.time

            println("$nw worker(s): $(round(avg[i],sigdigits=5))s ");flush(stdout)
        end

        hd  = data_import("$cachedir/hd-$(mn[1])-$(mn[2]).bson")
        plt0 = plot( wrk_samples, avg,
                title = "Worker Scaling - $(size(hd.h_ops.h,1))",
                xguide = "Cores",
                yguide = "Time(s)",
                label  = "",
                margin = 10Plots.mm,
                color  = :red,
                xlims  = (wrk_samples[1],wrk_samples[end]),
                frame  = :box,
            );
        Plots.pdf(plt0, "$rootdir/$model-$datatype-scaling-test-$(mn)")

        plt1 = plot( 1 ./ wrk_samples, avg,
                title = "Core Scaling - $(size(hd.h_ops.h,1)) ",
                xguide = "1/Ncores",
                yguide = "Time(s)",
                label  = "",
                margin = 10Plots.mm,
                color  = :red,
                frame  = :box
            );
        plot!(plt1, 1 ./ wrk_samples, max(avg...) ./ wrk_samples, label  = "",color=:black);
        Plots.pdf(plt1, "$rootdir/$model-$datatype-scaling-test-$(mn)-inv")

        plot!(plt,  1 ./ wrk_samples, avg ./ max(avg...), label  = "$(size(hd.h_ops.h,1))", legend = :topleft);
        avgs[:,k] .= avg
    end

    plot!(plt, 1 ./ wrk_samples, 1 ./ wrk_samples, label  = "ref", color = :black, legend = :topleft);
    Plots.pdf(plt,"$(rootdir)/$model-$datatype-scaling-test-$(nworkers())")
    bson("$rootdir/$model-$datatype-timings.bson",Dict(:comargs=>comargs,:times=>avgs,:wrks=>wrk_samples))

    nothing
end
