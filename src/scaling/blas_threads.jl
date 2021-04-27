function blas_thread_map_test(model, datatype, indices, priors, base, args...; comargs, cachedirr, scriptdir, blas_thread_max, kargs...)
    plt = plot(
        title = "BLAS Thread Scaling: $model $datatype",
        xguide = "BLAS Threads",
        yguide = "Evaluation Timing",
    )

    avgs = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))
    stds = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))

    println("");flush(stdout)
    println("BLAS Thread Scaling with vendor $(LinearAlgebra.BLAS.vendor()) for $model $datatype");flush(stdout)
    for (mn_id,mn) ∈ enumerate(comargs)
        k = [0.0,0.0];
        asd = BSON.load("$cachedirr/$model/asd-$(mn[1])-$(mn[2]).bson")
        hd  = data_import("$cachedirr/$model/hd-$(mn[1])-$(mn[2]).bson")
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
    bson("$(rootdir)/blas-scaling-$model-$datatype-$blas_thread_max.bson", Dict(:avg=>avgs,:std=>stds,:comargs=>comargs,:datatype=>datatype,:model=>model))
    Plots.pdf(plt, "$(rootdir)/blas-scaling-$model-$datatype-$blas_thread_max")

    nothing
end

function blas_thread_integral_test(model, datatype, indices, priors, base, args...; comargs, scriptdir, blas_thread_max, Neval, pool, kargs...)
    plt = plot(
        title  = "BLAS Thread Scaling: $(fieldname(datatype,1))",
        xguide = "BLAS Threads",
        yguide = "Evaluation Timing",
    )
    cachedir = mkpath("$cachedirr/$model");
    rootdir  = mkpath("$scriptdir/R$(length(readdir(mkpath("$scriptdir")))+1)")

    avgs = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))
    stds = Matrix{Float64}(undef, (blas_thread_max,length(comargs)))

    println("");flush(stdout)
    println("BLAS Thread Scaling with vendor $(LinearAlgebra.BLAS.vendor()) for $model $datatype");flush(stdout)
    for mn ∈ comargs
        k       = [0.0,0.0];
        outdir  = mkpath("$rootdir/$model-$(mn[1])-$(mn[2])-$(round(cθ(mn...)*180/π,digits=2))")
        avg     = Vector{Float64}(undef, blas_thread_max)
        std = Vector{Float64}(undef, blas_thread_max)
        # Thread Testing
        println("Testing Model $mn");flush(stdout)
        for i ∈ 1:blas_thread_max
            LinearAlgebra.BLAS.set_num_threads(i)
            integrate(datatype,indices,priors,base,length(pool),mn,pool,cachedir,outdir)

            stat   = @benchmark integrate($datatype,$indices,$priors,$base,$Neval,$mn,$pool,$cachedir,$outdir)
            times  = sort(sort(stat.times))[10:end-10]
            avg[i] = mean(times)
            std[i] = sqrt(sum((times.-avg[i]).^2)/length(times))

            println("$i Thread(s)------AVG: $(round(avg[i],sigdigits=4))-----");flush(stdout)
        end

        avgs[:,mn_id] .= avg
        stds[:,mn_id] .= std
        hd = data_import("$cachedir/hd-$(mn[1])-$(mn[2]).bson")
        plot!(plt, 1:blas_thread_max, avg/avg[1], label  = "$(size(hd.h_ops.h,1))");
    end

    bson("$(rootdir)/blas-scaling-$model-$datatype-$blas_thread_max.bson", Dict(:avg=>avgs,:std=>stds,:comargs=>comargs,:datatype=>datatype,:model=>model))
    Plots.pdf("$(rootdir)/scaling-test-integral-$blas_thread_max")
    nothing
end
