function map_dim_scaling(model, datatype, indices, priors, base, args...; comargs, cachedirr, scriptdir, plotdir, kargs...)
    # Inputs
    cachedir = mkpath("$cachedirr/$model");
    rootdir  = mkpath("$scriptdir");

    dims = zeros(Float64,length(comargs))
    avg  = zeros(Float64,length(comargs))
    std  = zeros(Float64,length(comargs))

    println("");flush(stdout)
    println("Dimension Scaling $model $datatype");flush(stdout)
    for (i,mn) ∈ enumerate(comargs)
        k = rand(2)
        asd = BSON.load("$cachedir/asd-$(mn[1])-$(mn[2]).bson")
        hd  = data_import("$cachedir/hd-$(mn[1])-$(mn[2]).bson")
        dm  = DataMap(datatype,asd,hd,indices,priors,base)

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
    ttl  = "dim scaling: $model $datatype"

    plt = plot(
    legend = :topleft,
    margin = 5Plots.mm,
    xguide = "Hamiltonian Dimension",
    yguide = "time (ns)",
    frame  = :box,
    title  = ttl
    )

    plot!(plt, dims, avg, ribbon=std, label = "")
    plot!(plt, dimx, avgx, label = "" )
    annotate!((dims[1]+1/3*(dims[end]-dims[1]),avgx[end-10],Plots.text("T = τ dᵅ \n τ: $(round(fitd[:p][1],sigdigits=4)) \n α: $(round(fitd[:p][2],sigdigits=4)) ")))

    Plots.pdf(plt, "$plotdir/$model-$datatype-fit.pdf")
    bson("$rootdir/$model-$datatype-fit.bson",Dict(:dims=>dims, :avg=>avg))

    nothing
end
