function integral_convergence(model, datatype, indices, priors, base, f; comargs, cachedirr, scriptdir, pool = default_worker_pool(), evals = 1000:20000:300000, kargs...)
    cachedir = "$cachedirr/$model";
    rootdir0 = mkpath("$scriptdir");
    rootdir  = mkpath("$(rootdir0)/R$(length(readdir(rootdir0))+1)")

    for (k,mn) ∈ enumerate(comargs)
        outdir = mkpath("$rootdir/$model-$(mn[1])-$(mn[2])-$(round(cθ(mn...)*180/π,digits=2))")
        integrate(datatype, indices, priors, base, evals, mn, pool, cachedir, outdir)

        for outdir ∈ readdir("$rootdir",join=true)
            di_plot("$(fieldname(datatype,1))"; title = "$model $datatype", f=f, rootdir=outdir);
        end
    end

    nothing
end
