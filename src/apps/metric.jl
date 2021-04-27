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
