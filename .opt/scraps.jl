

function abs_angle_widget(ds::DataSection, di::DataIntegral; outdir=(@__DIR__))

    plt = plot()
    savefunction(i) = begin
        Plots.pdf(plt,mkpath("$outdir/plot")*"/capture-$i.pdf")
    end

    savebutton = button("save")
    on(savefunction,savebutton)

    buttons = [savebutton]

    selector = @manipulate for i ∈ 1:size(ds.chart.data,5)
        n1axis = getindex.(getindex.(ds.chart.base,1),1)[:,1]
        n2axis = getindex.(getindex.(ds.chart.base,2),1)[1,:]

        data_angle = angle.(ds.chart.data[1,1,1,1,i,:,:])
        plt_hm_angle = heatmap(n1axis,n2axis,data_angle,
            frame=:box,
            margins=Plots.mm,
            aspectratio=n1axis[end]/n2axis[end],
            xlims = (n1axis[1],n1axis[end]),
            ylims = (n2axis[1],n2axis[end]),
            color = cgrad(:bilbao),
            colorbar=false
        )

        data_abs = abs.(ds.chart.data[1,1,1,1,i,:,:])
        plt_hm_abs = heatmap(n1axis,n2axis,data_abs,
            frame=:box,
            margins=1Plots.mm,
            aspectratio=n1axis[end]/n2axis[end],
            xlims = (n1axis[1],n1axis[end]),
            ylims = (n2axis[1],n2axis[end]),
            color = cgrad(:bilbao),
            colorbar=false
        )

        base = getindex.(getfield(di.dm.chart,1).base,1)
        plt_spectra = plot(base, abs.(di.data[di.evals.==max(di.evals...)][1])[:],
            label="",
            frame=:box
            )
        plot!(plt_spectra,[base[i],base[i]],[0.0,max(abs.(di.data[di.evals.==max(di.evals...)][1])[:]...)],
            label="",
            )
        plot!(plt_spectra,[base[i],base[i]],[min(abs.(di.data[di.evals.==min(di.evals...)][1])[:]...),max(abs.(di.data[di.evals.==max(di.evals...)][1])[:]...)],
            color=:black,
            label="",
            )

        l = @layout [ grid(1,2) ; a{0.3h}]
        plt = plot(plt_hm_abs,plt_hm_angle,plt_spectra,
            layout=l
            )
    end

    hbox(
        vbox(buttons...),
        selector
    )
end


####
get_asd(args) = SolidStateApps.CommensurateASD(eval(Expr(:call,args[:asd])),args[:mn])

function ds_plot(comps; fs, idx, clims0, layout, plotsize,i)
    # Load the information from inputs
    args = comps[][:args]
    ds = comps[][:ds]
    di = comps[][:di]
    asd = get_asd(args)

    # Make the section integral plot
    int_plt = SolidStateApps.ds_integral_plot(di,fs[1],
        title="$(args[:dtype])",
        size = (300,300)
    )
    SolidStateApps.add_linecut!(int_plt,di,i)

    # Make the Slice Plots
    slice_plts = Vector{typeof(plot())}(undef,length(fs))
    for (i,f) ∈ enumerate(fs)
        ds_max = max(f.(ds.chart.data)...)
        ds_min = min(f.(ds.chart.data)...)
        clims = clims0 == :full ? (ds_min,ds_max) : clims0
        slice_plts[i] = SolidStateApps.ds_slice_plot(ds,f,idx;
            # colorbar = :right,
            clims = clims,
            size = (300,300),
            title="$f"
            # margins = 1Plots.mm
        )
    end

    # Add point labels
    asdg = asd|>SolidState.ASDGeometry
    pt_keys = ["K1","M1","Γ","K'2","M'3"]
    M = inv(asdg["dxtal"][1])
    pts = Dict(pt_keys.=>[M*v for v ∈ getindex.(Ref(asdg["bz_hs"]),pt_keys)])

    map(keys(pts),values(pts)) do k,v
        scatter!(slice_plts[1],[v[1]],[v[2]],color = :black,label="")
        annotate!(slice_plts[1],[(v[1]-0.1,v[2]+0.1,Plots.text(k))],label="")
        plot!(slice_plts[1],[v[1]-0.02,v[1]-0.06],[v[2]+0.02,v[2]+0.06],label="",color=:black)
    end

    #Plot all together
    plot(int_plt,slice_plts...,layout=layout,size=plotsize,margins=1Plots.mm)
end
