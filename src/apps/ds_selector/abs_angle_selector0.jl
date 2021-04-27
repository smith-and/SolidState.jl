function abs_angle_ds_plot(ds::DataSection,idx)
    n1axis = getindex.(getindex.(ds.chart.base,1),1)[:,1]
    n2axis = getindex.(getindex.(ds.chart.base,2),1)[1,:]

    data_angle = angle.(ds.chart.data[idx...])
    plt_hm_angle = heatmap(n1axis,n2axis,data_angle,
        frame=:box,
        margins=Plots.mm,
        aspectratio=n1axis[end]/n2axis[end],
        xlims = (n1axis[1],n1axis[end]),
        ylims = (n2axis[1],n2axis[end]),
        color = cgrad(:bilbao),
        # colorbar=false,
        size = (300, 300),
    )

    data_abs = abs.(ds.chart.data[idx...])
    plt_hm_abs = heatmap(n1axis,n2axis,data_abs,
        frame=:box,
        margins=1Plots.mm,
        aspectratio=n1axis[end]/n2axis[end],
        xlims = (n1axis[1],n1axis[end]),
        ylims = (n2axis[1],n2axis[end]),
        color = cgrad(:bilbao),
        # colorbar=false,
        size = (300, 300),
    )

    (plt_hm_abs,plt_hm_angle)
end

function ds_integral_plot(di::DataIntegral,i)
    base = getindex.(getfield(di.dm.chart,1).base,1)
    plt_spectra = plot(base, abs.(di.data[di.evals.==max(di.evals...)][1])[:],
        label="",
        frame=:box,
        size = (300, 300),
        )
    plot!(plt_spectra,[base[i],base[i]],[min(abs.(di.data[di.evals.==min(di.evals...)][1])[:]...),max(abs.(di.data[di.evals.==max(di.evals...)][1])[:]...)],
        color=:black,
        label="",
        )

    plt_spectra
end

function abs_angle_base_callback(comps)
    i = comps[:selector][]
    ds,di = comps[:dsdi]

    idx = (1,1,1,1,i,:,:)
    (plt_hm_abs, plt_hm_angle) = abs_angle_ds_plot(ds,idx)

    plt_spectra = ds_integral_plot(di,i)

    plt = plot(plt_spectra, plt_hm_abs, plt_hm_angle,
        layout=(@layout grid(1,3)),
        size = (900, 300)
    )
end

function abs_angle_widget(args, ds::DataSection; plotdir=pwd())

    components = Dict(
        :dsdi       => (ds,integrate(ds)),
        :update     => button("update")
    )

    push!(components,:selector   => slider(1:size(ds.chart.data,5)))
    on(components[:selector]) do i
        components[:update][]+=1
    end

    output = map(components[:update]) do i
        abs_angle_base_callback(components)
    end

    w = Widget(components, output=output)
    @layout! w vbox(:selector,observe(_))

    # node(:div,
    #     node(:div, dom"h1"("$(args[:asd]), $(args[:mn]), $(args[:dtype]), $(args[:N])")),
    #     node(:div, w)
    # )
    w
end

function abs_angle_selector(args; force=false)
    if isfile("$(args[:outdir])/$(args[:asd])/out/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])/ds-$(args[:N]).bson") && !force
        ds = SolidState.data_import("$(args[:outdir])/$(args[:asd])/out/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])/ds-$(args[:N]).bson")
    else
        ds = ds_calc(; args...)
    end

    abs_angle_widget(args, ds; plotdir="$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])")
end


function scope_selector((sym_x,min_x,max_x,N_x))
    comps = Dict(
        :update=>button("update"),
        :sym => textbox(sym_x ; value = sym_x),
        :min => spinbox(min_x ; value = min_x),
        :max => spinbox(max_x ; value = max_x),
        :N   => spinbox(N_x   ; value = N_x)
    )

    map(comps|>keys|>collect) do key
        if key!=:update
            map(comps[key]) do value
                comps[:update][]+=1
            end
        end
    end

    output = map(comps[:update]) do i
        (
            Symbol(comps[:sym][]),
            comps[:min][],
            comps[:max][],
            Int(floor(comps[:N][]))
        )
    end

    w = Widget(comps,output=output)

    @layout! w hbox(hbox(:sym,:min,:max,:N), observe(_))

    w
end

function scope_collector()

    comps = Dict(
            :update => button("update"),
            :add_scope => button("+scope"),
            :scopes => OrderedDict{Int,Any}(),
            :values => Dict(),
    )

    on(comps[:add_scope]) do i
        push!(comps[:scopes],i=>scope_selector(("x",0.0,0.0,1)))
        comps[:update][]+=1
    end

    output = map(comps[:update]) do i
        vbox(collect(values(comps[:scopes]))...)
    end

    w = Widget(comps, output=output)

    @layout! w vbox(:add_scope,observe(_))
end

# ps = [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)]
# vbox(scope_selector.(ps)...)

function index_selector()
    comps = Dict(
        :idx=>textbox(value = "2,2,2"),
        :update=>button("update")
    )

    output = map(comps[:idx]) do x
        char_idx = filter((comps[:idx][]...,)) do c
            try
                parse(Int,c)
                true
            catch
                false
            end
        end
        parse.(Int,char_idx)
    end


    w = Interact.Widget(comps,output=output)

    @layout! w hbox(:idx,observe(_))
end



function widget_input()
    chart_types = vcat((InteractiveUtils.subtypes(SolidState.DataChart).|>InteractiveUtils.subtypes)...)
    comps = Dict(
        :update => button("update"),
        :cachedir => textbox(value="/home/smitha/Dropbox/Graduate/packages/.scripts/local/.cache"),
        :outdir => textbox(value="$(pwd())"),
        :asd     => dropdown([:ASD2,:ASD4,:ASD2B]),
        :mn      => dropdown([[(1,1)],[(1,2)]]),
        :dtype   => dropdown(chart_types),
        :N       => dropdown([10,100,200]),
        :indices => dropdown([[(1,)],[(2,2,2)],[(1,1)]]),
        :priors  => dropdown([[(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)]]),
        :base    => dropdown([[(:ω,0.0,5.0,300)],[(:ω,0.0,8.0,300)],(:ω,:bandwidth,300)]),
        :datafile => dropdown([:none]),
    );

    map(keys(comps),values(comps)) do k,v
        if k == :update
            nothing
        else
            on(v) do i
                comps[:update][]+=1
            end
        end
    end

    output = map(comps[:update]) do i
        out = Dict(keys(comps).=>getindex.(values(comps)))
        in(:mn,keys(out)) ? out[:mn] = out[:mn][1] : nothing
        out
    end

    w = Widget(comps,output=output)

    @layout! w hbox(:asd,:mn,:dtype,:N,:widgetf,:indices,:priors,:base)
end

function ds_calc(; asd::Symbol, mn, dtype, indices, priors, base, N, offset = (0,0), cachedir::String=:none, outdir=:auto, kargs...)
    dom     = [(:n1,-0.5+offset[1],0.5+offset[1],N),(:n2,-0.5+offset[2],0.5+offset[2],N)]

    if cachedir==:none
        dm = DataMap(dtype, SolidState.CommensurateASD(asd(),mn), indices, priors, base)
    else
        dm = DataMap(asd, dtype, indices, priors, base; mn=mn, cachedir=cachedir)
    end
    ds = DataSection(dm,dom);
    ds()
    SolidState.data_export(mkpath("$outdir/$asd/out/$dtype/$asd-$(mn[1])-$(mn[2])")*"/ds-$N.bson",ds)

    ds

end
