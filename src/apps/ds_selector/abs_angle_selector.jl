
### Input widget

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
        parsed_input = strip.(isequal('('),strip.(isequal(')'),split(x,isequal(','))))

        y = map(parsed_input) do c
            try
                parse(Int,c)
            catch
                try
                    parse(Float64,c)
                catch
                    "nah"
                end
            end
        end
        (y...,)
    end


    w = Interact.Widget(comps,output=output)

    @layout! w hbox(:idx,observe(_))
end

function arg_asd(arg)
    SolidState.CommensurateASD(eval(Expr(:call,Symbol(arg[:asd]))),arg[:mn])
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

function ds_input()
    chart_types = vcat((InteractiveUtils.subtypes(SolidState.DataChart).|>InteractiveUtils.subtypes)...)
    comps = Dict{Symbol,Any}(
        :update => button("update"),
        :load   => button("load"),
        :force  => checkbox("force"),
        :cachedir => textbox(value="/home/smitha/Dropbox/Graduate/packages/.scripts/local/.cache"),
        :outdir => textbox(value="$(pwd())"),
        :asd     => dropdown([:ASD2,:ASD4,:ASD2B]),
        :mn      => dropdown([[(1,1)],[(1,2)]]),
        :dtype   => dropdown(chart_types, value = DOS),
        :N       => dropdown([10,100,200]),
        :indices => dropdown([[(1,)],[(2,2,2)],[(1,1)]]),
        :priors  => dropdown([[(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)]]),
        :base    => dropdown([[(:ω,0.0,5.0,300)],[(:ω,0.0,8.0,300)],(:ω,:bandwidth,300)]),
        :datafile => dropdown([:none]),
    );

    on(comps[:load]) do i
        arg_keys = filter!(x->((x!=:ds)&(x!=:update)&(x!=:args)&(x!=:di)),keys(comps)|>collect)
        args = Dict{Symbol,Any}(arg_keys.=>getindex.(getindex.(Ref(comps),arg_keys)))
        in(:mn,keys(args)) ? args[:mn] = args[:mn][1] : nothing

        push!(comps ,:args => args)
        if isfile("$(args[:outdir])/$(args[:asd])/out/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])/ds-$(args[:N]).bson") && !(args[:force])
            println("Loaded Calculation")
            push!(comps ,:ds =>  SolidState.data_import("$(args[:outdir])/$(args[:asd])/out/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])/ds-$(args[:N]).bson"))
            push!(comps,:di=>SolidState.integrate(comps[:ds]))
        else
            println("Did Calculation")
            push!(comps , :ds => ds_calc(; args...))
            push!(comps,:di=>SolidState.integrate(comps[:ds]))
        end
        comps[:update][]+=1
        nothing
    end

    map(keys(comps),values(comps)) do k,v
        if k == :update || k == :ds || k == :args || k == :di
            nothing
        else
            on(v) do i
                comps[:update][]+=1n
            end
        end
    end

    output = map(comps[:update]) do i
        try
            comps
        catch
            "nothing loaded"
        end
    end

    w = Widget(comps,output=output)

    @layout! w vbox(hbox(:force,:load,:N,:asd,:mn,),hbox(:dtype,:indices,:priors,:base))
end


### Plotting


function ds_slice_plot(ds::DataSection,f,idx; args...)
    n1axis = getindex.(getindex.(ds.chart.base,1),1)[:,1]
    n2axis = getindex.(getindex.(ds.chart.base,2),1)[1,:]
    # pyplot()

    data = f.(ds.chart.data[idx...])
    plt_hm = heatmap(n1axis,n2axis,data;
        frame=:box,
        margins=Plots.mm,
        aspectratio=n1axis[end]/n2axis[end],
        xlims = (n1axis[1],n1axis[end]),
        ylims = (n2axis[1],n2axis[end]),
        color = cgrad(:bilbao),
        xticks = :none,
        yticks = :none,
        size = (300, 300),
        args...
    )

    plt_hm
end


function ds_integral_plot(di::DataIntegral,f; args...)
    # pyplot()

    base = getindex.(getfield(di.dm.chart,1).base,1)
    plt_spectra = plot(base, f.(di.data[di.evals.==max(di.evals...)][1])[:];
        label="",
        frame=:box,
        size = (300, 300),
        args ...
        )

    plt_spectra
end

function add_hsp_labels!(slice_plts, asdg, pt_keys)
    # Add point labels
    M = inv(asdg["dxtal"][1])
    pts = Dict(pt_keys.=>[M*v for v ∈ getindex.(Ref(asdg["bz_hs"]),pt_keys)])

    map(keys(pts),values(pts)) do k,v
        scatter!(slice_plts[1],[v[1]],[v[2]],color = :black,label="")
        annotate!(slice_plts[1],[(v[1]-0.1,v[2]+0.1,Plots.text(k))],label="")
        plot!(slice_plts[1],[v[1]-0.02,v[1]-0.06],[v[2]+0.02,v[2]+0.06],label="",color=:black)
    end

    slice_plts
end


### Sample Selector

function ds_base_callback(comps)
    i = comps[:selector][]
    ds,di = comps[:dsdi]

    idx = (1,1,1,1,i,:,:)
    (plt_hm_abs, plt_hm_angle) = ds_slice_plot(ds,idx)

    plt_spectra = ds_integral_plot(di,i)

    l = Plots.@layout [ a [ b ; c] ]
    plt = plot(plt_spectra, plt_hm_abs, plt_hm_angle,
        layout=l,
        size = (900, 300)
    )
end

function abs_angle_widget(args, ds::DataSection)

    components = Dict{Symbol,Any}(
        :input   => ds_input(),
    )

    on(components[:input]) do input_obs
        args = components[:input][][:args]
        ds = components[:input][][:ds]

        push!(components, :plotdir => "$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])")
        push!(components, :dsdi    => (ds,integrate(ds)))
        push!(components, :update  => button("update"))

    end

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
