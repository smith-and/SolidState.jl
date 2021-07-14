module Apps

using Distributed
using OrderedCollections, LinearAlgebra, Distributed, PackageCompiler, BSON, SharedArrays
using Plots, BenchmarkTools, LsqFit, Dates, Measures
using Mux, WebIO, Interact, InteractiveUtils

using ..SolidState
using ..SolidState: cθ, CommensurateASD,  integrate, cointegrate, data_import

###########################
#### State Projector
###########################


function plot_callback(comps)

    k = comps[:k_input][]

    bz_plot = bz_hs_points(comps[:asdg], comps[:kpts], k)

    comps[:hd](k)
    site_weights = collect(eachslice(abs.(eigvecs(Hermitian(comps[:hd].h_ops.h))),dims=2))

    uc_xy = uc_face_plot(comps[:asdg])

    uc_plots = face_plots(comps[:asdg]; weights = site_weights[comps[:state_index][]])

    path_plot = SolidStateApps.band_broken_plot(comps[:bp_dict],10;
#                 size= (600,600)
        );
    #(1,comps[:state_index][])
    state_coords = (
        [comps[:k_odi_leg][comps[:k_input][]]],
        [comps[:energies][comps[:k_idx_leg][comps[:k_input][]]][comps[:state_index][]]]
    )
    scatter!(path_plot, state_coords..., color=:red, label="")

    l = @layout [ [ [ a b ] ; c ] [ d ; e ; f;] ]

    plot(
        bz_plot,uc_xy,path_plot,uc_plots[:xy],uc_plots[:yz],uc_plots[:xz],
        layout = l
    )
end

function state_projector_widget(asd,hd)

    asdg = asd|>SolidState.ASDGeometry
    bp_dict = SolidState.band_plot(asd,hd)

    comps = Dict{Symbol,Any}()
    comps[:hd] = hd
    comps[:asdg] = asdg
    comps[:bp_dict] = bp_dict

    xlabels = bp_dict[:args][:xticks][2]
    xticks = bp_dict[:args][:xticks][1]

    k_ledger = OrderedDict(getindex.(Ref(comps[:asdg]["bz_hs"]),comps[:bp_dict][:path]).=>comps[:bp_dict][:path])
    kpts = OrderedDict(comps[:bp_dict][:path].=>getindex.(Ref(comps[:asdg]["bz_hs"]),comps[:bp_dict][:path]))
    k_odi_leg = OrderedDict(getindex.(Ref(kpts),xlabels).=>xticks);
    k_idx_leg = OrderedDict(getindex.(Ref(kpts),xlabels).=>1:length(xticks));

    hs_odimeter_pts = OrderedDict([xlabels[i]=>findfirst(z->z==pt,comps[:bp_dict][:odimeter][1]) for (i,pt) ∈ enumerate(xticks) ]...)
    energies = collect(eachslice(comps[:bp_dict][:Es][values(hs_odimeter_pts)|>collect,:],dims=1))

    # Form the widget
    comps[:kpts] = kpts
    comps[:k_idx_leg] = k_idx_leg
    comps[:k_odi_leg] = k_odi_leg
    comps[:energies] = energies

    push!(comps, :update      => button("update") )
    push!(comps, :state_index => dropdown(max(1,Int(floor(size(energies[1],1)/2)-10)):min(Int(floor(size(energies[1],1))/2+10),size(energies[1],1)); width =10))
    push!(comps, :k_input     => dropdown(kpts) )

    on(comps[:state_index]) do x
        comps[:update][]+=1
    end

    on(comps[:k_input]) do x
        comps[:update][]+=1
    end

    output = map(comps[:update]) do args
        plot_callback(comps)
    end

    w = Widget(comps, output = output)

    @layout! w vbox(hbox(:state_index,:k_input),observe(_))
end

function state_projector_widget(;asd,mn,cachedir, kargs...)
    asd0 = BSON.load("$cachedir/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")
    state_projector_widget(asd0,hd)
end

function state_projector_widget(args)
    state_projector_widget(; args...)
end



###########################
#### Abs Angle Explorer
###########################

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

###########################
#### Snipe Machine
###########################

function snipe_input(executables)

    w = Dict{Symbol,Any}(
        :update => button("update"),
        :command   => dropdown(["run","walk","slurp"],value="slurp"),
        :comp_type => dropdown(["jit","sys"]),
        :scripts   => dropdown(executables),
        :workers => dropdown(string.([:big,collect(1:100)...])),
        :script_call => textbox(value = "\'main(some code goes here...)\'"),
        :slurm_args  => textbox(value = ""),
        :run         => button("run")
    )

    on(w[:workers]) do i
        w[:update][]+=1
    end

    on(w[:comp_type]) do i
        w[:update][]+=1
    end

    on(w[:script_call]) do i
        w[:update][]+=1
    end

    on(w[:slurm_args]) do i
        w[:update][]+=1
    end

    on(w[:command]) do i
        w[:update][]+=1
    end

    on(w[:scripts]) do dir
        w[:update][]+=1
    end

    on(w[:run]) do i
        run(`sh -c $(w[])`)
    end

    output = map(w[:update]) do i
        string((getindex.(getindex.(Ref(w),[:command,:comp_type,:scripts,:workers,:script_call,:slurm_args])).*" ")...)
    end

    w = Widget(w,output=output)

    @layout! w vbox(hbox(:run,:command,:comp_type,:scripts,:workers), :script_call, :slurm_args, observe(_))
end


function snipe_machine()

    w = OrderedDict{Symbol,Any}(
        :machines => dropdown(["local","dirac","bridges"])
    )

    output = map(w[:machines]) do dir
         scripts = filter(readdir("/home/smitha/Dropbox/Graduate/packages/.scripts/$(dir)")) do file
            file[1]=='.' && file!=".cache" && file!=".side"
        end|> x -> isempty(x) ? [""] : x

        snipe_input(scripts)
    end

    w = Widget(w,output=output)

    @layout! w vbox(:machines,observe(_))
end


###########################
#### Extraction Plotter
###########################
function ext_args(asd,RN,f,brng;)
   #Extraction Parameters

   #Collection Importing and Post-Process Path
   args = (
      outdir = mkpath("/home/smitha/Dropbox/Graduate/packages/.extract/bridges2/$(asd)/.plot/$(RN)/Ext"),
      local_dir  = "/home/smitha/Dropbox/Graduate/packages/.extract/bridges2",
      mount_dir  = "/ocean/projects/phy190028p/asmithc/scripts",
      asd = asd,
      RN  = RN,
      f   = f,
      brng= brng,
      dtype = "Ext",
   )

end

function chart_main_widget(asd,RN,f,brng; plotargs...)
   #Extraction Parameters
   kargs = (
      local_dir  = "/home/smitha/Dropbox/Graduate/completion/bridges2",
      mount_dir  = "/ocean/projects/phy190028p/asmithc/scripts",
      asd = asd,
      RN  = RN,
      f   = f,
      brng= brng
   )

   #Collection Importing and Post-Process Path
   args = (
      outdir = mkpath("$(kargs.local_dir)/$(kargs.asd)/.plot/$(kargs.RN)/Ext"),
      chnl = import_collection(; kargs...),
      dtype = "Ext",
      asd = asd,
      plotargs ...
   )

   #Run Plot Methods on Data
    pdict = Dict{Symbol,Any}(
        :heatmapf=>heatmapf(      kargs.f; args...),
       :cutheatmapf=>cutheatmapf(      kargs.f; xlims=(1.0,3.1), args...),
       :waterfall=>waterfall(     kargs.f; args...),
       :isowaterfall=>isowaterfall(  kargs.f; args...),
       :overlay=>overlay(       kargs.f; args...),
       :scaled_overlay=>scaled_overlay(kargs.f; args...),
       :angle_dep=>angle_dep(     kargs.f, kargs.brng;  args...),
    )

    comps = Dict(
        :selector => dropdown(collect(keys(pdict))),
        :args => args
        )

    output = map(comps[:selector]) do plt
        pdict[plt]
    end

    w = Widget(comps,output=output)

    @layout! w vbox(:selector, observe(_))

    w
end

function extraction_widget(kargs...)
    comps = Dict(
        :extract => button("extract"),
        :asd     => dropdown([:ASD2,:ASD4,:ASD2B]),
        :dtype   => textarea("Run"; value = "SHG-50k", cols =10, rows=1),
        :f       => dropdown([abs,real,angle,imag]),
        :brng    => dropdown([(0.0,1.0),(1.0,3.1),(0.0,4.0)])
    )

    output = map(comps[:extract]) do i
        args = ext_args(Symbol(comps[:asd][]),comps[:dtype][],comps[:f][],comps[:brng][])
        SolidStateApps.slide_wrapper(
            chart_main_widget(Symbol(comps[:asd][]),comps[:dtype][],comps[:f][],comps[:brng][])
            ; args...)
    end

    w = Widget(comps, output=output)

    @layout! w vbox(hbox(:extract,:asd,:f,:brng,:dtype), observe(_))

end


###########################
#### Slide Wrapper
###########################
# Add Save Button Component
function add_plot_save!(d,args)
    push!(d[:buttons],:savebutton)

    push!(d,:savebutton  => button("Save"))
    push!(d,:save_path   => textbox(; value = "plot"))

    p_path = pwd()

    try
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])")
    catch
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])")
    end

    on(d[:savebutton]) do b
        Plots.pdf(d[:w0][],"$(p_path)/"*d[:save_path].output[]*"-$b")
        Plots.png(d[:w0][],"$(p_path)/"*d[:save_path].output[]*"-$b")
        d[:update][]+=1
    end

    nothing
end

# Add Figure Tex Component
function add_plot_figure!(d,args)
    push!(d,:fig_caption => Widgets.textarea("figure caption";rows=5,cols=30))

    p_path = pwd()
    try
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])")
    catch
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])")
    end

    push!(d,:fig_tex => Widgets.textarea("caption"))
    on(d[:fig_caption]) do b
        d[:fig_tex].output[] = "\\begin{figure}
            \\includegraphics[width=1\\textwidth]{$(p_path)/$(d[:save_path].output[])-$(d[:savebutton][])}
            \\caption{\\label{test} $(d[:fig_caption].output[])}
        \\end{figure}"
        d[:update][]+=1
    end

    nothing
end

# Add Beamer Frame Tex Component
function add_fig_slide!(d)
    push!(d,:slide_text => Widgets.textarea("slide text";rows=5,cols=30))
    push!(d,:slide_tex  => Widgets.textarea("tex output", rows=10,cols=30))
    on(d[:slide_text]) do b
        d[:slide_tex].output[] = "
\\begin{frame}
    \\begin{columns}
        \\begin{column}{0.5\\textwidth}
            $(d[:slide_text].output[])
        \\end{column}
        \\begin{column}{0.5\\textwidth}
            $(d[:fig_tex].output[])
        \\end{column}
    \\end{columns}
\\end{frame}
        "
        d[:update][]+=1
    end

    nothing
end

#Add Callback and return a Widget

#Slide Creation tools
function slide_wrapper(w0; args...)
    # Creates a set of components which need to be
    components = OrderedDict{Symbol,Any}(
            :buttons => Symbol[],
            :update  => button("update"),
            :w0      => w0
            )

    add_plot_save!(components, args)
    add_plot_figure!(components, args)
    add_fig_slide!(components)

    output = map(components[:update]) do x
        components[:w0]
    end

    w = Widget(components, output=output)


    # Set layout of widget
    @layout! w hbox(
            observe(_),
            vbox(
                hbox(:savebutton,:save_path),
                :fig_caption,
                :slide_text
            ),
            :slide_tex
        )


    w
end

###########################
#### Component Deck
###########################
function component_deck()
    comps = Dict{Symbol,Any}(
            :update           => button("update"),
            :add_tab          => button("add tab"),
            :clear_tab        => button("clear tab"),
            :add_row          => button("add row"),
            :clear_component  => button("clear"),
            :export           => button("export"),
            :widgetf          => dropdown([abs_angle_selector,state_projector_widget,extraction_widget]),
            :components       => OrderedDict{Int64,Any}()
    )

    on(comps[:add_tab]) do i
        push!(comps[:components], i => vbox())
        comps[:tabs][] = length(comps[:components]|>keys)
        comps[:update][]+=1
    end

    # on(comps[:inputs]) do i
    #     comps[:update][]+=1
    # end

    on(comps[:clear_component]) do i
        empty!(comps[:components])
        comps[:add_tab][]=1
        comps[:update][]+=1
    end

    comps[:tabs]=tabulator(comps[:components])

    on(comps[:clear_tab]) do i
        tabN = comps[:tabs][]
        comps[:components][tabN] = vbox()
        comps[:update][]+=1
    end

    on(comps[:add_row]) do i
        tabN = comps[:tabs][]
        if comps[:inputs][][:widgetf] == extraction_widget
            comps[:components][tabN] = vbox(comps[:components][tabN],comps[:inputs][][:widgetf](comps[:inputs][]))
        else
            comps[:components][tabN] = vbox(comps[:components][tabN],comps[:widgetf][]())
        end
        comps[:update][]+=1
    end

    on(comps[:export]) do i

    end

    output = map(comps[:update]) do i
        tabN = comps[:tabs][]
        comps[:tabs]=tabulator(comps[:components])
        comps[:tabs][] = tabN

        comps[:tabs]
    end

    w = Widget(comps, output=output)

    @layout! w vbox(hbox(:add_tab,:clear_tab,:add_row,:clear_component,:widgetf),observe(_))

    w
end

###########################################
#### Serving Methods
###########################################

### Widget Composition
"""
    component_map(inputs)

Pass an OrderedDict containing dictionaries which are as widget descriptions and can be called as

    inputs[i][:widget](intputs[i])

The widgets are then displayed vertically in the order they are listed, in one vbox
"""
function component_map(inputs)
    components = Dict{Symbol,Any}()
    for key in keys(inputs)
        push!(components, key => inputs[key][:widgetf](inputs[key]))
    end

    components
end

### Mux Utilities.
### Currently really only at a sinlge page

"""
    appserve(app::Function)

Starts a firefox tab with a single page webio server displaying the output of app
"""
function appserve(app::Function)
    port = rand(8000:9000)
    webio_serve(page("/", app), port)
    run(`firefox localhost:$port`)
end



"""
    widget_serve(widget)

Serve widget in firefox tab on single page
"""
function widget_serve(widget)
    SolidStateApps.appserve() do req
        widget
    end
end

end
