function side_view_uc(asdg,i1,i2; weights = :none, colors = :auto, kargs...)
    face_symbol = Dict(1=>:x, 2=>:y, 3=>:z)
    box_size= max(max(asdg["xtal"][2]...)...)*1.2
    if colors == :auto
        site_colors = asdg["colors"]
    else
        site_colors = fill(:black,length(asdg["xtal"][2]));
    end

    if weights == :none
        site_weights = fill(1.0,length(asdg["xtal"][2]))
    else
        site_weights = weights
    end

    scatter(getindex.(asdg["xtal"][2],i1),getindex.(asdg["xtal"][2],i2);
        frame = :box,
        ylims = (-box_size,box_size),
        xlims = (-box_size,box_size),
        label = "",
        xguide = "$(face_symbol[i1]) (nm)",
        yguide = "$(face_symbol[i2]) (nm)",
#         bottom_margin = 5Plots.mm,
#         top_margin = 5Plots.mm,
        yrotation = 60,
        xrotation = 20,
        color = site_colors,
        opacity = site_weights,
        size = (300,300),
        m = 5,
        kargs...
    )
end

function face_plots(asdg; kargs...)
    plots = Dict(
        :xy => side_view_uc(asdg,1,2; kargs...),
        :yz => side_view_uc(asdg,2,3; kargs...),
        :xz => side_view_uc(asdg,1,3; kargs...),
    )
end

function face_plot(asd; kargs...)
    asdg = SolidState.ASDGeometry(asd)
    plots = Dict(
        :xy => side_view_uc(asdg,1,2; kargs...),
        :yz => side_view_uc(asdg,2,3; kargs...),
        :xz => side_view_uc(asdg,1,3; kargs...),
    )
end

function bz_hs_points(asdg, pts, k_selected)
    k_scale = max(max(values(asdg["bz_c"])...)...)
    bz_plt = plot(collect(eachslice(hcat(asdg["bz_c"]...,asdg["bz_c"][1]),dims=1))[1:2]...;
        aspectratio=1.0,
        color=:black,
        label="",
        frame=:none,
#         top_margin = 5Plots.mm,
#         right_margin = 5Plots.mm,
        size= (300,300)
    )

    map(keys(pts),values(pts)) do k,v
        annotate!((v[1:2].+0.1.*k_scale.*(1,1))...,k)
    end

    map(keys(pts),values(pts)) do k,v
        if ((v-k_selected)|>norm) < 1e-8
            scatter!(bz_plt,[v[1]],[v[2]], label="", color=:red)
        else
            scatter!(bz_plt,[v[1]],[v[2]], label="", color=:black)
        end
    end
    bz_plt
end


function uc_face_plot(asdg)
    uc_xy = side_view_uc(asdg,1,2,
        frame = :none,
        size= (300,300)
#         bottom_margin = 1Plots.mm,
#         top_margin = 1Plots.mm,
    )

    plot!(uc_xy, collect(eachslice(hcat(asdg["uc_c"]...,asdg["uc_c"][1]),dims=1))[1:2]...;
        aspectratio=1.0,
        color=:black,
        label="",
        frame=:none,
        left_margin = 5Plots.mm,
        size= (300,300)
    )


end

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
    bp_dict = SolidStateApps.band_plot(asd,hd)

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


### Commensurate Arg Plot
function hull_comargs_plot(dcut)

        #need way to relate these to the dcut
        (mm,sm) = (200,200)
        θid = θindex(mm,sm)
        θca = θcomargs(mm,sm)
        θsp = θspectra(mm,sm)
        θLM = LMspectra(mm,sm)
        θse = θseries(mm,sm)

        unqθ = (union(round.(θsp[θLM .< dcut ],digits=4)))
        θspargs = [findall(x->round(x,digits=4)==θ,θsp[θLM .< dcut ]) for θ∈unqθ]
        θspcopies=getindex.(Ref(θsp[θLM .< dcut ]),θspargs)
        θLMcopies=getindex.(Ref(θLM[θLM .< dcut ]),θspargs)
        θcacopies=getindex.(Ref(θca[θLM .< dcut ]),θspargs)

        θhullca = Vector{typeof(θca[1])}(undef,length(unqθ))
        θhullsp = Vector{typeof(θsp[1])}(undef,length(unqθ))
        θhullLM = Vector{typeof(θLM[1])}(undef,length(unqθ))
        for (i,LMs) ∈ enumerate(θLMcopies)
                θhullca[i] = θcacopies[i][argmin(LMs)]
                θhullsp[i] = θspcopies[i][argmin(LMs)]
                θhullLM[i] = θLMcopies[i][argmin(LMs)]
        end


        round.(sort(θhullsp),digits=4)

        hullpatch = scatter(θsp,θLM, group=θse, legendloc=:topleft,label="",ylabel="",opacity=0.75, ylim=[-dcut/10,2*dcut])
        plot!(hullpatch,[0,60],[dcut,dcut],label="dimension mulitplier",color=:black)
        Plots.scatter!(hullpatch, θhullsp,θhullLM, legendloc=:topleft,color=:black, label="")
        scatter!(hullpatch,θhullsp,[-(dcut/10).+0.0.*θhullsp], color=:black, label="", opacity=0.5)

        return hullpatch
end
