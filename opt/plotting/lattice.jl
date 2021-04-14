#TO DO
#Plot the structure factors as weights in a sum of Gaussian densities the unit cell sites, or just with opacity and a fixed glyph
#Lattice Graphics
#plot outline of Wigner Seitz Unit Cell
function  add_uc_outline(uc_plot, asdg, a=1, b=2, α=1,θ=0, color=:black)
    ucc = [R3D(α,θ)*v for v∈ asdg["uc_c"]]
    uc_c_a = [getindex.(ucc,Ref(a))...] |> x->push!(x,x[1])
    uc_c_b = [getindex.(ucc,Ref(b))...] |> x->push!(x,x[1])

    plot!(uc_plot,uc_c_a,uc_c_b, color=color, leg=false, label=false)
end

function  add_bz_outline(bz_plot, asdg, a=1, b=2, α=1,θ=0, color=:black)
    bzc = [R3D(α,θ)*v for v∈ asdg["bz_c"]]
    bz_c_a = [getindex.(bzc,Ref(a))...] |> x->push!(x,x[1])
    bz_c_b = [getindex.(bzc,Ref(b))...] |> x->push!(x,x[1])

    plot!(bz_plot,bz_c_a,bz_c_b, color=color, leg=false, label=false)
end

function uc_label_plot(asdg)
    #Form a plot and add uc outline
    uc_plot = plot( framestyle=:none, grid=:no, axes=false, leg=false, aspectratio=1)

    #plot each atom in the site list of the ASD
    row_pos             = [2*(i-(length(asdg["xtal"][2])/2+1)) for i=1:length(asdg["xtal"][2])]
    marker_col_pos      = [-1 for _=1:length(asdg["colors"])]
    marker_annotations  = [(1, row_pos[i], Plots.text(asdg["labels"][i], 6)) for i=1:length(asdg["colors"])]

    scatter!(uc_plot, marker_col_pos, row_pos, color=asdg["colors"], marker=asdg["glyph"], opacity=.75, leg=false, annotations = marker_annotations)
    xlims!(-1.75, 1.2*xlims()[2])
    ylims!(ylims().*1.2)

    uc_plot
end

function uc_side_plot(asdg, a, b, title, ws_outline=true)
    #Form a plot and add uc outline
    #marker_annotations  = [(1, row_pos[i], Plots.text(asdg["labels"][i], 6)) for i=1:length(asdg["colors"])]
    #scatter!(uc_plot, marker_col_pos, row_pos, color=asdg["colors"], marker=asdg["glyph"], opacity=.75, leg=false, annotations = marker_annotations)

    uc_plot = plot(framestyle=:none, grid=:no, title=title, titlefontsize = 8, axes=false, leg=false, titlefontvalign=:bottom, title_location=:right    )
    ws_outline ? add_uc_outline(uc_plot,asdg,a,b) : nothing

    #plot each atom in the site list of the ASD
    atom_a_pos = getindex.(asdg["xtal"][2],Ref(a))
    atom_b_pos = getindex.(asdg["xtal"][2],Ref(b))

    scatter!(uc_plot, atom_a_pos,atom_b_pos, color=asdg["colors"], marker=asdg["glyph"], opacity=.75)

    if     a==1 && b==2
        xlims!((xlims().-sum(xlims())/2) .* 1.2 .+ sum(xlims())/2)
        ylims!((ylims().-sum(ylims())/2) .* 1.2 .+ sum(ylims())/2)
    elseif a==2 && b==3
        xlims!((xlims().-sum(xlims())/2) .* 1.2 .+ sum(xlims())/2)
        ylims!((ylims().-sum(ylims())/2) .* 1.2 .+ sum(ylims())/2)
    elseif a==1 && b==3
        xlims!((xlims().-sum(xlims())/2) .* 1.2 .+ sum(xlims())/2)
        ylims!((ylims().-sum(ylims())/2) .* 1.2 .+ sum(ylims())/2)

    end
    uc_plot
end

function uc_plot_panels(ex)
    asdg = ASDGeometry(ex["asd"])
    plot_names = [:XY,:YZ,:XZ,:labels]
    plots = [
        uc_side_plot(asdg, 1,2, "XY", true),
        uc_side_plot(asdg, 2,3, "YZ", true),
        uc_side_plot(asdg, 1,3, "XZ", true),
        uc_label_plot(asdg)
    ]

    Dict{Symbol,AbstractPlot}(Pair.(plot_names,plots))
end

function uc_XY_plot(ex)
    #Form a plot and add uc outline
    #marker_annotations  = [(1, row_pos[i], Plots.text(asdg["labels"][i], 6)) for i=1:length(asdg["colors"])]
    #scatter!(uc_plot, marker_col_pos, row_pos, color=asdg["colors"], marker=asdg["glyph"], opacity=.75, leg=false, annotations = marker_annotations)
    asdg0 = eval(quote $(ex["model_input"].asd)() end )|>ASDGeometry
    asdg  = ASDGeometry(ex["asd"])
    ϕ=((ex["model_input"].s!=0) ? Lattice.cθ(ex["model_input"].m,ex["model_input"].m+ex["model_input"].s) : 0.0)
    uc_plot = plot(framestyle=:none, grid=:no, axes=false, leg=false,aspectratio=1)
    add_uc_outline(uc_plot,asdg,1,2)
    #add_uc_outline(uc_plot,asdg0,1,2)
    add_uc_outline(uc_plot,asdg0,1,2,1,ϕ/2,:red)
    add_uc_outline(uc_plot,asdg0,1,2,1,-ϕ/2,:purple)

    #plot each atom in the site list of the ASD
    atom_a_pos = getindex.(asdg["xtal"][2],Ref(1))
    atom_b_pos = getindex.(asdg["xtal"][2],Ref(2))

    scatter!(uc_plot, atom_a_pos,atom_b_pos, color=asdg["colors"], marker=asdg["glyph"], opacity=.75)

    xlims!((xlims().-sum(xlims())/2) .* 1.2 .+ sum(xlims())/2)

    uc_plot
end

function bz_XY_plot(ex)
    #Form a plot and add bz outline
    #marker_annotations  = [(1, row_pos[i], Plots.text(asdg["labels"][i], 6)) for i=1:length(asdg["colors"])]
    #scatter!(bz_plot, marker_col_pos, row_pos, color=asdg["colors"], marker=asdg["glyph"], opacity=.75, leg=false, annotations = marker_annotations)
    asdg0 = eval(quote $(ex["model_input"].asd)() end )|>ASDGeometry
    asdg  = ASDGeometry(ex["asd"])
    ϕ=((ex["model_input"].s!=0) ? Lattice.cθ(ex["model_input"].m,ex["model_input"].m+ex["model_input"].s) : 0.0)

    bz_plot = plot(framestyle=:none, grid=:no, axes=false, leg=false )
    add_bz_outline(bz_plot,asdg,1,2)
    #add_bz_outline(bz_plot,asdg0,1,2)
    add_bz_outline(bz_plot,asdg0,1,2,1,ϕ/2,:red)
    add_bz_outline(bz_plot,asdg0,1,2,1,-ϕ/2,:purple)

    #plot each atom in the site list of the ASD
    #atom_a_pos = getindex.(asdg["xtal"][2],Ref(1))
    #atom_b_pos = getindex.(asdg["xtal"][2],Ref(2))
    #scatter!(bz_plot, atom_a_pos,atom_b_pos, color=asdg["colors"], marker=asdg["glyph"], opacity=.75)

    #xlims!((xlims().-sum(xlims())/2) .* 1.2 .+ sum(xlims())/2)

    bz_plot
end


function lattice_panels(ex)
    plot_names = [:uc,:bz]
    plots = [
        uc_XY_plot(ex),
        bz_XY_plot(ex),
    ]

    Dict{Symbol,AbstractPlot}(Pair.(plot_names,plots))
end

function hull_comargs_plot(dcut)


        #need way to relate these to the dcut
        (mm,sm) = (200,200)
        θid = Lattice.θindex(mm,sm)
        θca = Lattice.θcomargs(mm,sm)
        θsp = Lattice.θspectra(mm,sm)
        θLM = Lattice.LMspectra(mm,sm)
        θse = Lattice.θseries(mm,sm)

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
        plot!(hullpatch,[0,60],[dcut,dcut],label="dimension cut",color=:black)
        Plots.scatter!(hullpatch, θhullsp,θhullLM, legendloc=:topleft,color=:black, label="")
        scatter!(hullpatch,θhullsp,[-(dcut/10).+0.0.*θhullsp], color=:black, label="", opacity=0.5)
        return hullpatch
end


#=
#Path Plotting
function bz_plot0(ex)
    #Start by composing the path plot
    plt= plot(yrotation=60, xtickfonthalign=:center, xtickfontsize=6, yguide ="E(eV)", yguidefontsize=6)# left_margin=1cm, right_margin=2cm, bottom_margin=1cm,top_margin=1cm)#, inset_subplots=(1, bbox(0, 0, 0.1, 0.1, :top,:left)))
    plot!(title=title, title_pos=:right, titlefontvalign=:center)
    plot!(    minorgrid=false, framestyle=:box, grid=:x, gridalpha=.7 , gridlinewidth=1, gridstyle=:dashdot)
    xticks!(  tick_odimeters, tick_symbols)
    xlims!(   odimeter[1],odimeter[end])
    ylims!(min(band_data...),max(band_data...))
    #ylim!=:all && ylims!(ylim...)
    plot!(    odimeter, transpose(band_data), label=false)
    #=#println("odimeter\n", odimeter,"band_data\n", transpose(band_data))
    #Make an inset plot and then add series data
    plot!(plt, size=size, inset_subplots = [(1, bbox(0.01,0.01,0.3w,0.3w,:bottom,:right))],  bg_color_inside=RGBA(1,1,1,0), framestyle=:none, grid=:no, axes=false, subplot=2)
    plot!(plt[2], collect(eachslice(hcat(bz_c...);dims=1))..., color=:black, label=false)
    plot!(plt[2], collect(eachslice(hcat(corners...);dims=1))..., line=(:red), label=false)
    =#
    plt
end


function section_function(type::Val{:lattice},  program::Dict{String,Any}, ex_idx::Int64, scope_idx::Int64=1; kwargs...)
    ex = program["executables"][ex_idx];

    uc_plots = Analysis.uc_plot_panels(ex)
    l = @layout [ a b; c d ]
    uc_face_plots = plot(uc_plots[:XY], uc_plots[:XZ], uc_plots[:YZ], uc_plots[:labels],
                    layout = l, framestyle=:none, grid=:no, axes=false, leg=false)
    merge!(uc_plots,Dict(:faces=>uc_face_plots))
end

function section_function(type::Val{:shiftplot},  programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; kwargs...)
    ex = programs["executables"][ex_idx]
    asdb = ASDBasics(ex["asd"])
    asdg = ASDGeometry(ex["asd"])

    uc_c  = getindex.(asdg["uc_c"],Ref(1:2))
    uc_c  = push!([uc_c...],uc_c[1])

    l_colors = map(l_idx-> l_idx==1 ? :blue : :red, asdg["layer"])
    vr=vrd(10)
    lbase = LatticePointsSym((transpose(asdg["xtal"][1]),([0, 0, 0 ],)),4*[0 1 0 ;-1 0 0 ; 0 0 0 ],3)[1]
    println("did shift plot for cover size ",length(lbase))
    s_plot = let lbase2=lbase
        @animate for s_i=-50:150
            shift_plot = plot()
            plot!(shift_plot, framestyle=:none, grid=:no, aspectratio=1, axes=false)
            plot!(shift_plot, collect(eachslice(hcat(uc_c...);dims=1))..., color=:black, leg=false)
            for j = 1:length(lbase2)
                scatter!(shift_plot, collect(eachslice(hcat(getindex.([v[3]==0.0 ? v+lbase2[j] : v+lbase2[j]+ sin((1/200)*2π*s_i)*(lbase2[3]+0.0*lbase2[2]) for v∈asdg["xtal"][2]],Ref(1:2))...);dims=1))..., color=l_colors, marker=asdg["glyph"], opacity=.5 )
            end
            xlims!(-LinearAlgebra.norm(lbase2[2]),LinearAlgebra.norm(lbase2[2]))
            ylims!(-LinearAlgebra.norm(lbase2[2]),LinearAlgebra.norm(lbase2[2]))
        end
    end
    gif(s_plot, "shift.gif", fps = 12)
end

=#
