module Plot

using Revise
using Distributed, Dates, OrderedCollections, BSON, Plots
using LinearAlgebra, SharedArrays, StaticArrays
using CubicSplines, Roots, SpecialFunctions, HCubature
using Base: Threads
using ..SolidState

######################################################################
### Crystal Geometry Plots
######################################################################

θindex(mmax,smax) = [(m,s) for m∈-mmax:mmax for s∈-smax:smax]
θcomargs(mmax,smax) = [(m:m,s:s)  for m∈-mmax:mmax for s∈-smax:smax]
θseries(mmax,smax) = [s  for m∈-mmax:mmax for s∈-smax:smax]
θspectra(mmax,smax) = [SolidState.cθ(m,m+s)*180/π  for m∈-mmax:mmax for s∈-smax:smax]
LMspectra(mmax,smax) = [(SolidState.TwistedTriangularGeometry(SolidState.ASD2()["blv"],(m,m+s))["blv"]|>det)/(SolidState.ASD2()["blv"]|>det)  for m∈-mmax:mmax for s∈-smax:smax]

### Commensurate Arg Plot
function bulkhead_indices(dcut)
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
    θsecopies=getindex.(Ref(θse[θLM .< dcut ]),θspargs)

    θhullca = Vector{typeof(θca[1])}(undef,length(unqθ))
    θhullsp = Vector{typeof(θsp[1])}(undef,length(unqθ))
    θhullLM = Vector{typeof(θLM[1])}(undef,length(unqθ))
    θhullse = Vector{typeof(θse[1])}(undef,length(unqθ))

    for (i,LMs) ∈ enumerate(θLMcopies)
        θhullca[i] = θcacopies[i][argmin(LMs)]
        θhullsp[i] = θspcopies[i][argmin(LMs)]
        θhullLM[i] = θLMcopies[i][argmin(LMs)]
        θhullse[i] = θsecopies[i][argmin(LMs)]
    end

    round.(sort(θhullsp),digits=4)
    #
    hullpatch = scatter(θsp,θLM,
        legendloc=:topleft,
        label="",
        ylabel="",
        opacity=0.25,
        # group=θse,
        color=:black,
        ylim=[-dcut/10,1.3*dcut],
        xlim=[0.0,120.0]
        )
    # unique indices
    scatter!(hullpatch, θhullsp,θhullLM,
        group = θhullse,
        legendloc=:topleft,
        xlim=[0.0,120.0],
        label=""
        )
    #cutoff line
    plot!(hullpatch,[0,120],[dcut,dcut],
        label="dim cut",
        color=:black,
        xlim=[0.0,120.0]
        )
    plot!(hullpatch,[120,120],[0,dcut],
        label="",
        color=:black,
        xlim=[0.0,120.0]
        )
    # projection
    scatter!(hullpatch,θhullsp,[-(dcut/10).+0.0.*θhullsp],
        color=:black, label="", opacity=0.5,
        xlim=[0.0,120.0]
    )

       return hullpatch
end

##################################################
#### Crystal Faces
##################################################
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
        color = site_colors,
        opacity = site_weights,
        kargs...
    )
    # size = (300,300),
    # m = 5,

    plot!(collect(eachslice(hcat(asdg["uc_c"]...,asdg["uc_c"][1]),dims=1))[[i1,i2]]...;
        color=:black,
        label="",
        frame=:none,
        aspectratio=1.0,
    )
    # left_margin = 5Plots.mm,
    # size= (300,300)

end

function crystalfaces(asd; kargs...)
    asdg = SolidState.ASDGeometry(asd())
    plots = Dict(
        :xy => side_view_uc(asdg,1,2; kargs...),
        :yz => side_view_uc(asdg,2,3; kargs...),
        :xz => side_view_uc(asdg,1,3; kargs...),
    )
end

##################################################
#### 3D Unit Cell Plot
##################################################

function unitcell3d(asd; kargs...)
    asdb = asd|>SolidState.ASDBasics
    asdg = asd|>SolidState.ASDGeometry
    xtal = asdb["xtal"]

    atom_colors = Dict(
        "B" => RGB(1.0,0.0,0.0),
        "N" => RGB(0.0,0.0,1.0)
    )

    site_colors = map(asdb["sld"]) do site
        atom_colors[site[1][1]]
    end

    points  = hcat(getindex.(asd["sites"],3)...)

    plt = plot(collect(eachslice(points,dims=1))...;
        frame = :grid,
        label = "",
        ticks = :none,
        seriestype=:scatter,
        camera = (45,45),
        aspectratio = 1,
        color = site_colors,
        xlims = (-2,2).*asdg["Lmag"],
        ylims = (-2,2).*asdg["Lmag"],
        zlims = (-2,2).*asdg["Lmag"],
        kargs...
        )

    plt
end


##################################################
#### Crystal Bond Line Plots
##################################################

struct BondLine{V <: AbstractArray}
    a::V
    b::V
end

function blshift!(bl::BondLine,v::AbstractArray)
    bl.a .+= v
    bl.b .+= v
    bl
end

function blshift(bl0::BondLine,v::AbstractArray)
    bl = deepcopy(bl0)
    bl.a .+= v
    bl.b .+= v
    bl
end

function blrotate!(bl::BondLine,θ::Float64)
    bl.a .= SolidState.R2D(1.0,θ)*bl.a
    bl.b .= SolidState.R2D(1.0,θ)*bl.b
    bl
end

function getline(bl::BondLine)
    eachslice(hcat(bl.a,bl.b)',dims=2)|>collect
end

function plotline!(plt::AbstractPlot, bl::BondLine; kargs...)
    plot!(plt, getline(bl)...;kargs...)
end

struct BondTree{BL <: BondLine, G<:AbstractArray}
    branches::Vector{BL}
    leafs::Vector{BL}
    generators::G
end

function treeshift(bt0,shift)
    bt = deepcopy(bt0)
    blshift!.(bt.branches,Ref(shift))
    blshift!.(bt.leafs,Ref(shift))
    bt
end

function treeshift!(bt,shift)
    blshift!.(bt.branches,Ref(shift))
    blshift!.(bt.leafs,Ref(shift))
    bt
end

function treerotate!(bt,θ)
    blrotate!.(bt.branches,θ)
    blrotate!.(bt.leafs,θ)
    bt
end

function treerotate(bt0,θ)
    bt = deepcopy(bt0)
    blrotate!.(bt.branches,θ)
    blrotate!.(bt.leafs,θ)
    bt
end

function treeplot!(plt,bt::BondTree; kargs...)
    plotline!.(Ref(plt),bt.branches;
        label="",
        aspectratio=1,
        color=:black,
        kargs...
        )
    plt
end

function grow(bt,nvec::AbstractArray)
    push!.(Ref(bt.branches),blshift.(bt.leafs,Ref(bt.generators*nvec)))
end

function grow(bt,N::Int)
    mesh = [ [n1,n2] for n1 in -N:N, n2 in -N:N]
    grow.(Ref(bt),mesh)
    bt
end

function basic_tree(asd,N)
    asd0 = asd()
    asdb = SolidState.ASDBasics(asd0)
    xtal = asdb["xtal"][1][1:2,1:2]
    slv = asdb["xtal"][2]
    slv[1][1:2]
    blength = (xtal[1,:]|>norm)/sqrt(3)/2

    linesA = [BondLine([0.0,0.0],[reim(-blength*exp(im*2π*n/3+im*π/2))...]) for n in 0:2]
    blshift!.(linesA,Ref(-slv[1][1:2]))
    linesB = [BondLine([0.0,0.0],[reim(blength*exp(im*2π*n/3+im*π/2))...]) for n in 0:2]
    blshift!.(linesB,Ref(slv[1][1:2]))
    lines = vcat(linesA,linesB)

    bt = BondTree(eltype(lines)[],lines,xtal')
    grow(bt,N)
    bt
end

function twist_plot(asd,θ,N)
    bt = basic_tree(asd,N)
    bt_1 = treerotate(bt,θ)
    bt_2 = treerotate(bt,-θ)
    plt = plot(
        frame=:none,
        aspectratio = 1,
        )
    treeplot!(plt,bt_1; color=:black);
    treeplot!(plt,bt_2; color=:red);
    plt
end

function shift_plot(asd,shifts,N)
    bt = basic_tree(asd,N)
    bt_1 = treeshift(bt,shifts[1])
    bt_2 = treeshift(bt,shifts[2])
    plt = plot(
        frame=:none,
        aspectratio = 1,
        )
    treeplot!(plt,bt_1; color=:black);
    treeplot!(plt,bt_2; color=:red);
    plt
end


function resize_box!(plt,box)
    plot!(plt;
        xlim = (-box,box),
        ylim = (-box,box),
        )
end

##################################################
#### Symmetric Unit Cell
##################################################

const SO3Defining = [ [ 0 0 0 ;
  0 0 -1 ;
  0 1 0 ], [  0 0 1 ;
   0 0 0 ;
  -1 0 0 ], [ 0 -1 0 ;
  1 0 0 ;
  0 0 0 ],
]

axisangle(θ,n) = exp(θ.*sum(n.*SO3Defining))

function group_generate(G0)
        G = deepcopy(G0)
        for g∈G0, h∈G
                if (Ref(h*g) .- G).|>norm.|>(x->x>1e-8)|>all
                        push!(G,h*g)
                end
        end
        if length(G)==length(G0)
                return G0
        else
                group_generate(G)
        end
end

function parallel_project(bt0,direction)
    bt = deepcopy(bt0)
    filter!(bt.branches) do branch
        bond = branch.a .- branch.b
        abs(acos(dot(direction,bond)/(norm(bond)*norm(direction)))) < 1e-5
    end
    bt
end

function branch_parallel_Q(branch,direction)
    bond = branch.a .- branch.b
    abs(acos(dot(direction,bond)/(norm(bond)*norm(direction)))) < 1e-5
end

function branch_vertical_Q(branch)
    (abs(branch.a[1] -  branch.b[1]) < 1e-5) && (abs(branch.a[2]-branch.b[2]) < 1e-5)
end

function plane_bonds(asd)
    asdg = asd|>SolidState.ASDGeometry
    boundaryQ = SolidState.BoundaryEdgeQ((asdg["uc_c"]...,asdg["uc_c"][1],asdg["uc_c"][2]))
    bt = BondTree(typeof(BondLine(rand(3),rand(3)))[],typeof(BondLine(rand(3),rand(3)))[],zeros(3,3))
    for (i,site1) in enumerate(getindex.(asd["sites"],3))
        for site2 in getindex.(asd["sites"],3)
            if boundaryQ((site1+site2)/2) && (abs(site1[3]-site2[3]) < 1e-5)
                push!(bt.branches,BondLine(copy(site1),copy(site2)))
            end
        end
    end
    bt
end

function vertical_bonds(asd)
    asdg = asd|>SolidState.ASDGeometry
    bt = BondTree(typeof(BondLine(rand(3),rand(3)))[],typeof(BondLine(rand(3),rand(3)))[],zeros(3,3))
    for (i,site1) in enumerate(getindex.(asd["sites"],3))
        for site2 in getindex.(asd["sites"],3)
            push!(bt.branches,BondLine(copy(site1),copy(site2)))
        end
    end
    btz = parallel_project(bt,[0.0,0.0,1.0])

    btz
end

##################################################
#### Brillouin Zone Plots
##################################################


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

##################################################
#### Plotting for High Symmetry Spectra
##################################################

function com_windows(comargs)
    angles  = Vector{Float64}(undef,length(comargs))
    for (i,mn) ∈ enumerate(comargs)
        angles[i] = SolidState.cθ(mn...)*180/π
    end

    comargs .= comargs[sortperm(angles)];
    angles .= angles[sortperm(angles)];

    leftwindows = angles .+  vcat([0.0],reverse(diff(reverse(angles))))./2
    rightwindows = angles .+ vcat(diff(angles),[0.0])./2

    windows = eachslice(hcat(leftwindows,rightwindows),dims=1)|>collect

    OrderedDict((angles.=>windows)...)
end

function com_windows_loose(comargs;α=0.1)
    angles  = Vector{Float64}(undef,length(comargs))
    for (i,mn) ∈ enumerate(comargs)
        angles[i] = SolidState.cθ(mn...)*180/π
    end

    comargs .= comargs[sortperm(angles)];
    angles .= angles[sortperm(angles)];

    leftwindows = (angles  .+ vcat([0.0],reverse(diff(reverse(angles))))./2)*(1+α)
    rightwindows = angles .+ vcat(diff(angles),[0.0])./2*(1-α)

    windows = eachslice(hcat(leftwindows,rightwindows),dims=1)|>collect

    OrderedDict((angles.=>windows)...)
end

function plot_twist_spectra_levels(datafile; α,  outdir, saving = true, kargs...)
    data      = BSON.load(datafile)
    ksymbols  = data[:ksymbols]
    kenergies = data[:spectra]
    comargs   = data[:comargs]
    windows   = com_windows_loose(comargs; α=α)

    #Initialize data structures for plots and spectra
    kplots = copy(OrderedDict("hi"=>plot()))|>empty
    for k ∈ ["Γ","K1"]
        push!(kplots,
            k=>plot(;
                label="",
                frame = :box,
                title = "$k twist spectrum",
                xguide = "twist angle (deg)",
                yguide = "energy (eV)",
                color=:black,
            )
        )
        for key ∈ kenergies[k]|>keys
            println("Plotting $k $key");flush(stdout)
            plot!(kplots[k],
                fill(windows[key],length(kenergies[k])),
                collect(eachslice(hcat(kenergies[k][key],kenergies[k][key]),dims=1));
                label="",
                color=:black,
                kargs...
            )
            display(kplots[k])
        end
    end
    saving && println("Saving Plots");flush(stdout)
    saving && Plots.pdf.(values(kplots),"$outdir/".*keys(kplots))
    saving && Plots.png.(values(kplots),"$outdir/".*keys(kplots))

    kplots|>values
end

function plot_twist_spectra_points(datafile; α,  outdir, saving = true, kargs...)
    data      = BSON.load(datafile)
    ksymbols  = data[:ksymbols]
    kenergies = data[:spectra]
    comargs   = data[:comargs]
    windows   = com_windows_loose(comargs; α=α)

    #Initialize data structures for plots and spectra
    kplots = copy(OrderedDict("hi"=>plot()))|>empty
    for k ∈ ["Γ","K1"]
        push!(kplots,
            k=>plot(;
                label="",
                frame = :box,
                title = "$k twist spectrum",
                xguide = "twist angle (deg)",
                yguide = "energy (eV)",
                color=:black,
            )
        )
        for key ∈ kenergies[k]|>keys
            scatter!(kplots[k],
            fill(key,length(kenergies[k][key])),kenergies[k][key];
            label="",
            color=:black,
            marker=:circle,
            kargs...
            )
        end
        display(kplots[k])
    end
    saving && println("Saving Plots");flush(stdout)
    saving && Plots.pdf.(values(kplots),"$outdir/".*keys(kplots))
    saving && Plots.png.(values(kplots),"$outdir/".*keys(kplots))

    kplots|>values
end

function spectral_density(energies,damping,kenergies)
    d = length(energies)
    density = Array{Float64,2}(undef,(length(energies),length(keys(kenergies))))
    for (i,en) in enumerate(energies)
        for (j,kens) in enumerate(collect(values(kenergies)))
            z= 0.0
            for ken in kens
                z += imag(1.0 / Complex(ken-en,damping))
            end
            density[i+d*(j-1)] = -z
        end
    end
    return density
end

function plot_twist_spectral_density(datafile; cachedir, damping = 0.02, N=10, outdir, saving = true, kargs...)
    data      = BSON.load(datafile)
    ksymbols  = data[:ksymbols]
    kenergies = data[:spectra]
    comargs   = data[:comargs]
    windows   = com_windows(comargs)

    ucvols = Vector{Float64}(undef,length(comargs))
    for (i,mn)∈ enumerate(comargs)
        ucvols[i] = det(BSON.load("$cachedir/$(data[:asd])/asd-$(mn[1])-$(mn[2]).bson")["blv"][1:2,1:2])
    end

    #Initialize data structures for plots and spectra
    kplots = copy(OrderedDict("hi"=>plot()))|>empty
    for k ∈ ["Γ","K1"]

        emin = min([min(ens...) for ens in values(kenergies[k])]...)
        emax = max([max(ens...) for ens in values(kenergies[k])]...)

        energies = emin:((emax-emin)/N):emax
        angles = keys(kenergies[k])|>collect
        density = spectral_density(energies,damping,kenergies[k])
        spa     = sortperm(angles)
        vols    = ucvols[spa]


        println("Plotting $k");flush(stdout)
        kplots[k] = heatmap(angles[spa],energies,(density./(ucvols'))[:,spa];
            title = "$k twist spectrum",
            xguide = "twist angle (deg)",
            yguide = "energy (eV)",
            kargs...
            )
        display(kplots[k])
    end

    saving && println("Saving Plots");flush(stdout)
    saving && Plots.pdf.(values(kplots),"$outdir/".*keys(kplots))
    saving && Plots.png.(values(kplots),"$outdir/".*keys(kplots))

    kplots|>values
end

######################################################################
### Plotting Simple Data
######################################################################

##################################################
#### Plotting for State Projection
##################################################

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

function state_projector_widget(;asd,mn,cachedir, kargs...)
    ASD = BSON.load("$cachedir/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")

    asdg = ASD|>SolidState.ASDGeometry
    bp_dict = SolidState.band_plot(ASD,hd)

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

    plot_callback(comps)
end




##################################################
#### Plotting for Single Bandstructure
##################################################

"""
    band_region_plot(dict,(bmin,bmax); args...)
"""
function band_region_plot(dict,(bmin,bmax); args...)
    dim = size(dict[:Es],2)

    if bmin|>typeof <: Int
        if bmin > 0
            bm = bmin
        elseif bmin <= 0
            bm = Int(floor(dim/2))+bmin
        end
    elseif bmin == :half
        bm = Int(floor(dim/2))
    elseif bmin == :cond
        bm = Int(floor(dim/2))+1
        bx = (bmax == :end) ? dim : min(Int(floor(dim/2))+bmax,dim)
    elseif bmin == :val
        bx = Int(floor(dim/2))
        bm = (bmax == :end) ? 1 : max(Int(floor(dim/2))+1-bmax,1)
    end

    if all(bmin.!=(:val,:cond))
        if bmax|>typeof <: Int
            if bmax > 0
                bx = bmax
            elseif bmax <= 0
                bx = Int(floor(dim/2))+1-bmax
            end
        elseif bmax == :half
            bx = Int(floor(dim/2))+1
        elseif bmax == :end
            bx = dim
        end
    end
    bm = max(1,bm)
    bx = min(dim,bx)

    Emin = min(dict[:Es][:,bm:bx]...)|>x->(x < 0 ? 1.01*x : 0.99*x)
    Emax = max(dict[:Es][:,bm:bx]...)|>x->(x > 0 ? 1.01*x : 0.99*x)

    plt = Plots.plot(dict[:odimeter], dict[:Es] ; dict[:args]..., ylims=(Emin,Emax), opacity=0.3, args...)
    Plots.plot!(dict[:odimeter], dict[:Es][:,bm:bx] ; dict[:args]..., ylims=(Emin,Emax), args...)
end

function band_broken_plot(dict, nbands; args...)
    plot(
        band_region_plot(dict,(:cond,nbands),
            frame=:axis, xaxis = false, xtickfontsize=1,margin=0.0Plots.mm),
        band_region_plot(dict,(:val,nbands),
            frame=:axis,margin=0.0Plots.mm,title=""),
        layout = grid(2,1),
        args...
    )
end

function band_step_gif(dict; plotdir, fps=2.5, args...)
    anim = Animation()
    for nbands ∈ 1:Int(floor(size(dict[:Es],2)/2))
        plt = band_broken_plot(dict,nbands)
        #display(plt)
        frame(anim,plt)
        if (nbands==1)||(nbands==Int(floor(size(dict[:Es],2)/2)))
            for _∈1:10
                frame(anim,plt)
            end
        end
    end
    gif(anim,"$plotdir/band_step.gif",fps=fps)
end

function collection_gif(dicts::Vector{AbstractDict},NBs; plotdir, fps=2.5, args...)
    anims = Vector{typeof(Animation())}(undef, length(dicts))
    for i∈1:length(dicts)
        anims[i] = Animation()
    end

    plotbook0 = band_plotbook(dict[1],NBs,plotdir=plotdir)
    # plotseries = Vector{typeof(values(plotbook0)[1])}(undef,(length(dicts),length(values(plotbook0))))
    for (i,dict) ∈ enumerate(dicts)
        plotbook = band_plotbook(dict)
        for (j,plt) ∈ enumerate(values(plotbook))
            frame(anims[j],plt)
        end
    end

    for (i,anim) ∈ enumerate(anims)
        gif(anim,"$plotdir/$(keys(plotbook0)[i]).gif",fps=fps)
    end
end

"""
    band(; mn::Tuple{Int,Int}, NBs = [1,2], asd::Symbol, plotdir, scriptdir, cachedirr, pathlist=["K1","Γ","M1","K'1"], Npath = 300, plotting=true,args...)

Calculate the band structure of a model along the high symmetry points listed in
"""
function bands(dict, NBs; args...)
    plotbook = Dict(
        "band-reg-all"=>band_region_plot(       dict,(1,:end); args...),
        "band-val-all"=>band_region_plot(       dict,(:val,:end); args...),
        "band-cond-all"=>band_region_plot(      dict,(:cond,:end); args...),
        "band-broken-all"=>band_broken_plot(    dict,:end; args...),
    )
    for NB ∈ NBs
        merge!(plotbook,Dict(
        "band-reg-$NB"=> band_region_plot(      dict,(-NB,-NB); args...),
        "band-val-$NB"=>band_region_plot(       dict,(:val,NB); args...),
        "band-cond-$NB"=>band_region_plot(      dict,(:cond,NB); args...),
        "band-broken-$NB"=>band_broken_plot(    dict,NB; args...),
        ))
    end
    plotbook
end

##################################################
#### Plotting for Data Integral
##################################################

function spectral_section(; plotdir, asd,mn,asdg,patch_rng,dω,re1,re2,Δ1,Δ2,kargs...)

    plts = Dict{Symbol,typeof(plot())}()

    plts[:dω] = heatmap(patch_rng,patch_rng,abs.(dω)',
        color = :bilbao,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        # title = "dω"
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )


    plts[:re1] = heatmap(patch_rng,patch_rng,abs.(re1)',
        color = :bilbao,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        # title = "rₑ₁"
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )

    ###############

    plts[:re2] = heatmap(patch_rng,patch_rng,abs.(re2)',
        color = :bilbao,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        # title = "rₑ₂"
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )

    ###############

    plts[:Δ1] = heatmap(patch_rng,patch_rng,abs.(Δ1)',
        color = :bilbao,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        # title = "Δ₁"
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )


    ###############

    plts[:Δ2] = heatmap(patch_rng,patch_rng,abs.(Δ2)',
        color = :bilbao,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        # title = "Δ₂"
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )

    #############################################
    plts[:panel]=plot(plts[:dω],plts[:Δ1],plts[:re1],plts[:re2],
        layout = (@layout grid(2,2)),
        margins = 1Plots.mm
        )

    SolidState.Main.book_save(plts,"$plotdir/$asd-$(mn[1])-$(mn[2])"|>mkpath)

    plts
end

#############################################
#### SHG Direct Sections
#############################################
function insetplot_example()
    p=plot()
    plot!(rand(100))
    bar!(rand(5), inset_subplots = [(1, bbox(0.7,0.025,0.3,0.3))], subplot=2, label="")
    bar!(p[2],-rand(5), label="")
end

function shg_section(args=(frame=:none,); plotdir, asd,mn,asdg,patch_rng,T1d,T2d,kargs...)

    θ = round(SolidState.cθ(mn...)*180/π,digits=2)

    plts = Dict{Symbol,typeof(plot())}()

    plts[:T1d] = heatmap(patch_rng,patch_rng,abs.(T1d)';
        color = :bilbao,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        title = "$(θ)ᵒ T1d",
        args...
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )

    plts[:T2d] = heatmap(patch_rng,patch_rng,abs.(T2d)';
        color = :bilbao,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        title = "$(θ)ᵒ T1d",
        args...
        # title = "rₑ₁"
        );

    # plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
    #     label="",
    #     color = :black,
    #     linestyle=:dot
    #     )

    plts[:T1d_phase] = heatmap(patch_rng,patch_rng,angle.(T1d)',
        color = :twilight,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        # title = "dω"
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )

    plts[:T2d_phase] = heatmap(patch_rng,patch_rng,angle.(T2d)',
        color = :twilight,
        aspectratio = 1,
        frame = :none,
        xlims = (patch_rng[1],patch_rng[end]),
        ylims = (patch_rng[1],patch_rng[end]),
        # title = "rₑ₁"
        );

    plot!([getindex.(asdg["bz_c"],1)...,asdg["bz_c"][1][1]],[getindex.(asdg["bz_c"],2)...,asdg["bz_c"][1][2]],
        label="",
        color = :black,
        linestyle=:dot
        )

    SolidState.Main.book_save(plts,"$plotdir/$asd-$(mn[1])-$(mn[2])"|>mkpath)

    plts
end

function shifted_shg(; RN, asd, Nevals, steps, kargs...)
    cache = BSON.load("$(ENV["scriptdir"])/out/$RN/shifted-$asd-$Nevals-$steps.bson")
    data = cache[:data]
    bdom = cache[:bdom]
    plt = plot()
    rmax = 0.0
    map(values(data)|>collect) do spectra
        plot!(plt,bdom,abs.(spectra),label="")
        rmax = max(rmax,abs.(spectra[1]))
    end

    lens!([0.0,0.5], [0.0,2rmax],
        inset = (1, bbox(0.7, 0.0, 0.3, 0.4)),
        yrot = 60
        )
    mkpath("$(ENV["scriptdir"])/plot/$RN")
    Plots.pdf(plt, "$(ENV["scriptdir"])/plot/$RN/shifted-$asd-$Nevals-$steps")
    plt
end

function shifted_shg(RN,asd,Nevals,steps, kargs...)
    shifted_shg(RN=RN,asd=asd,Nevals=Nevals,steps=steps)
end


##################################################
#### Plotting for Data Integral
##################################################

#Generic Plotting Method for Charts
function integral(f; evals, err, base, data, plotdir, handle, kargs...)
    plts = Vector{typeof(plot())}(undef,length(err))
    for i∈1:length(err)
        plts[i] = plot(base, f.(data[i][:]);
            ribbon = err,
            xlims  = (base[1],base[end]),
            legend = :topleft,
            label  = "$(evals[i])",
            frame  = :box,
            margins = 8Plots.mm,
            ylims = (min(0.0,(f.(data[i][:]))...),1.1*max(f.(data[i][:])...,1e-15)),
            )
        Plots.pdf(plts[i],"$plotdir/$handle-err-$i")
    end
    plts
end


#Generic Plotting Method for Charts
function integral!(plt::AbstractPlot, di::DataIntegral, f::Function, i::Int ; args...)
    plot!(plt,
        getindex.(getfield(di.dm.chart,1).base,1),
        f.(di.data[i][:]);
        ribbon = di.err[i],
        xlims  = (getindex.(getfield(di.dm.chart,1).base,1)[1],getindex.(getfield(di.dm.chart,1).base,1)[end]),
        legend = :topleft,
        label  = "$(di.evals[i])",
        frame  = :box,
        margins = 8Plots.mm,
        ylims = (min(0.0,(f.(di.data[i][:]))...),1.1*max(f.(di.data[i][:])...,1e-15)),
        args...
    )

    plt
end

function add_linecut!(plt_spectra,di,i,f;args...)
    base = getindex.(getfield(di.dm.chart,1).base,1)
    plot!(plt_spectra,[base[i],base[i]],[min(f.(di.data[di.evals.==min(di.evals...)][1])[:]...),max(f.(di.data[di.evals.==max(di.evals...)][1])[:]...)];
        label="",
        args...
        )

    plt_spectra
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


##################################################
#### Plotting for Data Section
##################################################

function add_hsp_labels!(sectionplot, asdg, pt_keys)
    # Add point labels
    M = inv(asdg["dxtal"][1])
    pts = Dict(pt_keys.=>[M*v for v ∈ getindex.(Ref(asdg["bz_hs"]),pt_keys)])

    map(keys(pts),values(pts)) do k,v
        scatter!(sectionplot,[v[1]],[v[2]],color = :black,label="")
        annotate!(sectionplot,[(v[1]-0.1,v[2]+0.1,Plots.text(k))],label="")
        plot!(sectionplot,[v[1]-0.02,v[1]-0.06],[v[2]+0.02,v[2]+0.06],label="",color=:black)
    end

    sectionplot
end

function section(ds::DataSection,f,idx,hsp=true;  args...)
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
        args...
    )
end

######################################################################
### Plotting Channeled Data
######################################################################

###################################
### Band Series Plots
###################################

function band_series(chnl; idxs)
    data  = chnl.data[idxs]

    angles = getindex.(data,:angle)
    odi = data[1][:data][:odimeter]
    cond_band_data = zeros(size(data[1][:data][:Es],1),length(data))
    val_band_data  = zeros(size(data[1][:data][:Es],1),length(data))

    valence_plt = plot()
    conduct_plt = plot()
    bandwidth_plt = plot(title = "log(bandwidth)", frame = :box)
    gap_plt = plot(title = "bandgap (eV)", frame = :box)

    for (i, datum) in enumerate(data)
        dim = size(datum[:data][:Es],2)/2|>Int
        odi = datum[:data][:odimeter][1]
        valence_band = datum[:data][:Es][:,dim]
        conduct_band = datum[:data][:Es][:,dim+1]

        cond_band_data[:,i] .= conduct_band./max(abs.(conduct_band)...)
        val_band_data[:,i] .= valence_band./valence_band[1]


        plot!(valence_plt, odi, valence_band.-valence_band[1],#.+datum[:angle]/50 ,
            label="",
            color=:black,
            )
        plot!(conduct_plt, odi, (conduct_band.-conduct_band[1])./max(abs.(conduct_band)...),#.+datum[:angle]/50 ,
            label="",
            color=:black,
            )

        scatter!(gap_plt, [datum[:angle]],[min(conduct_band...)-max(valence_band...)],
            color = :black,
            label = "",
            )

        scatter!(bandwidth_plt, [datum[:angle]],[log(max(conduct_band...)-min(conduct_band...))],
            color = :blue,
            label = "",
            )

        scatter!(bandwidth_plt, [datum[:angle]],[log(max(valence_band...)-min(valence_band...))],
            color = :red,
            label = "",
            )
    end

    cond_heatmap = heatmap(angles,odi,cond_band_data)
    vale_heatmap = heatmap(angles,odi,val_band_data)

    l = @layout [ a b ; c d ; e{.2h} f ]

    panel = plot(cond_heatmap, conduct_plt,
         vale_heatmap, valence_plt,
         bandwidth_plt, gap_plt,
        layout = l,
        size   = (800, 800)
        )

    Dict(
        :cond_heatmap => cond_heatmap,
        :conduct_plt => conduct_plt,
        :vale_heatmap => vale_heatmap,
        :valence_plt => valence_plt,
        :bandwidth => bandwidth_plt,
        :gap => gap_plt,
        :panel => panel
    )

end

###########################################
### Plotting for a Series of Data Integral
###########################################

function scaled_overlay(f::Function; chnl::AbstractChannel, outdir=:none, kargs...)
    println("Plotting Scaled Overlay")
    plt = plot(; frame = :box, legend=:outertopright, kargs...)
    for exdict ∈ chnl.data

        di = exdict[:data]
        data = f.((di.data)[1][:])./exdict[:ucvol]
        dom  = getindex.(getfield(di.chart,1).base,1)

        plot!(plt, dom, data,
            #label="$(round(exdict[:angle],digits=3))"
            label=""
            )
        #display(plt)
    end
    (outdir != :none) && Plots.pdf(plt, "$(mkpath(outdir))/scaled-overlay-$f.pdf")
    (outdir != :none) && Plots.png(plt, "$(mkpath(outdir))/scaled-overlay-$f.pdf")

    plt
end

function waterfall(f::Function; chnl::AbstractChannel, outdir=:none, kargs...)
    println("Plotting Waterfall")
    datamax = 0.0
    for (i,exdict) ∈ enumerate(chnl.data)
        datamax = max(datamax, (f.(exdict[:data].data[1][:])./exdict[:ucvol])...)
    end

    plt = plot(; frame = :box, legend=:outertopright, kargs...)
    for (i,exdict) ∈ enumerate(reverse(chnl.data))

        di = exdict[:data]
        data = f.(di.data[1][:])./exdict[:ucvol]
        dom  = getindex.(getfield(di.chart,1).base,1)

        plot!(plt, dom, data .+ exdict[:angle]/60.0*10datamax,
            label = "",
            xlims = (dom[1]-0.075(dom[end]-dom[1]),dom[end]),
            annotate=(-0.0375(dom[end]-dom[1]),exdict[:angle]/60.0*10datamax,Plots.text("$(round(exdict[:angle],digits=3))ᵒ",6))
            )
        #display(plt)
    end
    (outdir != :none) && Plots.pdf(plt, "$(mkpath(outdir))/waterfall-$f")
    (outdir != :none) && Plots.png(plt, "$(mkpath(outdir))/waterfall-$f")

    plt
end

##################################################
#### Scaling Plots
##################################################

function dm_scaling(datafile::String, plotdir)
    dm_scaling(; BSON.load(datafile)..., plotdir=plotdir)
end

function dm_scaling(dict::AbstractDict)
    dm_scaling(; dict...)
end

function dm_scaling(; dims, avg, std, asd, datatype, plotdir, handle, avgx, fitd, dimx, kargs...)

    ttl  = "dim scaling: $asd $datatype"

    plt = plot(
        legend = :topleft,
        margin = 5Plots.mm,
        xguide = "Hamiltonian Dimension",
        yguide = "node days ()",
        frame  = :box,
        title  = ttl,
        right_margin=5Plots.mm,
    )
    # ns/pt -> s/pt -> hr/pt -> hr -> days -> node days / core
    units = 1.0 / 1e9 / 3600 * 5e4 / 24 / 128

    scatter!(plt, dims, avg.*units, m=3, yerror=std.*units, label = "")
    plot!(plt, dimx, avgx.*units, label = "" )


    annotate!((dims[1]+1/3*(dims[end]-dims[1]),avgx[end-10]*units,Plots.text("T = τ dᵅ \n τ: $(round(fitd[:p][1]*units,sigdigits=4)) \n α: $(round(fitd[:p][2],sigdigits=4)) \n cost: $(round(sum(avg.*units),digits=3)) ")))

    Plots.pdf(plt, "$plotdir/$handle.pdf")

    plt

end

function core_scaling(dict)
    core_scaling(; dict...)
end

function core_scaling(dict,plotdir)
    core_scaling(;dict...,plotdir=plotdir)
end

function core_scaling(; avgs, wrk_samples, asd, datatype, comargs, cachedir, plotdir, pool = default_worker_pool(), kargs...)

    plt = plot(
            title = "Core Scaling",
            xguide = "1/Ncores",
            yguide = "Tₙ/T₀",
            margin = 10Plots.mm
        );
    for (k,mn) ∈ enumerate(comargs)
        avg = avgs[:,k]

        hd  = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")
        plt0 = plot( wrk_samples, avg,
                title = "Worker Scaling - $(size(hd.h_ops.h,1))",
                xguide = "Cores",
                yguide = "Time(s)",
                label  = "",
                margin = 10Plots.mm,
                color  = :red,
                xlims  = (wrk_samples[1],wrk_samples[end]),
                frame  = :box,
            );
        Plots.pdf(plt0, "$plotdir/$asd-$datatype-scaling-test-$(mn)")

        plt1 = plot( 1 ./ wrk_samples, avg,
                title = "Core Scaling - $(size(hd.h_ops.h,1)) ",
                xguide = "1/Ncores",
                yguide = "Time(s)",
                label  = "",
                margin = 10Plots.mm,
                color  = :red,
                frame  = :box
            );
        plot!(plt1, 1 ./ wrk_samples, max(avg...) ./ wrk_samples, label  = "",color=:black);
        Plots.pdf(plt1, "$plotdir/core-scaling-$asd-$datatype-scaling-test-$(mn)-inv")

        plot!(plt,  1 ./ wrk_samples, avg ./ max(avg...), label  = "$(size(hd.h_ops.h,1))", legend = :topleft);
    end

    plot!(plt, 1 ./ wrk_samples, 1 ./ wrk_samples, label  = "ref", color = :black, legend = :topleft);
    Plots.pdf(plt,"$plotdir/core-scaling-$asd-$datatype-scaling-test-$(nworkers())")

    plt
end

function mem_scaling(; com_dict, handle, plotdir, comargs, args...)
    plt = scatter(getindex.(com_dict|>values|>collect,:dim),getindex.(com_dict|>values|>collect,:top);
        label="",
        margin=8Plots.mm
        )
    # SolidState.make_models(asd, comargs..., cachedir=cachedir)
    Plots.pdf(plt,"$plotdir/$handle")

    plt
end

function integral_convergence()

end

function resource_estimation()

end


end
