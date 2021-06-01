module Plot

using Revise
using Distributed, Dates, OrderedCollections, BSON, Plots
using LinearAlgebra, SharedArrays, StaticArrays
using Mux, WebIO, Interact, InteractiveUtils
using CubicSplines, Roots, SpecialFunctions, HCubature
using Base: Threads
using ..SolidState

######################################################################
### Plotting Simple Data
######################################################################

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

#Generic Plotting Method for Charts
function integral(f=identity ; data::DataIntegral, args...)
    plts = Vector{typeof(plot())}(undef,length(data.err))
    for i∈1:length(data.err)
        plts[i] = plot(getindex.(getfield(data.dm.chart,1).base,1), f.(data.data[i][:]);
            ribbon = data.err[i],
            xlims  = (getindex.(getfield(data.dm.chart,1).base,1)[1],getindex.(getfield(data.dm.chart,1).base,1)[end]),
            legend = :topleft,
            label  = "$(data.evals[i])",
            frame  = :box,
            margins = 8Plots.mm,
            ylims = (min(0.0,(f.(data.data[i][:]))...),1.1*max(f.(data.data[i][:])...,1e-15)),
            args...
            )
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

function angle_dep(chnl::AbstractChannel; f::Function, rng, kargs...)
    (bmin,bmax) = rng
    println("Plotting Angle Dependence")
    #Prepare Data
    angles = getindex.(chnl.data,:angle)
    base = getindex.(getfield(chnl.data[1][:data].chart,1).base,1)
    hm_data = Matrix{Float64}(undef,length.((base,angles)))
    for (i,dict) ∈ enumerate(chnl.data)
        hm_data[:,i] .= f.(dict[:data].data[1][:])./dict[:ucvol]
    end

    #Plot
    d_max = max(hm_data[bmin .< base .< bmax, :]...)
    anim = Animation()
    plt_all = plot(;
        frame = :box,
        # margin=10Plots.mm,
        xlim=(angles[2],angles[end]),
        # title = "Angle Dependence",
        # xaxis = "twist angle (deg.)",
        kargs...
        )

    for (i,b) ∈ enumerate(base)
        if bmin < b < bmax
            Plots.plot!(plt_all, angles,hm_data[i,:];
                color = get(reverse(cgrad(:rainbow)),(b-bmin)/(bmax-bmin)),
                label = "",
                kargs...
            )
        end
    end

    plt_all
end

function cutheatmapf(chnl::AbstractChannel; f::Function, outdir=:none, xlims, kargs...)
    println("Plotting Heatmap")
    angles = getindex.(chnl.data,:angle)
    base = getindex.(getfield(chnl.data[1][:data].chart,1).base,1)
    hm_data = Matrix{Float64}(undef,length.((base,angles)))

    for (i,dict) ∈ enumerate(chnl.data)
        hm_data[:,i] .= f.(dict[:data].data[1][:])./dict[:ucvol]
    end


    plt = Plots.heatmap(base[xlims[1] .< base .< xlims[2]],angles,hm_data[xlims[1] .< base .< xlims[2],:]';
        color  = cgrad(:bilbao),
        frame  = :box,
        # margin = 10Plots.mm,
        # yguide = "twist angle (degree)",
        # xguide = "driving frequency (eV)",
        kargs...
    )
    #display(plt)

    (outdir != :none) && Plots.pdf(plt, "$(mkpath(outdir))/cutheatmap-$f")
    (outdir != :none) && Plots.png(plt, "$(mkpath(outdir))/cutheatmap-$f")

    plt
end

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


    annotate!((dims[1]+1/3*(dims[end]-dims[1]),avgx[end-10]*units,Plots.text("T = τ dᵅ \n τ: $(round(fitd[:p][1]*units,sigdigits=4)) \n α: $(round(fitd[:p][2],sigdigits=4)) ")))

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
    plt = scatter(getindex.(com_dict|>values|>collect,:dim),getindex.(com_dict|>values|>collect,:top),label="")
    # SolidState.make_models(asd, comargs..., cachedir=cachedir)
    Plots.pdf(plt,"$plotdir/$handle")

    plt
end

function integral_convergence()

end

function resource_estimation()

end


end
