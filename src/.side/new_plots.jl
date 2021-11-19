
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
function bands(dict, NBs)
    plotbook = Dict(
        "band-reg-all"=>band_region_plot(       dict,(1,:end)),
        "band-val-all"=>band_region_plot(       dict,(:val,:end)),
        "band-cond-all"=>band_region_plot(      dict,(:cond,:end)),
        "band-broken-all"=>band_broken_plot(    dict,:end),
    )
    for NB ∈ NBs
        merge!(plotbook,Dict(
        "band-reg-$NB"=> band_region_plot(      dict,(-NB,-NB)),
        "band-val-$NB"=>band_region_plot(       dict,(:val,NB)),
        "band-cond-$NB"=>band_region_plot(      dict,(:cond,NB)),
        "band-broken-$NB"=>band_broken_plot(    dict,NB),
        ))
    end
    plotbook
end

function bands(NBs; dict, plotdir, kargs...)
    plotbook = bands(dict,NBs)
    SolidState.Main.book_save(plotbook,plotdir)
    plotbook
end
