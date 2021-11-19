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
