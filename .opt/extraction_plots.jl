#Overlay
function overlay(f::Function; chnl::AbstractChannel, outdir=:none, kargs...)
    println("Plotting Overlay")
    plt = plot(;frame = :box, legend=:outertopright, kargs...)
    for exdict ∈ chnl.data

        di = exdict[:data]
        data = f.(di.data[1][:])
        dom  = getindex.(getfield(di.chart,1).base,1)

        plot!(plt, dom, data,
            #label="$(round(exdict[:angle],digits=3))"
            label=""
        )
        #display(plt)
    end

    (outdir != :none) && Plots.pdf(plt, "$(mkpath(outdir))/overlay-$f.pdf")
    (outdir != :none) && Plots.png(plt, "$(mkpath(outdir))/overlay-$f.pdf")

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

#GIF
function scaled_overlay_gif(f::Function; chnl::AbstractChannel, outdir=:none, kargs...)
    println("Making Overlay GIF")
    anim = Animation()
    for exdict ∈ chnl.data

        di = exdict[:data]
        data = f.(di.data[1][:])./exdict[:ucvol]
        dom  = getindex.(getfield(di.chart,1).base,1)

        plt = plot(dom, data,
            legend=:outertopright,
            #label="$(round(exdict[:angle],digits=3))",
            label="",
            xlim = (dom[1],dom[end]),
            ylim = (0.0,0.003),
            frame = :box,
            kargs...
        )
        frame(anim,plt)
    end
    (outdir != :none) && gif(anim,"$(mkpath(outdir))/scaled-overlay-$f.gif",fps=1.5)

    nothing
end

#Waterfall - Offset overlay


function isowaterfall(f::Function; chnl::AbstractChannel, outdir=:none, kargs...)
    println("Plotting Equally Spaced Waterfall")
    datamax = 0.0
    for (i,exdict) ∈ enumerate(chnl.data)
        datamax = max(datamax, (f.(exdict[:data].data[1][:])./exdict[:ucvol])...)
    end

    plt = plot(; frame = :box, legend=:outertopright, kargs...)
    for (i,exdict) ∈ enumerate(reverse(chnl.data))

        di = exdict[:data]
        data = f.(di.data[1][:])./exdict[:ucvol]
        dom  = getindex.(getfield(di.chart,1).base,1)

        plot!(plt, dom, data .+ (length(chnl.data)-i)*datamax/10,
            label="",
            xlims = (dom[1]-0.075(dom[end]-dom[1]),dom[end]),
            annotate=(-0.0375(dom[end]-dom[1]),(length(chnl.data)-i)*datamax/10,Plots.text("$(round(exdict[:angle],digits=3))ᵒ",6))
        )
        #display(plt)
    end
    (outdir != :none) && Plots.pdf(plt, "$(mkpath(outdir))/iso-waterfall-$f")
    (outdir != :none) && Plots.png(plt, "$(mkpath(outdir))/iso-waterfall-$f")

    plt
end

#Heatmap
function cutheatmapf(f::Function; chnl::AbstractChannel, outdir=:none, xlims, kargs...)
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
        margin = 10Plots.mm,
        yguide = "twist angle (degree)",
        xguide = "driving frequency (eV)",
        kargs...
    )
    #display(plt)

    (outdir != :none) && Plots.pdf(plt, "$(mkpath(outdir))/cutheatmap-$f")
    (outdir != :none) && Plots.png(plt, "$(mkpath(outdir))/cutheatmap-$f")

    plt
end


#Heatmap
function heatmapf(f::Function; chnl::AbstractChannel, outdir=:none, kargs...)
    println("Plotting Heatmap")
    angles = getindex.(chnl.data,:angle)
    base = getindex.(getfield(chnl.data[1][:data].chart,1).base,1)
    hm_data = Matrix{Float64}(undef,length.((base,angles)))

    for (i,dict) ∈ enumerate(chnl.data)
        hm_data[:,i] .= f.(dict[:data].data[1][:])./dict[:ucvol]
    end

    plt = Plots.heatmap(base,angles,hm_data';
        color  = cgrad(:bilbao),
        frame  = :box,
        margin = 10Plots.mm,
        yguide = "twist angle (degree)",
        xguide = "driving frequency (eV)",
        kargs...
    )
    #display(plt)

    (outdir != :none) && Plots.pdf(plt, "$(mkpath(outdir))/heatmap-$f")
    (outdir != :none) && Plots.png(plt, "$(mkpath(outdir))/heatmap-$f")

    plt
end

function angle_dep(f::Function, (bmin,bmax); chnl::AbstractChannel, outdir=:none, kargs...)
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
        margin=10Plots.mm,
        xlim=(angles[2],angles[end]),
        title = "Angle Dependence",
        xaxis = "twist angle (deg.)",
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

    (outdir != :none) && Plots.pdf(plt_all, "$(mkpath(outdir))/angle-dep-$f")
    (outdir != :none) && Plots.png(plt_all, "$(mkpath(outdir))/angle-dep-$f")

    plt_all
end


function angle_dep_gif(f::Function, (bmin,bmax); chnl::AbstractChannel, outdir=:none, kargs...)
println("Making Angle Dependence GIF")
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
    plt_all = plot(; frame = :box, margin=10Plots.mm, xlim=(angles[1],angles[end]), kargs...)

    for (i,b) ∈ enumerate(base)
        if bmin < b < bmax
            Plots.plot!(plt_all, angles,hm_data[i,:],
                label = "",
            )
            #display(plt_all)

            plt = Plots.plot(angles,hm_data[i,:],
                label = "$b",
                margin=10Plots.mm,
                xlim  = (angles[2],angles[end]),
                ylim  = (0.0,d_max)
            )
            frame(anim, plt)
        end
    end

    (outdir != :none) && gif(anim, "$outdir/angle-dep-$f.gif", fps = 10)
    nothing
end
