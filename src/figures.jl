module Figure

using Revise
using Distributed, Dates, OrderedCollections, BSON, Plots
using LinearAlgebra, SharedArrays, StaticArrays
using Mux, WebIO, Interact, InteractiveUtils
using CubicSplines, Roots, SpecialFunctions, HCubature
using Base: Threads
using ..SolidState
using SolidState.Plot

########################################
#### Plot Panels
########################################

function shg_band_panel(seriesdata, sectiondata; f, xlims, rng, size, idx)
    plt2 = Dict(
        :bandbook => Plot.band_series(seriesdata[:bands]; idxs=2:length(seriesdata[:bands].data)),
        :heatm    => Plot.cutheatmapf(seriesdata[:shg]; f=f, xlims=xlims),
        :andep    => Plot.angle_dep(seriesdata[:shg]; f=f, rng=rng),
    )
    bands = plt2[:bandbook]

    sectionplot = Plot.section.(sectiondata,f,Ref(idx))
    angsectionplot = Plot.section.(sectiondata,angle,Ref(idx))

    plot(
        bands[:conduct_plt], bands[:valence_plt], bands[:gap], bands[:bandwidth], plt2[:heatm], plt2[:andep],sectionplot..., angsectionplot...,
        layout = (@layout [ [ b11 ; b12 ; b21 b22 ] [ c1 ; c2 ] grid(length(sectionplot),2) ]),
        size = size
        )

end

# Dict(
#     :cond_heatmap => cond_heatmap,
#     :conduct_plt => conduct_plt,
#     :vale_heatmap => vale_heatmap,
#     :valence_plt => valence_plt,
#     :bandwidth => bandwidth_plt,
#     :gap => gap_plt,
#     :panel => panel
# )


end
