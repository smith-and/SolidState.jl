function shg_chart_plot(chart::TensorChart; δ_idx, T_idx, μ_idx, priors)
    ωs        =  getindex.(chart.base[:,1],1)
    title     =  title
    title=" SHG"
    ylabel="χ(ω,ω)"
    xlabel="ω"

    legendtitle = "δ="*string(round(chart.base[1,δ_idx][2];digits=3))

    t_plot = plot(
                size=(300,300), framestyle=:box, xlabel=xlabel, ylabel=ylabel,yrotation=60,
                title=title, titlefontvalign=:bottom, title_location=:right,
                legendtitlefonthalign=:left,legendtitlefontsize=6, legendtitle=legendtitle,
                legend=:best, legendfontsize=5, fglegend=RGBA(1.0,1.0,1.0,0.25), bglegend=RGBA(1.0,1.0,1.0,0.25)
                )

    #ylims!(0.0,1.1*max(max.(real.(chart.data[:,T_idx,μ_idx,:,δ_idx])...),1e-5))

    for t_idx=1:length(chart.indices)
        plot!(abs_plot, ωs,  abs.(chart.data[t_idx,T_idx,μ_idx,:,δ_idx]),label=chart.indices[t_idx])
        #plot!(arg_plot, ωs,  arg.(chart.data[t_idx,T_idx,μ_idx,:,δ_idx]),label=chart.indices[t_idx])
        #plot!(re_plot,  ωs, real.(chart.data[t_idx,T_idx,μ_idx,:,δ_idx]),label=chart.indices[t_idx])
        #plot!(im_plot,  ωs, imag.(chart.data[t_idx,T_idx,μ_idx,:,δ_idx]),label=chart.indices[t_idx])
    end
    #println((chart.data[:,T_idx,:,δ_idx]))
    t_plot
end

function cpge_chart_plot(chart::TensorChart; δ_idx, T_idx, μ_idx, priors)
    ωs        =  getindex.(chart.base[:,1],1)
    title     =  title

    title=" CPGE"
    ylabel="χ(ω,-ω)"
    xlabel="ω"


    legendtitle = "δ="*string(round(chart.base[1,δ_idx][2];digits=3))
    t_plot = plot(
                size=(300,300), framestyle=:box, xlabel=xlabel, ylabel=ylabel, yrotation=60,
                title=title, titlefontvalign=:bottom, title_location=:right,
                legendtitlefonthalign=:left,legendtitlefontsize=6, legendtitle=legendtitle,
                legend=:best, legendfontsize=5, fglegend=RGBA(1.0,1.0,1.0,0.25), bglegend=RGBA(1.0,1.0,1.0,0.25)
                )
    ylims!(0.0,max(max.(abs.(chart.data[:,T_idx,μ_idx,:,δ_idx])...),1e-5))
    for t_idx=1:length(chart.indices)
        plot!(t_plot, ωs,abs.(chart.data[t_idx,T_idx,μ_idx,:,δ_idx]),label=chart.indices[t_idx],title_location=:right)
    end
    #println((chart.data[:,T_idx,:,δ_idx]))
    t_plot
end

function optσ_chart_plot(chart::TensorChart; δ_idx, T_idx, μ_idx, priors)
    ωs        =  getindex.(chart.base[:,1],1)
    title     =  title

    title=" Dyn. Cond."
    ylabel="σ(ω)"
    xlabel="ω"

    legendtitle = "δ="*string(round(chart.base[1,δ_idx][2];digits=3))
    t_plot = plot(
                size=(300,300), framestyle=:box, xlabel=xlabel, ylabel=ylabel,yrotation=60,
                title=title, titlefontvalign=:bottom, title_location=:right,
                legendtitlefonthalign=:left,legendtitlefontsize=6, legendtitle=legendtitle,
                legend=:topright, legendfontsize=5, fglegend=RGBA(1.0,1.0,1.0,0.25), bglegend=RGBA(1.0,1.0,1.0,0.25)
                )
    #ylims!(0.0,max(max.(abs.(chart.data[:,T_idx,:,δ_idx])...),1.0))
    #println(chart.data)
    for t_idx=1:length(chart.indices)
        plot!(t_plot, ωs,abs.(chart.data[t_idx,T_idx,μ_idx,:,δ_idx]),label=chart.indices[t_idx],title_location=:right)
    end
    #println((chart.data[:,T_idx,:,δ_idx]))
    t_plot
end

function plot_style(title,ylabel,xlabel,legendtitle="")
    (   size=(300,300), framestyle=:box, xlabel=xlabel, ylabel=ylabel, yrotation=60,
        title=title, titlefontvalign=:bottom, title_location=:right,
        legendtitlefonthalign=:left,legendtitlefontsize=6, legendtitle=legendtitle,
        legend=:best, legendfontsize=5, fglegend=RGBA(1.0,1.0,1.0,0.25), bglegend=RGBA(1.0,1.0,1.0,0.25)
    )
end

function chart_plot(chart::TensorChart; b_idx, p_idx, priors, title, ylabel, xlabel, factor=1)
    #Apply scaling for unit conversion
    chart.data .*= factor
    #Get x-axis range
    ωs        =  getindex.(chart.base[:,1],1)
    legendtitle = "δ="*string(round(chart.base[1,b_idx[2]][2];digits=3))

    abs_plot = plot(;plot_style("ABS"*title, ylabel,xlabel,legendtitle)...)
    ylims!( 1.1*min(min.(abs.(chart.data[:,p_idx..., b_idx...])...),1e-5),
    1.1*max(max.(abs.(chart.data[:,p_idx..., b_idx...])...),1e-5))

    arg_plot = plot(;plot_style("ARG"*title, ylabel,xlabel,legendtitle)...)
    ylims!( 1.1*min(min.(angle.(chart.data[:,p_idx..., b_idx...])...),1e-5),
    1.1*max(max.(angle.(chart.data[:,p_idx..., b_idx...])...),1e-5))

    re_plot = plot(;plot_style( "RE"*title, ylabel,xlabel,legendtitle)...)
    ylims!( 1.1*min(min.(real.(chart.data[:,p_idx..., b_idx...])...),1e-5),
    1.1*max(max.(real.(chart.data[:,p_idx..., b_idx...])...),1e-5))

    im_plot = plot(;plot_style( "IM"*title, ylabel,xlabel,legendtitle)...)
    ylims!( 1.1*min(min.(imag.(chart.data[:,p_idx..., b_idx...])...),1e-5),
            1.1*max(max.(imag.(chart.data[:,p_idx..., b_idx...])...),1e-5))

    for t_idx=1:length(chart.indices)
        plot!(abs_plot, ωs, abs.(chart.data[t_idx,p_idx...,:,b_idx...]),    label=chart.indices[t_idx])
        plot!(arg_plot, ωs, angle.(chart.data[t_idx,p_idx...,:,b_idx...]),  label=chart.indices[t_idx])
        plot!(re_plot,  ωs, real.(chart.data[t_idx,p_idx...,:,b_idx...]),   label=chart.indices[t_idx])
        plot!(im_plot,  ωs, imag.(chart.data[t_idx,p_idx...,:,b_idx...]),   label=chart.indices[t_idx])
    end

    (abs_plot,arg_plot,re_plot,im_plot)
end
