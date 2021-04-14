
#Spectral Integral Functions
include("ldos.jl")
include("jdos.jl")
include("shdos.jl")


struct MoreDosData{Int_Type <: Int, Float_Type <: Real} <: QuadratureData{Int_Type, Float_Type}
    jdos::TensorChart{Int_Type,  1, 1, Float_Type, 2, Float_Type, 2,  Complex{Float_Type}, 5, 1}
    ldos::TensorChart{Int_Type,  1, 1, Float_Type, 2, Float_Type, 2,  Complex{Float_Type}, 5, 1}
    shdos::TensorChart{Int_Type,  1, 1, Float_Type, 2, Float_Type, 2,  Complex{Float_Type}, 5, 2}
end

function integrand(t::Type{Analysis.MoreDosData{Int64,Float64}}, cb_input::ChartSectionInput, priors::Prior, Hdim::Int64)::MoreDosData{Int64, Float64}
    MoreDosData(TensorChart.(cb_input, Ref(priors), Ref(Hdim))...)
end

#Spectral Quadrature

#=
function plot_horz_dos(chart::TensorChart; δ_idx, T_idx, μ_idx, priors, title, xlabel, ylabel)
    gr()
    #Form the plot object
    legendtitle = "δ = "*string(round(chart.base[1,δ_idx][2];digits=3))
    t_plot = plot(
        size=(600,600), framestyle=:box, xlabel=xlabel,ylabel=:none, yrotation=60,
        title=title, titlefontvalign=:bottom, title_location=:right,
        legendtitlefonthalign=:left,legendtitlefontsize=6, legendtitle=legendtitle,
        legend=:bottomright, legendfontsize=5, fglegend=RGBA(1.0,1.0,1.0,0.25), bglegend=RGBA(1.0,1.0,1.0,0.25)
    )

    #Plot the projected densities
    ωs = getindex.(chart.base[:,1],1)
    for t_idx=1:length(chart.indices)
        plot!(t_plot, ωs, imag.(chart.data[t_idx,T_idx,μ_idx,:,δ_idx]), label=chart.indices[t_idx],title_location=:right)
    end
    #Return the Plot
    t_plot
end
=#

function section_plot(name, section::MoreDosData{Int64,Float64}, ex::Dict{String,Any}, q_idx::Int64)::Nothing
    q_tag = "ex-"*string(name)

    δ_idx=1;
    T_idx=1;
    μ_idx=1;
    angle_title = string(round(180/π*Lattice.cθ(ex["input"]["asd"].m, ex["input"]["asd"].m+ex["input"]["asd"].s);digits=2))*" "*string(ex["input"]["asd"].asd)

    #Spectral Integration
    plot_names = [:jdos,:shdos,:all]
    plots0 = [
        #plot_horz_dos(section.ldos; title="LDOS",  xlabel="ω", ylabel="", δ_idx=δ_idx, T_idx=T_idx, μ_idx=μ_idx, priors = section.ldos.priors[T_idx,μ_idx]),
        plot_horz_dos(section.jdos; title="JDOS", xlabel="ω", ylabel="(1/nm)²", δ_idx=δ_idx, T_idx=T_idx, μ_idx=μ_idx, priors = section.jdos.priors[T_idx,μ_idx]),
        plot_horz_dos(section.shdos; title="SHDOS", xlabel="ω", ylabel="(1/nm)²", δ_idx=δ_idx, T_idx=T_idx, μ_idx=μ_idx, priors = section.shdos.priors[T_idx,μ_idx])
    ]

    plots = [
        plots0...,
        plot(plots0..., size=(length(plots0)*600,600),layour=grid(2,1))
    ]

    #Saving
    hdf5_paths = string.(Ref(ex["export/path"]*"/hdf5/md-"), string.(plot_names), Ref("-s"*string(q_idx)*".hdf5"))
    pdf_paths = string.(Ref(ex["export/path"]*"/md-"), string.(plot_names), Ref("-s"*string(q_idx)))

    for i = 1:length(plots)
        Plots.pdf(plots[i], pdf_paths[i])
        #hdf5_plot_pdf(hdf5_paths[i],pdf_paths[i])
    end

    #Dictionary Output
    ex[q_tag]["plot"][q_idx] = Dict{Symbol,String}(Pair.(plot_names,pdf_paths));

    nothing
end


function section_export(name, section::MoreDosData{Int64,Float64}, ex::Dict{String,Any}, q_idx::Int64; )::Nothing
    data_path=ex["export/path"]*"/spectral-s$(q_idx)"
    TensorCharts.write(section, data_path)
    ex["ex-"*string(name)]["result"][q_idx]   =  data_path
    nothing
end



function more_dos_plots(program, ex_id)
    gr()

    ex          = program["executables"][ex_id]
    exportpath  = program["exportpath"]

    dos_h_plot= ex["ex-Analysis.DosData{Int64,Float64}"]["plot"][end][:dos_h]|>Plots.hdf5plot_read
    jdos_plot= ex["ex-Analysis.DosData{Int64,Float64}"]["plot"][end][:jdos]|>Plots.hdf5plot_read
    shdos_plot= ex["ex-Analysis.DosData{Int64,Float64}"]["plot"][end][:shdos]|>Plots.hdf5plot_read
    #Dos Plots
    l3 = @layout [ x y z ]
    moredos_overview = plot(
        dos_h_plot,
        jdos_plot,
        shdos_plot,
        layout=l3, size=(1000,500)
    )
    mkpath(exportpath*"/moredos")
    Plots.pdf(moredos_overview, exportpath*"/moredos/$(ex["model_id"])")

end
