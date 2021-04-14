function model_overview_plot(ex_id, program)::Nothing
    gr()

    ex          = program["executables"][ex_id]
    exportpath  = program["exportpath"]
    #High Symmetry Band Structure
    path1 = ex["sampling:path"]["plot"][end][:path1]
    path_min = min((x->min(x...)).(getindex.(getfield.(Plots.hdf5plot_read(path1).series_list[1:(end-2)],:plotattributes),:y))...)
    path_max = max((x->max(x...)).(getindex.(getfield.(Plots.hdf5plot_read(path1).series_list[1:(end-2)],:plotattributes),:y))...)
    #ylims!(path1, (path_min,path_max))

    #Density of States Plot
    dos_plot= ex["quad:Analysis.SpectralIntegral"]["plot"][end][:dos]|>Plots.hdf5plot_read
    dos_h_plot= ex["quad:Analysis.SpectralIntegral"]["plot"][end][:dos_h]|>Plots.hdf5plot_read
    ylims!(dos_plot, (path_min,path_max))
    jdos_plot= ex["quad:Analysis.SpectralIntegral"]["plot"][end][:jdos]|>Plots.hdf5plot_read
    shdos_plot= ex["quad:Analysis.SpectralIntegral"]["plot"][end][:shdos]|>Plots.hdf5plot_read

    #Plot of Lattice Unit Cell
    uc_plots = Analysis.uc_plot_panels(ex)
    l3 = @layout [x ; y ; z ]
    moredos_overview = plot(
        dos_h_plot,
        jdos_plot,
        shdos_plot,
        layout=l3, size=(500,1000)
    )
    #display(moredos_overview)

    #Composing the Overview Alignment
    l1 = @layout [ grid(4,1) e g f g ]
    spectra_overview = plot(
        uc_plots[:XY],
        uc_plots[:XZ],
        uc_plots[:YZ],
        uc_plots[:labels],
        ex["sampling:path"]["plot"][end][:path1]|>Plots.hdf5plot_read,
        ex["sampling:path"]["plot"][end][:path2]|>Plots.hdf5plot_read,
        ex["sampling:path"]["plot"][end][:path3]|>Plots.hdf5plot_read,
        dos_plot,
        layout=l1, size=(1500,1000)
    )
    #display(spectra_overview)

    l2=@layout grid(2,4)
    optical_overview = plot(
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:optcondABS]|>Plots.hdf5plot_read,
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:optcondRE]|>Plots.hdf5plot_read,
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:optcondIM]|>Plots.hdf5plot_read,
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:optcondARG]|>Plots.hdf5plot_read,
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:shgABS]|>Plots.hdf5plot_read,
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:shgRE]|>Plots.hdf5plot_read,
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:shgIM]|>Plots.hdf5plot_read,
        ex["quad:Analysis.OpticalResponse"]["plot"][end][:shgARG]|>Plots.hdf5plot_read,
        #ex["quad:Analysis.OpticalResponse"]["plot"][end][:cpgeABS],
        #ex["quad:Analysis.OpticalResponse"]["plot"][end][:cpgeRE],
        #ex["quad:Analysis.OpticalResponse"]["plot"][end][:cpgeIM],
        #ex["quad:Analysis.OpticalResponse"]["plot"][end][:cpgeARG],
        layout=l2, size=(1500,1000)
    )
    #display(optical_overview)

    #This is the export step
    mkpath(exportpath*"/spectra")
    mkpath(exportpath*"/optical")
    mkpath(exportpath*"/moredos")

    Plots.pdf(moredos_overview, exportpath*"/moredos/$(ex["model_id"])")
    Plots.pdf(spectra_overview, exportpath*"/spectra/$(ex["model_id"])")
    Plots.pdf(optical_overview, exportpath*"/optical/$(ex["model_id"])")
    display.([spectra_overview,optical_overview,moredos_overview])


    nothing
end

#Path and Quadrature Sampling
function section_function(type::Val{:overview}, program::Dict{String,Any}, ex_idx::Int64=1, s_idx::Int64=1; kwargs...)::Nothing
    Analysis.do_sampling(:path, program, ex_idx, s_idx)
    Analysis.do_quad(Val(:BZ), :spectral, program, ex_idx, s_idx)
    Analysis.do_quad(Val(:BZ), :optical,  program, ex_idx, s_idx)

    Analysis.model_overview_plot(ex_idx, program)
end

function post_process_models(program)::Nothing
    for ex_id ∈ eachindex(program["executables"])
        program["executables"][ex_id]["plot?"] && model_overview_plot(ex_id, program);
        nothing
    end
    nothing
end
