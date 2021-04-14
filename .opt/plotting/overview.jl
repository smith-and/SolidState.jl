function model_overview_plot(ex_id, program)::Nothing
    gr()

    ex          = program["executables"][ex_id]
    exportpath  = program["exportpath"]

    #Density of States Plot
    dos_plot= ex["ex-Analysis.DosData{Int64,Float64}"]["plot"][end][:dos]|>Plots.hdf5plot_read
    ylims!(dos_plot, (ex["bandmin"],ex["bandmax"]))

    #Plot of Lattice Unit Cell
    uc_plots = Analysis.uc_plot_panels(ex)
    lattice_plot = plot(
        uc_plots[:XY],
        uc_plots[:XZ],
        uc_plots[:YZ],
        uc_plots[:labels],
        layout=grid(2,2), size=(1000,1000)
    )

    #Composing the Overview Alignment
    l1 = @layout [ e g f g ]
    spectra_overview = plot(
        ex["sampling:path"]["plot"][end][:path1]|>Plots.hdf5plot_read,
        ex["sampling:path"]["plot"][end][:path2]|>Plots.hdf5plot_read,
        ex["sampling:path"]["plot"][end][:path3]|>Plots.hdf5plot_read,
        dos_plot,
        layout=l1, size=(1500,1000)
    )

    l2=@layout grid(2,4)
    optical_overview = plot(
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:optcondABS]|>Plots.hdf5plot_read,
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:optcondRE]|>Plots.hdf5plot_read,
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:optcondIM]|>Plots.hdf5plot_read,
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:optcondARG]|>Plots.hdf5plot_read,
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:shgABS]|>Plots.hdf5plot_read,
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:shgRE]|>Plots.hdf5plot_read,
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:shgIM]|>Plots.hdf5plot_read,
        ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:shgARG]|>Plots.hdf5plot_read,
        layout=l2, size=(1500,1000)
    )

    #ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:cpgeABS],
    #ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:cpgeRE],
    #ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:cpgeIM],
    #ex["ex-Analysis.OpticalData{Int64,Float64}"]["plot"][end][:cpgeARG],
    #display(optical_overview)

    #This is the export step
    mkpath(exportpath*"/lattice")
    mkpath(exportpath*"/spectra")
    mkpath(exportpath*"/optical")

    Plots.pdf(lattice_plot, exportpath*"/lattice/$(ex["model_id"])")
    Plots.pdf(spectra_overview, exportpath*"/spectra/$(ex["model_id"])")
    Plots.pdf(optical_overview, exportpath*"/optical/$(ex["model_id"])")
    display.([spectra_overview,optical_overview])

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
