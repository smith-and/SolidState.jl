#Plotting and Exporting
include("../plotting/equilibrium.jl")

#Equilibrium Expectations
function quad_plot(v::Val{:GCE}, name, section, ex, q_idx)
    δ_idx=1;
    T_idx=1;

    #optical responses
    filling_plot     = plot(rand(10))
    energy_plot      = plot(rand(10))
    free_energy_plot = plot(rand(10))

    plot_names = [:filling, :energy, :free_energy]
    plots      = [filling_plot,energy_plot,free_energy_plot]

    plot_paths = string.(Ref(ex["export/path"]*"/plots/"*string(name)*"/"), string.(plot_names), Ref("-s"*string(q_idx)*".pdf"))
    Plots.pdf.(plots, plot_paths)

    Dict{Symbol,AbstractPlot}(Pair.(plot_names,plots))
end

function quad_export(v::Val{:GCE}, name, section, ex, q_idx)
    data_path=ex["export/path"]* "/data/GCE"
    section.ParticleNumber( data_path* "/particle_number-s"*string(q_idx))
    section.Energy(data_path*"/energy-s"*string(q_idx))
    section.FreeEnergy(data_path*"/free_energy-s"*string(q_idx))

    return data_path
end
