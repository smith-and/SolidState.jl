#Plotting
include("../plotting/path.jl")
include("../plotting/mesh.jl")

#Exporting

#Selected Plot alignments
#Path Plotting
function sampling_plot(v::Val{:path}, type, sd, ex, s_idx)::Nothing
    s_tag = "sampling:"*string(type)
    plots = many_path_plot(ex,sd, s_idx)
    plot_path=ex["export/path"]
    plot_leg=Dict{Symbol,String}()
    a_plot = plot()
    mkpath(plot_path*"/hdf5")
    for i=1:length(plots)
        a_plot = plots[i]
        #hdf5_path = string(plot_path*"/hdf5/path-", i, "-s-"*string(s_idx)*".hdf5")
        pdf_path = string(plot_path*"/path-", i, "-s-"*string(s_idx))
        Plots.pdf(a_plot, pdf_path)
        #hdf5_plot_pdf(hdf5_path,pdf_path)
        push!(plot_leg, Symbol("path"*string(i))=>pdf_path)
    end
    ex[s_tag]["plot"][s_idx] = plot_leg
    nothing
end

#Mesh Plotting
function sampling_plot(v::Val{:mesh}, sampling_data, ex, s_idx)
    #do mesh data analysis
    #do fermi surface slices at different fillings
    #make 3D surface plot
    #make vector plots of the velocity
    #plot heatmaps of band mass
    Dict{Symbol,String}(:all=>"of it")
end

#Save Data and Plots
function sampling_export(name, sd, ex, s_idx)::Nothing
    type = name==:path ? "/paths-s$(s_idx)" : "/meshes-s$(s_idx)"
    data_path=ex["export/path"]*type
    ex["sampling:"*string(name)]["result"][s_idx]   =  data_path
    #Write the chart data for each field to a
    TensorCharts.write(sd,data_path)

#    sd.bands(     data_path*     "/bands-s"*string(s_idx))
#    sd.band_grads(data_path*"/band_grads-s"*string(s_idx))
#    sd.inv_mass(  data_path*  "/inv_mass-s"*string(s_idx))
#    sd.structure( data_path* "/structure-s"*string(s_idx))
    nothing
end
