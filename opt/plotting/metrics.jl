export matrixplot
function matrixplot(f::Function,A::AbstractArray)
    heatmap(collect.(UnitRange.(1,size(A)))...,f.(A),aspect_ratio=1, color=:viridis, yaxis=(:flip), xmirror=true)
    ylims!(1,size(A,2))
    xlims!(1,size(A,1))
end

#Metric Exporting
#Named Tuple to BSON conversion
function write_metrics(path,inputs, name="input")
    inputs==:none ? nothing : bson(  path*"/"*name*".bson"; keys = keys(inputs[1]), values=values.(inputs))
end

#Program Metric Analysis and Plotting
function metric_plots(section_info)
    (scatter((getindex.(section_info["metric"],:t)[:]./3600),     title="timing (h)"),
     scatter((getindex.(section_info["metric"],:bytes)[:])./1e9, title="bytes (GB)"),
     scatter((getindex.(section_info["metric"],:bytes)[:]).*(getindex.(section_info["metric"],:t)[:])./1e9, title="alloc resource (GBs)"),
     scatter(getindex.(section_info["metric"],:gc)[:],    title="garbage collection (bytes)")
    ) |> plots-> Dict{Symbol,AbstractPlot}(
    :time     => plots[1],
    :bytes    => plots[2],
    :resource => plots[3],
    :gc       => plots[4],
    )
end

function metric_export(ex, metric_path, name, i=0)
    write_metrics(metric_path,ex["metric"],"name-"*string(i))

    ex["metric:plot"] = metric_plots(ex)
    for key ∈ keys(ex["metric:plot"])
        Plots.pdf(ex["metric:plot"][key],metric_path*"/"*name*"-"*string(key)*"-"*string(i)*".pdf")
    end
end

function basic_metric_analysis(program)
    metric_export(program["build:metrics"], program["exportpath"]*"/metrics/build", "asd", 0)
    metric_export(program["model:metrics"], program["exportpath"]*"/metrics/model", "model", 0)

    for i=eachindex(program["executables"])
        metric_export(program["executables"][i]["sampling:path"], program["exportpath"]*"/metrics/sampling" , "path", i)
        metric_export(program["executables"][i]["quad:Analysis.SpectralIntegral"], program["exportpath"]*"/metrics/quad" ,     "spectral",  i)
        metric_export(program["executables"][i]["quad:Analysis.OpticalResponse"],  program["exportpath"]*"/metrics/quad" ,     "optical",  i)
    end

end
