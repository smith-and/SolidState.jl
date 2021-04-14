
#Path Plotting
function path_plot(title, odimeter, band_data, tick_odimeters, tick_symbols ,bz_c, corners; ylim=:all, size=(600,300))::(T where T <: AbstractPlot)
    #Start by composing the path plot
    plt= plot(yrotation=60, xtickfonthalign=:center, xtickfontsize=6, yguide ="E(eV)", yguidefontsize=6)# left_margin=1cm, right_margin=2cm, bottom_margin=1cm,top_margin=1cm)#, inset_subplots=(1, bbox(0, 0, 0.1, 0.1, :top,:left)))
    plot!(title=title, title_pos=:right, titlefontvalign=:center)
    plot!(    minorgrid=false, framestyle=:box, grid=:x, gridalpha=.7 , gridlinewidth=1, gridstyle=:dashdot)
    xticks!(  tick_odimeters, tick_symbols)
    xlims!(   odimeter[1],odimeter[end])
    ylims!(min(band_data...),max(band_data...))
    #ylim!=:all && ylims!(ylim...)
    plot!(    odimeter, transpose(band_data), label=false)
    #=#println("odimeter\n", odimeter,"band_data\n", transpose(band_data))
    #Make an inset plot and then add series data
    plot!(plt, size=size, inset_subplots = [(1, bbox(0.01,0.01,0.3w,0.3w,:bottom,:right))],  bg_color_inside=RGBA(1,1,1,0), framestyle=:none, grid=:no, axes=false, subplot=2)
    plot!(plt[2], collect(eachslice(hcat(bz_c...);dims=1))..., color=:black, label=false)
    plot!(plt[2], collect(eachslice(hcat(corners...);dims=1))..., line=(:red), label=false)
    =#
    plt
end

#Information for the path plot
#titletag=" "
#ex["val"]==:com && (titletag = ("θ = "*string(round(Lattice.cθ(Dict(ex["kwargs"])[:mn]...)*180/π; digits=2))*titletag))
#ex["v"]==:hss && (titletag = ("\\theta="*string(round(Lattice.cθ(ex["kwargs"].mn...)*180/π; digits=2))*titletag))
function many_path_plot(ex::Dict{String,Any}, sd::SpectralData{Int64,Float64}, s_idx::Int64; ylim=:all)::(Vector{T} where T<: AbstractPlot)
    #Constructing the tick marks for each plot
    n_path                            = ex["sampling:path"]["inputs"][s_idx].nbz;
    hs_bz_leger                       = ex["sampling:path"]["base"].info(ex, n_path)[3]::Dict{String,Vector{Float64}};
    if  ex["input"]["asd"].s == 0
        title                             = string(round(180/π*Lattice.cθ(1,1);digits=2))*" "*string(ex["input"]["asd"].asd)
    else
        title                             = string(round(180/π*Lattice.cθ(ex["input"]["asd"].m, ex["input"]["asd"].m+ex["input"]["asd"].s);digits=2))*" "*string(ex["input"]["asd"].asd)
    end
    (paths, corner_indices, odimeter) = Analysis.path_points(hs_bz_leger, ex["input"]["paths"], n_path)
    tick_odimeters                    = [getindex.(Ref(odimeter[i]),corner_indices[i])        for i=1:length(corner_indices)]
    tick_symbol                       = [map(x-> (x==="Γ" ? "Γ" : x), ex["input"]["paths"][i])   for i=1:length(corner_indices)]

    #Information for the inset plot
    asdg=ASDGeometry(ex["asd"])
    bz_c = getindex.(asdg["bz_c"],Ref(1:2))
    bz_c = push!([bz_c...],bz_c[1])
    leg = asdg["bz_hs"]
    corners = [getindex.(get.(Ref(leg), path,"oof"),Ref(1:2)) for path ∈ ex["input"]["paths"]]

    #Call the Path Plot routine for each path in the list
    median_i = floor(size(sd.bands.data,1)/2)
    #ylimit=[[min(sd.bands.data[max(Int(median_i-ylim[1]),1),:,i]...),max(sd.bands.data[min(Int(median_i+ylim[2]),size(sd.bands.data,1)),:,i]...)]   for i=1:length(corner_indices)]
    #ylimit=[[min(sd.bands.data[max(Int(median_i-ylim[1]),1),:,i]...),max(sd.bands.data[min(Int(median_i+ylim[2]),size(sd.bands.data,1)),:,i]...)]   for i=1:length(corner_indices)]
    #ylim=ylimit[i]
    p_plots = [path_plot(title, odimeter[i], sd.bands.data[:,:,i],tick_odimeters[i], tick_symbol[i], bz_c, corners[i] )   for i=1:length(corner_indices)]::(Vector{T} where T<: AbstractPlot)

    #display(tick_odimeters[1])
    #display(tick_symbol[1])

    p_plots
end

function velo_path_plot(ex::Dict{String,Any}, sd::SpectralData{Int64,Float64}, s_idx::Int64; ylim=:all)
    #Information for the path plot
    n_path                            = ex["sampling:path"]["inputs"][s_idx].nbz;
    hs_bz_leger                       = ex["sampling:path"]["base"].info(ex, n_path)[3]::Dict{String,Vector{Float64}};
    title                             = ("θ = "*string(round(Lattice.cθ(ex["mn"]...)*180/π; digits=2))*" "*ex["asd-base"])
    (paths, corner_indices, odimeter) = Analysis.path_points(hs_bz_leger, ex["paths"], n_path)
    tick_odimeters                    = [getindex.(Ref(odimeter[i]),corner_indices[i])        for i=1:length(corner_indices)]
    tick_symbol                       = [map(x-> (x==="Γ" ? "Γ" : x), ex["paths"][i])   for i=1:length(corner_indices)]
    #Information for the inset plot
    asdg=ASDGeometry(ex["asd"])
    bz_c = getindex.(asdg["bz_c"],Ref(1:2))
    bz_c = push!([bz_c...],bz_c[1])
    corners = [getindex.(get.(Ref(asdg["bz_hs"]), path,"oof"),Ref(1:2)) for path ∈ ex["paths"]]
    #Call the Path Plot routine for each path in the list
    ax = 2
    #median_i = (Int(floor(size(sd.band_grads.data,1)/2)))::Int64
    #ylimit=[[min(abs2.(sd.band_grads.data[:,i,max(Int(median_i-ylim[1]),1),ax])...),max(abs2.(sd.band_grads.data[:,i,min(Int(median_i+ylim[2]),size(sd.band_grads.data,3)),ax])...)]   for i=1:length(corner_indices)]
    p_plots = [path_plot(title, odimeter[i], sqrt.(abs2.(sd.band_grads.data[:,i,:,1]).+abs2.(sd.band_grads.data[:,i,:,2])),tick_odimeters[i], tick_symbol[i], bz_c, corners[i] ; ylim=:all)   for i=1:length(corner_indices)]
end

function velo_path_angle_plot(ex::Dict{String,Any}, sd::SpectralData{Int64,Float64}, s_idx::Int64; ylim=:all)
    #Information for the path plot
    n_path                            = ex["sampling:path"]["inputs"][s_idx].nbz;
    hs_bz_leger                       = ex["sampling:path"]["base"].info(ex, n_path)[3]::Dict{String,Vector{Float64}};
    title                             = ("θ = "*string(round(Lattice.cθ(ex["mn"]...)*180/π; digits=2))*" "*ex["asd-base"])
    (paths, corner_indices, odimeter) = Analysis.path_points(hs_bz_leger, ex["paths"], n_path)
    tick_odimeters                    = [getindex.(Ref(odimeter[i]),corner_indices[i])        for i=1:length(corner_indices)]
    tick_symbol                       = [map(x-> (x==="Γ" ? "Γ" : x), ex["paths"][i])   for i=1:length(corner_indices)]
    #Information for the inset plot
    asdg=ASDGeometry(ex["asd"])
    bz_c = getindex.(asdg["bz_c"],Ref(1:2))
    bz_c = push!([bz_c...],bz_c[1])
    corners = [getindex.(get.(Ref(asdg["bz_hs"]), path,"oof"),Ref(1:2)) for path ∈ ex["paths"]]
    #Call the Path Plot routine for each path in the list
    ax = 2
    #median_i = (Int(floor(size(sd.band_grads.data,1)/2)))::Int64
    #ylimit=[[min(abs2.(sd.band_grads.data[:,i,max(Int(median_i-ylim[1]),1),ax])...),max(abs2.(sd.band_grads.data[:,i,min(Int(median_i+ylim[2]),size(sd.band_grads.data,3)),ax])...)]   for i=1:length(corner_indices)]
    p_plots = [path_plot(title, odimeter[i], (180/π).*atan.(abs.(sd.band_grads.data[:,i,:,2])./abs.(sd.band_grads.data[:,i,:,1])),tick_odimeters[i], tick_symbol[i], bz_c, corners[i] ; ylim=:all)   for i=1:length(corner_indices)]
end



#Path Sampling Primitive
function section_function(type::Val{:path},     programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; kwargs...)
    ylim=[10,10];
    sd = Program.do_sampling(:path, programs, ex_idx, s_idx)
    ex = programs["executables"][ex_idx]
    p_plots = Analysis.many_path_plot(ex,sd, s_idx;ylim=ylim)
    plot(p_plots..., layout=(1,length(p_plots)))
end

function section_function(type::Val{:pathsep},     programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; kwargs...)
    ylim=[10,10];
    sd = Program.do_sampling(:path, programs, ex_idx, s_idx)
    ex = programs["executables"][ex_idx]
     Analysis.many_path_plot(ex,sd, s_idx;ylim=ylim)
end

function section_function(type::Val{:velopath},     programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; kwargs...)
    ylim=[10,10];
    sd = Program.do_sampling(:path, programs, ex_idx, s_idx)
    ex = programs["executables"][ex_idx]
    p_plots =  Analysis.velo_path_plot(ex,sd, s_idx;ylim=:all)
    plot(p_plots..., layout=(1,length(p_plots)),size=(1200,600))
end

function section_function(type::Val{:velopathangle},     programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; kwargs...)
    ylim=[10,10];
    sd = Program.do_sampling(:path, programs, ex_idx, s_idx)
    ex = programs["executables"][ex_idx]
    p_plots =  Analysis.velo_path_angle_plot(ex,sd, s_idx;ylim=:all)
    plot(p_plots..., layout=(1,length(p_plots)),size=(1200,600))
end
