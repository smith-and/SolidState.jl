#Region Mesh Plotting
function add_frame!(plt,coordinate_mesh; rgb=(0,0,0))
    #Adding the Boundaries of the Mesh Region
    meshx       = getindex.(coordinate_mesh,1)::Array{Float64,2}
    meshy       = getindex.(coordinate_mesh,2)::Array{Float64,2}
    plot!(plt, meshx[1,:],   meshy[1,:],   label=false, color=RGB(rgb...))
    plot!(plt, meshx[end,:], meshy[end,:], label=false, color=RGB(rgb...))
    plot!(plt, meshx[:,1],   meshy[:,1],   label=false, color=RGB(rgb...))
    plot!(plt, meshx[:,end], meshy[:,end], label=false, color=RGB(rgb...))

    plt
end

#region plot on a single value across the tensor indices
function tensor_levelset_plot(coordinate_mesh, data, coordinate_map, band_colors, lvlset)
    #Creating an empty plot
    plt = plot(frame=:none, aspectratio=2/sqrt(3.0))
    #Adding the Boundaries of the Mesh Region
    add_frame!(plt, coordinate_mesh)
    #Doing the marching squares calculation for each band and plotting
    for i=1:length(band_data)
        contour_levels = Contour.levels(Contour.contours(collect(0:(size(coordinate_mesh,1)-1)),collect(0:(size(coordinate_mesh,2)-1)), data[i,:,:], [lvlset]))
        for cl ∈ contour_levels
            for line in lines(cl)
                xs, ys = eachslice(hcat(map(coordinate_map, eachslice(hcat(coordinates(line)...);dims=1))...);dims=1) # coordinates of this line segment
                plot!(plt, xs, ys, color=band_colors[i],label=false) # pseuod-code; use whatever plotting package you prefer
            end
        end
    end

    plt

end

#Region Plot on a single tensor index across values of the domain
rgb_color_function(mn, mx, (r,g,b)) = lvl-> (abs((mx - lvl)/(mx - mn)) |> x-> (RGB(x*r,x*g,x*b)))

function tensor_topograph(coordinate_mesh, data, coordinate_map, band_color_arg)
    #Creating an empty plot
    plt = plot(frame=:none, aspectratio=2/sqrt(3.0))
    #Adding the Boundaries of the Mesh Region
    add_frame!(plt, coordinate_mesh)
    #Adding the contour lines weights with respect to the extrema provided
    c_f = rgb_color_function(band_color_arg...)
    contour_levels = Contour.levels(Contour.contours(collect(0:(size(coordinate_mesh,1)-1)),collect(0:(size(coordinate_mesh,2)-1)), data))
    for cl ∈ contour_levels
        lvl = level(cl) # the z-value of this contour level
        for line in lines(cl)
            xs, ys = eachslice(hcat(map(coordinate_map, eachslice(hcat(coordinates(line)...);dims=1))...);dims=1) # coordinates of this line segment
            plot!(plt, xs, ys, color=c_f(lvl),label=false) # pseuod-code; use whatever plotting package you prefer
        end
    end
    plt
end


#Mesh Sampling Primitive
function index_2_mesh(Λ::Array{Float64,2}, mesh_pt::(Tuple{N, Int64} where N))::Vector{Float64}
    Λ*[(mesh_pt.-1)...]
end

function index_2_mesh(Λ::Array{Float64,2}, mesh_pt::Vector{Float64})::Vector{Float64}
    Λ*(mesh_pt.-1)
end

function section_function(type::Val{:fermisurface}, programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; lvlset=0.0, kwargs...)
    sd = Program.do_sampling(:mesh, programs, ex_idx, s_idx)
    ex = programs["executables"][ex_idx]
    #Geometric Information
    nbz                  = ex["sampling:mesh"]["inputs"][s_idx].nbz
    Λ                    = ex["sampling:mesh"]["base"].info(ex,nbz)[1][1:2,1:2]::Array{Float64,2};
    coordinate_map(x)    = (Λ*x/nbz)
    #Extracting the coordinate grids
    integer_mesh    = (ex["sampling:mesh"]["base"].mesh)([nbz+1 for _=1:size(Λ,2)], 1, 0)[1]


    coordinate_mesh = index_2_mesh.(Ref(Λ/nbz), integer_mesh)
    plt = plot(getindex.(coordinate_mesh,1)[:],getindex.(coordinate_mesh,2)[:], frame=:none, axes=:none )
    #Data for contouring
    mx = max(sd.bands.data...)
    mn = min(sd.bands.data...)
    band_color_ints     = [rand(3) for i=1:size(sd.bands.data,1)]
    band_color_iso     = [(min(sd.bands.data[i]...),max(sd.bands.data[i]...),band_color_ints[i]) for i=1:size(sd.bands.data,1)]
    band_color_tot     = [(mn,mx,band_color_ints[i]) for i=1:size(sd.bands.data,1)]

     Analysis.tensor_topograph(coordinate_mesh, sd.bands.data[:,:,10], coordinate_map, band_color_iso[10])
     Analysis.tensor_levelset_plot(coordinate_mesh, data, coordinate_map, band_colors, lvlset)
end


function section_function(type::Val{:stackedcuts}, programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; lvlset=0.0, kwargs...)
end

function section_function(type::Val{:topograph}, programs::Dict{String,Any}, ex_idx::Int64, s_idx::Int64; kwargs...)
    sd = Program.do_sampling(:mesh, programs, ex_idx, s_idx)
    ex = programs["executables"][ex_idx]
    #Geometric Information
    nbz                  = ex["sampling:mesh"]["inputs"][s_idx].nbz
    Λ                    = ex["sampling:mesh"]["base"].info(ex,nbz)[1][1:2,1:2]::Array{Float64,2};
    coordinate_map(x)    = (Λ*x/nbz)
    #Extracting the coordinate grids
    integer_mesh    = (ex["sampling:mesh"]["base"].mesh)([nbz+1 for _=1:size(Λ,2)], 1, 0)[1]

    coordinate_mesh = index_2_mesh.(Ref(Λ/nbz), integer_mesh)

    #Data for contouring
    mx = max(sd.bands.data...)
    mn = min(sd.bands.data...)
    band_color_iso     = [(min(sd.bands.data[i]...),max(sd.bands.data[i]...),rand(3)) for i=1:size(sd.bands.data,1)]
    band_color_tot     = [(mn,mx,rand(3)) for i=1:size(sd.bands.data,1)]

    [ Analysis.tensor_topograph(coordinate_mesh, sd.bands.data[i,:,:], coordinate_map, band_color_iso[i]) for i=1:size(sd.bands.data,1)]
end
