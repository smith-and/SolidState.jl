#Piecewise-Linear Path (1D) Discretization
function path_points(path_corners::Vector{Vector{Float64}}, N0::Int64)::Array{Array{Float64,1},1}
    N=N0-1
    #obtaining the difference vectors between subsequent path points
    relative_vecs::Vector{Vector{Float64}}=[path_corners[i+1]-path_corners[i] for i=1:(length(path_corners)-1)]
    #Obtaining the number of points in each leg for approximately N overall
    p_lengths = LinearAlgebra.norm.(relative_vecs)
    p_length  = sum(p_lengths)
    leg_Ns = max.(floor.(p_lengths ./ p_length .* N ),Ref(1))
    i=0; while N-sum(leg_Ns) > 0
        leg_Ns[i%length(leg_Ns)+1] += 1
    end
    #Calculating the points along each leg of the path
    path_positions = push!(vcat([ [path_corners[leg_idx] + relative_vecs[leg_idx]*pt_idx/leg_Ns[leg_idx] for pt_idx=0:(leg_Ns[leg_idx]-1) ] for leg_idx=1:length(leg_Ns)]...),path_corners[end])

    path_positions
end

function findindex(pathlist,point)
    for i=1:length(pathlist)
        LinearAlgebra.norm(pathlist[i]-point) < 1e-10 ? (return i; break) : nothing
    end
    return 0 #this is bad
end

function path_points(path_list::Vector{Vector{Vector{Float64}}}, N::Int64)::(Tuple{Array{Vector{Float64},2},Array{Array{Int64,1},1},Array{Array{Float64,1},1}})
    paths           = Vector{Vector{Float64}}[];
    odometers       = Vector{Float64}[];
    corner_indices  = [zeros(Int64, length(path_list[i])) for i=1:length(path_list)];
    for i ∈ 1:length(path_list)
        positions = path_points(path_list[i], N)
        push!(paths,    positions)
        push!(odometers, vcat([0.0],cumsum(LinearAlgebra.norm.(diff(positions)))))
        #println(indexin.(hash.(path_list[i]), Ref(hash.(positions))))
        corner_indices[i] .= findindex.(Ref(positions), path_list[i])
    end

    hcat(paths...), corner_indices, odometers
end

function path_points(bz_pt_ledger::Dict{String,Array{Float64,1}}, path_list, N::Int64)::(Tuple{Array{Vector{Float64},2},Array{Array{Int64,1},1},Array{Array{Float64,1},1}})
    path_points([get.(Ref(bz_pt_ledger), path_corners, Ref([Inf,Inf])) for path_corners ∈ path_list], N)
end

function index_to_mesh(Λ::Array{Float64,2}, w::Vector{Float64})
    function mesh_coordinate(x::NTuple{N,Int64} where N)::Vector{Float64}
        Λ*[x...] .+ w
    end
    return mesh_coordinate
end

## Surface (2D) Meshes

#determine chunk sizes
function chunk_sizes(nw::Int64,dim::Int64)
    dims=[];
    for i=1:nw push!(dims,Int(floor(dim/nw)) + ((i-1)<(dim%nw) ? 1 : 0 ) ) end;
    dims
end

function chunk_sizes(nw::Int64,nbz::Vector{Int64},shift=0)
    dim = prod(nbz.+shift)
    dims=[];
    for i=1:nw push!(dims,Int(floor(dim/nw)) + ((i-1)<(dim%nw) ? 1 : 0 ) ) end;
    dims
end

#break up mesh into nw pieces
function mesh_chunks(mesh,nw)
    flat_mesh = mesh[:]
    mesh_fill = deepcopy(flat_mesh[1])
    pieces    = [fill( mesh_fill, dim) for dim ∈ chunk_sizes(nw,length(flat_mesh))];
    for i=1:length(pieces)
        for j=1:length(pieces[i])
            pieces[i][j] = pop!(flat_mesh)
        end
    end
    @assert length(flat_mesh)==0
    return pieces
end

#Hyper Cubic Lattice in first quadrant up to ndims values without zeros
function natural_mesh(ndims::Vector{Int64},nw::Int64)
    r_mesh    = collect(Iterators.product([1:ndims[i] for i=1:length(ndims)]...))
    return r_mesh, mesh_chunks(r_mesh, nw)
end

#Mesh for trigonal C3 system

function c3_star_patches(Λ::Array{Float64,2})::Vector{Array{Float64,2}}
    [Λ,R3D(1,2π/3)*Λ,R3D(1,4π/3)*Λ]
end

function c3_sector_mesh(ndims::Int64,nw::Int64)
    nbz = ndims[1]
    r_mesh    = collect(Iterators.product([0:nbz, 1:nbz]...))

    mesh_pts = mesh_chunks(r_mesh, nw)
    push!(mesh_pts[end], (0,0) )

    return r_mesh, mesh_pts
end

#Triangular Lattice Weight Maps for Trapezoidal Quadrature
function tt_wmap(N::Int64,dA::Float64)
    function wmap((i1,i2)::NTuple{2,Int64})::Float64
        if (i1<N) && (i2<N) && !(i1==i2==0)#bulk
            return 6.0/6.0*dA
        elseif ((i1==N) && (0<i2<N)) || ((i2==N) && (0<i1<N)) #edge
            return 3.0/6.0*dA
        elseif (i1==0 && i2==0) || ((i1==N) && (i2==N))  || ((i1==0) && (i2==N)) || ((i1==N) && (i2==0)) #corner
            return 2.0/6.0*dA
        else
            println((i1,i2))
            error("nope")
            return Inf
        end
    end
end

function c3_simplex_patch(Λ::Array{Float64,2})::Vector{Array{Float64,2}}
    [Λ]
end

function c3_sector_mesh_simplex(ndims::Int64,nw::Int64)
    nbz = ndims[1]
    r_mesh    = collect(Iterators.product([0:nbz, 0:nbz]...))

    mesh_pts = mesh_chunks(r_mesh, nw)
    push!(mesh_pts[end], (0,0) )

    return r_mesh, mesh_pts
end

#Triangular Lattice Weight Maps for Trapezoidal Quadrature
function tt_wmap_simplex(N::Int64,dA::Float64)
    function wmap((i1,i2)::NTuple{2,Int64})::Float64
        if (i1<N) && (i2<N) && !(i1==0|i2==0)#bulk
            return 6.0/6.0*dA
        elseif ((i1==N) && (0<i2<N)) || ((i2==N) && (0<i1<N)) #edge
            return 3.0/6.0*dA
        elseif (i1==0 && i2==0) || ((i1==N) && (i2==N))
            return 2.0/6.0*dA
        elseif ((i1==0) && (i2==N)) || ((i1==N) && (i2==0)) #corner
            return dA/6.0
        else
            println((i1,i2))
            error("nope")
            return Inf
        end
    end
end

#Lebesgue (Monte Carlo) integration
function lebesgue_mesh(ndims::Int64, nw::Int64)
    cdims = chunk_sizes(nw,ndims^2) .* 3
    pieces = [fill((0.0,0.0), dim) for dim ∈ cdims];

    int_C3 = ([ 1.0  0.0 ; 0.0  1.0 ],
              [ 0.0 -1.0 ; 1.0 -1.0 ],
              [ -1.0 1.0 ; -1.0 0.0]
    )

    for i=1:length(pieces)
        for j=1:length(pieces[i])
            pieces[i][j] = (int_C3[(j%3) + 1])*rand(2) |> x->(x...,)
        end
    end

    return [(0,0)], pieces
end

#Trivial Weight Map
function lebesgue_wmap(N::Int64,dA::Float64)
    function wmap(idx::NTuple{2,Float64})::Float64
        return dA
    end
end

#Mesh Visuals
weightcolor(w) = begin
    if w==6
        return :green
    elseif w==3
        return :yellow
    elseif w==2
        return :red
    else
        return :black
    end
end

function mesh_check(Λdks, mesh_piece, weights )
    x=Λdks[1]*[0.0,0.0]

    mesh_plot  = plot(size=(1500,1500), framestyle=:box, aspectratio=2/sqrt(3))
    area_check = 0.0

    for i=1:length(mesh_piece)
        #loop over patches to make the Voronio unit cell
        for j=1:length(Λdks)
            x .= Λdks[j]*mesh_piece[i]
            scatter!([x[1]],[x[2]], color=weightcolor(weights[i]),label="", m=10, opacity=.1)
            area_check += weights[i]
            display(mesh_plot)
        end
    end;
    println("the area was ", string(area_check))
end
