export SampleData
abstract type SampleData{Int_Type <: Int, Float_Type <: Real} end
abstract type BANDS end
#Data type to hold information about model that we want visualized over the BZ
struct SpectralData{Int_Type <: Int, Float_Type <: Real} <: SampleData{Int_Type, Float_Type}
    bands::TensorChart{      Int_Type, 1, 1, Float64, 1, Int_Type, 2, Float_Type, 3, 1}
    #band_grads::TensorChart{ Int_Type, 2, 2, Int_Type, 2, Complex{Float_Type}, 4, Float64, 1}
    #inv_mass::TensorChart{   Int_Type, 3, 3, Int_Type, 2, Complex{Float_Type}, 5, Float64, 1}
    #structure::TensorChart{  Int_Type, 2, 2, Int_Type, 2, Float_Type, 4,          Float64, 1}
end

#Zeroed Construtor
function SpectralData(dim::Int64, h_dim::Int64, mesh_ranges::Vector{Tuple{Symbol,Int64,Int64,Int64}})
    SpectralData(TensorChart.((
        ("bands"  ,    Float64, mesh_ranges, [(i,) for i=1:h_dim]),
        #("b_grads" ,   Complex{Float64}, mesh_ranges, collect(Iterators.product(1:h_dim,1:dim))),
        #("inv_mass" ,  Complex{Float64}, mesh_ranges, collect(Iterators.product(1:h_dim,1:dim,1:dim))),
        #("structure" , Float64, mesh_ranges, collect(Iterators.product(1:h_dim,1:h_dim)))
    ), Ref(1))...)
end

#Addition method on Sample Data to facilitate aggregation among workers
function Base.:(+)(sd1::SpectralData,sd2::SpectralData)::SpectralData
    sd1.bands.data      .= sd1.bands.data       .+ sd2.bands.data
    #sd1.band_grads.data .= sd1.band_grads.data  .+ sd2.band_grads.data
    #sd1.inv_mass.data   .= sd1.inv_mass.data    .+ sd2.inv_mass.data;
    #sd1.structure.data  .= sd1.structure.data   .+ sd2.structure.data;
    sd1
end

struct SpectralDensity{Int_Type <: Int, Float_Type <: Real, H_Type <: HamiltonianDensity{Int_Type, Float_Type}}
    k_obs::KinematicDensity{Int_Type,Float_Type, H_Type}
    data::SpectralData{Int_Type,Float_Type}
end

#External Constructor for SpectralDensity
function SpectralDensity(tbi::TightBindingInfo{Int64,Float64}, mesh_ranges::Vector{Tuple{Symbol,Int64,Int64,Int64}})
    SpectralDensity(KinematicDensity(tbi, [LinRange(0.0,0.0,1)]), SpectralData(2, size(tbi.index_rep,1), mesh_ranges))
end

#Sample Density Function routines
#Coroutine to store the energy eigenvalues
function band_function(pt::Array{Int64}, chart::TensorChart, energies::Vector{Float64})::Nothing
    setindex!(chart.data, energies,:,pt[1],pt[2])
    nothing
end

#Coroutine to store the band gradients
function band_gradient_function(pt::Array{Int64}, chart::TensorChart,vx::Vector{Complex{Float64}},vy::Vector{Complex{Float64}})::Nothing
    setindex!(chart.data,vx,:,1,pt[1],pt[2]);
    setindex!(chart.data,vy,:,2,pt[1],pt[2]);
    nothing
end

#Coroutine to store the band weights for each Bloch state
function structure_function(pt::Array{Int64}, chart::TensorChart, U::Array{Float64,2})::Nothing
    setindex!(chart.data,U,:,:,pt[1],pt[2])
    nothing
end

#Coroutines to calculate and store the inverse mass tensor
function inv_mass(α::Int64, dω::Array{Float64,2}, rej::Array{ComplexF64,2}, rei::Array{ComplexF64,2})::Complex{Float64}
    z = Complex(0.0)
    for β ∈ 1:size(dω,1)
        if β==α
            nothing
        else
            z +=  dω[β,α] * (rej[α,β] * rei[β,α] + rei[α,β] * rej[β,α])
        end
    end
    return z
end

function inv_mass(α::Int64, km::KinematicOperators{Float64}, tr::HamiltonianOperators{Float64},i::Int64,j::Int64)::Complex{Float64}
    tr.a[i,j][α,α] - inv_mass(α, km.dω, km.re[j], km.re[i])
end

function inv_mass_function(dim::Int64, pt::Array{Int64}, chart::TensorChart{Int64,3,3,Int64,2,Complex{Float64},5}, km::KinematicOperators{Float64}, tr::HamiltonianOperators{Float64})::Nothing
    setindex!(chart.data, inv_mass.(1:dim, Ref(km), Ref(tr), 1, 1), :,1,1,pt[1],pt[2])
    setindex!(chart.data, inv_mass.(1:dim, Ref(km), Ref(tr), 1, 2), :,1,2,pt[1],pt[2])
    setindex!(chart.data, inv_mass.(1:dim, Ref(km), Ref(tr), 2, 1), :,2,1,pt[1],pt[2])#setindex!(chart.data, getindex(chart.data,pt[1],pt[2],α,1,2),pt[1],pt[2],α,2,1)
    setindex!(chart.data, inv_mass.(1:dim, Ref(km), Ref(tr), 2, 2), :,2,2,pt[1],pt[2])
    nothing
end

#Function to call the different coroutines
function sample_density_function( sample::SpectralData{Int64,Float64}, km::KinematicOperators{Float64}, tr::HamiltonianOperators{Float64}, pt::Array{Int64,1}, dim::Int64)::Nothing
    band_function(pt,          sample.bands,             real.(diag(tr.h)))
    #band_gradient_function(pt, sample.band_grads,        diag(tr.v[1]), diag(tr.v[2]))
    #structure_function(pt,     sample.structure,         tr.structure)
    #inv_mass_function(dim, pt, sample.inv_mass,          km, tr)
    nothing
end

#Method to calculate the sampling sections
#path
function (s_d::SpectralDensity)(k_pt::Vector{Float64}, pt::Vector{Int64})::Nothing
    s_d.k_obs(k_pt)
    sample_density_function(s_d.data,  s_d.k_obs.k_m, s_d.k_obs.tbf.h_ops, pt, size(s_d.k_obs.aux_matrix,1))
    nothing
end

#mesh
function (s_d::SpectralDensity)(Λ::Array{Float64,2}, pt::Vector{Int64})::Nothing
    s_d.k_obs(Λ*pt)
    sample_density_function(s_d.data, s_d.k_obs.k_m, s_d.k_obs.tbf.h_ops , pt+[1,1], size(s_d.k_obs.aux_matrix,1))
    nothing
end

#Routines to do worker dispatch for sampling over the base values

#Coroutines
function sampling_coroutine(tbi::TightBindingInfo{Int64,Float64}, mesh_ranges::Vector{Tuple{Symbol,Int64,Int64,Int64}}, Λdk::Array{Vector{Float64},2}, mesh_piece::Vector{Vector{Int64}})::SpectralData{Int64,Float64}
    #Building the Response Density Function
    sdf = SpectralDensity(tbi, mesh_ranges)
    #Starting mesh work
    for pt ∈ mesh_piece
        sdf(Λdk[pt...], pt)
    end
    sdf.data
end

function sampling_agregate(worker_data::Vector{SpectralData{Int64,Float64}})::SpectralData{Int64,Float64}
    #Building the Sampling Data for Agregation
    sdf_data = worker_data[1];
    for data ∈ worker_data[2:end]
        sdf_data::SpectralData{Int64,Float64} = sdf_data + data
    end
    return sdf_data
end

#Executable Translations
function r_sampling_info(ex,nbz)
    asdg = ASDGeometry(ex["asd"])
    Λk=hcat(asdg["bz_hs"]["K3"],asdg["bz_hs"]["K1"])
    dA=abs(det(Λk[1:2,1:2]))/(2*9*nbz^2);
    ( Λk, dA, asdg["bz_hs"])
end

function sampling_prep(v::Val{:path}, ex::Dict{String,Any}, n_path::Int64, nw::Int64)::(Tuple{Array{Vector{Float64},2}, Vector{Tuple{Symbol, Int64, Int64, Int64}}, Vector{Vector{Vector{Int64}}},TightBindingInfo{Int64,Float64}})
    (Λ,dA,hs_bz_leger) = ex["sampling:path"]["base"].info(ex, n_path)::Tuple{Array{Float64,2},Float64, Dict{String,Vector{Float64}}};
    paths              = path_points(hs_bz_leger, ex["input"]["paths"], n_path)[1]
    path_ranges        = [(:n, 1, n_path, n_path),(:p, 1, length(paths), length(paths))]
    p_mesh, path_pieces        = (ex["sampling:path"]["base"].mesh)([n_path, length(ex["input"]["paths"])], nw)
    path_indices       = [(x->[x...]).(path_piece) for path_piece ∈ path_pieces]
    #Model Information
    tbi=ex["model"]::TightBindingInfo{Int64,Float64}

    (paths, path_ranges , path_indices, tbi)
end

function sampling_prep(type::Val{:mesh}, ex::Dict{String,Any}, nbz::Int64, nw::Int64)::(Tuple{Array{Vector{Float64},2}, Vector{Tuple{Symbol, Int64, Int64, Int64}}, Vector{Vector{Vector{Int64}}},TightBindingInfo{Int64,Float64}})
    ((v1_name,v2_name), shift) = ex["input"]["mesh"]
    asdg = ASDGeometry(ex["asd"])
    Λ=hcat(asdg["bz_hs"][v1_name],asdg["bz_hs"][v2_name])

    (mesh0, mesh_pieces) = (ex["sampling:mesh"]["base"].mesh)([nbz for _=1:size(Λ,2)], nw)
    mesh_indices        = [(x->[x...]).(mesh_piece) for mesh_piece ∈ mesh_pieces]
    mesh_ranges         = [(Meta.parse("n"*string(i)), 0, nbz, nbz+1) for i=1:size(Λ,2)]

    mesh                = (index_to_mesh(Λ/(nbz),shift)).(mesh0)::Array{Vector{Float64},2}

    #Model Information
    tbi=ex["model"]::TightBindingInfo{Int64,Float64}

    (mesh, mesh_ranges, mesh_indices, tbi)
end

#Public Interface
function sampling(type::Symbol, ex::Dict{String,Any}; nbz::Int64)::SpectralData{Int64,Float64}
    #Inputs for the Coroutine Call
    #Making the Mesh
    w_pool = ex["utilities"]["workers"];
    nw = length(w_pool);

    (type == :path) && ((samp_info, samp_ranges, samp_pts, tbi) = sampling_prep(Val(:path), ex, nbz, nw));
    (type == :mesh) && ((samp_info, samp_ranges, samp_pts, tbi) = sampling_prep(Val(:mesh), ex, nbz, nw));


    #Making Remote Calls and Storing the Futures
    worker_futures = Vector{Future}(undef, nw)
    for i=1:nw
        worker_futures[i] = remotecall(sampling_coroutine, w_pool, tbi, samp_ranges, samp_info, samp_pts[i])
    end;

    #Obtaining the data from the workers
    worker_results::(Vector{SpectralData{Int64,Float64}}) = fill(SpectralData(2, size(tbi.index_rep,2), samp_ranges), nw)::Vector{SpectralData{Int64,Float64}}
    for i=1:length(worker_results)
        worker_results[i] = fetch(worker_futures[i])::SpectralData{Int64,Float64}
    end

    #Aggregating the Worker Results and output data
    return sampling_agregate(worker_results)
end

#Functions for the Analysis Executable
#Sampling routines
function do_single_sampling(type::Symbol, ex::Dict{String,Any}, s_idx::Int64)::Nothing
    if !ex["sampling:"*string(type)]["done?"][s_idx]
        s_tag = "sampling:"*string(type)

        stats  = @timed sampling(type, ex; ex[s_tag]["inputs"][s_idx]...)
        stats2 = @timed (ex["utilities"]["plot?"] ? sampling_plot(Val(type), type, stats.value, ex, s_idx) : nothing )
        stats3 = @timed (ex["utilities"]["save?"] ? sampling_export(         type, stats.value, ex, s_idx) : nothing )

        ex["bandmin"]   = min(stats.value.bands.data[:]...)
        ex["bandmax"]   = max(stats.value.bands.data[:]...)
        ex["bandwidth"] = ex["bandmax"]-ex["bandmin"]

        for i=1:length(ex["input"]["Analysis.DosData{Int64,Float64}"])
            for j=1:length(ex["input"]["Analysis.DosData{Int64,Float64}"][i][3])
                if ex["input"]["Analysis.DosData{Int64,Float64}"][i][3][j][1]==:ω
                    ex["input"]["Analysis.DosData{Int64,Float64}"][i][3][j] = (:ω,ex["bandmin"]*1.1,ex["bandmax"]*1.1,ex["input"]["Analysis.DosData{Int64,Float64}"][i][3][j][4])
                end
            end
        end

        for i=1:length(ex["input"]["Analysis.MoreDosData{Int64,Float64}"])
            for j=1:length(ex["input"]["Analysis.MoreDosData{Int64,Float64}"][i][3])
                if ex["input"]["Analysis.MoreDosData{Int64,Float64}"][i][3][j][1]==:ω
                    ex["input"]["Analysis.MoreDosData{Int64,Float64}"][i][3][j] = (:ω,0.0,3.1, ex["input"]["Analysis.MoreDosData{Int64,Float64}"][i][3][j][4])
                end
            end
        end

        for i=1:length(ex["input"]["Analysis.OpticalData{Int64,Float64}"])
            for j=1:length(ex["input"]["Analysis.OpticalData{Int64,Float64}"][i][3])
                if ex["input"]["Analysis.OpticalData{Int64,Float64}"][i][3][j][1]==:ω
                    ex["input"]["Analysis.OpticalData{Int64,Float64}"][i][3][j] = (:ω,0.0,3.1,ex["input"]["Analysis.OpticalData{Int64,Float64}"][i][3][j][4])
                end
            end
        end


        #Writing and printing metrics
        ex[s_tag]["metric"][s_idx] = Base.structdiff(stats,(value=1,))
        ex[s_tag]["done?"][s_idx]  = true

        ptn = string(type==:path ? (ex[s_tag]["inputs"][s_idx].nbz*length(ex["input"]["paths"])) : (ex[s_tag]["inputs"][s_idx].nbz)^2*length(ex["paths"]))
        println("Analysis.SamplingPath"*"\t\t"*string(ptn)*"pts\t +"*string(round(stats.time;digits=3))*"s\t +"*string(round(stats2.time;digits=3))*"s\t +"*string(round(stats3.time;digits=3))*"s")
        flush(stdout)

        nothing
    end
end

#single executable scope call
function do_sampling(type::Symbol, program::Dict{String,Any}, ex_idx::Int64)::Nothing
    let ex = program["executables"][ex_idx]
        merge!( ex, Dict("sampling:$(string(type))"=>sampling_ex( ex["scope"]["$(string(type))"])))
        for s_idx ∈ eachindex(ex["scope"]["$(string(type))"])
            do_single_sampling(type, ex, s_idx)
        end
    end
end

#batch routine call
function do_sampling(type::Symbol, program::Dict{String,Any})::Nothing
    for ex_id ∈ 1:length(program["executables"])
        do_sampling(type,program,ex_id)
    end
    nothing
end

function sampling_ex(s_i)
   Dict{String,Any}(
        "base"        => (n=3, mesh=natural_mesh, info=r_sampling_info),
        "inputs"      => s_i,
        "metric"      => fill(Base.structdiff((@timed 1+1), (value=1.0,)),size(s_i)),
        "metric:plot" => Dict{Symbol,String}(),
        "result"      => fill("",size(s_i)),
        "plot"        => fill(Dict{Symbol,String}(),size(s_i)),
        "done?"       => fill(false,size(s_i))
    )
end
