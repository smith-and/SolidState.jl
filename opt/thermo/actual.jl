#Grand Canonical Ensemble Quadrature Section
abstract type EnsembleData{Int_Type <: Int, Float_Type <: Real} <: QuadratureData{Int_Type, Float_Type} end

#chemical potential and Temperature are the independent variables in the base
struct GrandCanonicalData{Int_Type <: Int, Float_Type <: Real} <: EnsembleData{Int_Type, Float_Type}
    particle_number::TensorChart{Int_Type, 2, 1, Float_Type, 2,  Complex{Float_Type}, 3, Float_Type, 1}
    energy::TensorChart{        Int_Type, 2, 1, Float_Type, 2,  Complex{Float_Type}, 3, Float_Type, 1}
    free_energy::TensorChart{    Int_Type, 2, 1, Float_Type, 2,  Complex{Float_Type}, 3, Float_Type, 1}
end

function GrandCanonicalData(priors::Vector{Tuple{Symbol,Float64,Float64,Int64}})
    GrandCanonicalData(TensorChart.([("ParticleNumber", Float64, priors, [(0,0),(1,0),(0,1)]),
                                     ("Energy",         Float64, priors, [(0,0),(1,0),(0,1)]),
                                     ("FreeEnergy",     Float64, priors, [(0,0),(1,0),(0,1)])])...)
end

#filling and Temperature are the independent variables in the base
struct CanonicalData{Int_Type <: Int, Float_Type <: Real} <: EnsembleData{Int_Type, Float_Type}
    ChemicalPotential::TensorChart{Int_Type, 2, 1, Float_Type, 2,  Complex{Float_Type}, 4, Float_Type, 1}
    Energy::TensorChart{           Int_Type, 2, 1, Float_Type, 2,  Complex{Float_Type}, 4, Float_Type, 1}
    FreeEnergy::TensorChart{       Int_Type, 2, 1, Float_Type, 2,  Complex{Float_Type}, 4, Float_Type, 1}
end

function CanonicalData(priors::Vector{Tuple{Symbol,Float64,Float64,Int64}})
    CanonicalData(TensorChart.([("ChemicalPotential", Float64, priors, [(0,0),(1,0),(0,1)]),
                                ("Energy",            Float64, priors, [(0,0),(1,0),(0,1)]),
                                ("FreeEnergy",        Float64, priors, [(0,0),(1,0),(0,1)])])...)
end

#Preforms the Convex Conjugation to the Canonical Ensemble
function (CE::CanonicalData)(GCE::GrandCanonicalData)
    #identify chemical potential range
    μ_min = GCE.ParticleNumber.base[1][2]
    μ_max = GCE.ParticleNumber.base[end][2]
    #Create a numerical interpolation of the filling
    N_spline = CubicSpline(GCE.ParticleNumber.base,GCE.ParticleNumber.data)

    #identify filling range
    ν_min = CE.ChemicalPotential.base[1][2]
    ν_max = CE.ChemicalPotential.base[end][2]
    ν_pts = LinRange(ν_min,ν_max,size(CE.ChemicalPotential.base,2))

    μ_pts = zeros(Float64, length(ν_pts))
    for i ∈ 1:length(ν_pts)
        μ_pts[i] = find_roots(x->(N_spline[x]-ν_pts[i]))
    end
end

#Integration of the thermal charts
#Worker Coroutine
function gce_particle_number(ϵ, μ, T, fϵ, ρϵ, Pϵ, Δϵ, ζ = 1)
    if T > 1e-10
        return (Pϵ, -ζ*(ϵ-μ)/T^2*(1-fϵ)*Pϵ, -ζ*(1-fϵ)*Pϵ)
    elseif (-Δϵ < (ϵ-μ) < Δϵ)
        return (Pϵ, 0.0, -ζ*ρϵ/8)
    else
        return (Pϵ, 0.0, 0.0)
    end
end

function gce_energy(ϵ, μ, T, fϵ, ρϵ, Pϵ, ζ = 1)
    if T > 1e-10
        return (Pϵ*(ϵ-μ), -ζ*Pϵ*((ϵ-μ)/T)^2*(1-fϵ), -Pϵ*(1+ζ*(1-fϵ)*(ϵ-μ)/T))
    else
        return (Pϵ*(ϵ-μ), 0.0, -Pϵ)
    end
end

function gce_free_energy(ϵ, μ, T, fϵ, ρϵ, Pϵ, ζ = 1)
    if T > 1e-10
        return (ρϵ*log(1+exp(-(ϵ-μ)/T)), 0.0, 0.0 )
    else
        return (Pϵ, 0.0, 0.0)
    end
end

function gce_quad_coroutine(dos::CubicSpline, energies::Vector{Float64}, base_enum_piece::Vector{Vector{Float64}}, priors::Vector{Tuple{Symbol,Float64,Float64,Int64}}, i_size::Int64)::QuadratureData{Int64,Float64}
    #Building the Response Density Function
    gce = GrandCanonicalData(priors)
    N = length(energies)
    Δϵ = energies[2]-energies[1]

    pn = [0.0, 0.0, 0.0]
    e  = [0.0, 0.0, 0.0]
    fe = [0.0, 0.0, 0.0]

    fϵ, ρϵ, Pϵ = 0.0, 0.0, 0.0
    #Starting mesh work
    for (b_i, (T,μ)) ∈ base_enum_piece
        pn .= 0.0; e  .= 0.0; fe .= 0.0
        for ϵ ∈ energies
            fϵ = fermi(ϵ,T,μ)
            ρϵ = dos[ϵ]
            Pϵ = fϵ*ρϵ

            pn .+= gce_particle_number(ϵ, μ, T, fϵ, ρϵ, Pϵ, Δϵ)
            e  .+= gce_energy(         ϵ, μ, T, fϵ, ρϵ, Pϵ)
            fe .+= gce_free_energy(    ϵ, μ, T, fϵ, ρϵ, Pϵ)
        end

        setindex!(gce.particle_number.data, pn/N, range(1+(b_i-1)*i_size,length=i_size) )
        setindex!(gce.energy.data,          e/N,  range(1+(b_i-1)*i_size,length=i_size) )
        setindex!(gce.free_energy.data,     fe/N, range(1+(b_i-1)*i_size,length=i_size) )
    end

    gce_data
end

function gce_quad_prep(ex, id, nw)
    priors = ex[:priors]::Vector{Tuple{Symbol,Float64,Float64,Int64}}

    base_pieces = mesh_chunks(collect(Iterators.Enumerate(Iterators.product(TensorCharts.prior_axis.(priors)...))), nw)

    dos_chart  = TensorChart(ex["quad:Spectral"]["result"][id]*"/jdos")
    dos_splines = Vector{CubicSpline}(undef, size(dos_chart.base,2))

    for δ_idx = size(dos_chart.base,2):size(dos_chart.base,2)
        dos_splines[δ_idx] = CubicSpline(dos_chart.base[:,δ_idx],dos_chart.data[1,1,:,δ_idx])
    end

    (priors, dos_splines, base_pieces)
end

function quadrature(v::Val{:GCE}, QDataType::Symbol, ex::Dict{String,Any}; id::Int64,  nbz::Int64)
    w_pool = CachingPool(workers());
    nw = length(w_pool);

    (priors, dos_splines, base_pieces)  = gce_quad_prep(ex, id, nw)
    #Making Remote Calls and Storing the Futures
    worker_data = Vector{Future}(undef, nw)
    for i=1:nw
        worker_data[i] = remotecall(gce_quad_coroutine, w_pool, dos_splines[end], energy_mesh, base_pieces, priors, 3)
    end
    #Aggregate the Worker Results and output data
    return quad_agregate(fetch.(worker_data))
end

#Quadrature routines
include("gce.jl")

function do_single_quad(val::Val{:GCE}, qd::DataType, ex::Dict{String,Any}, q_idx::Int64)::String
    q_tag = "quad:"*string(quad_type)
    ex[q_tag]["done?"][q_idx] ? ex[q_tag]["result"][q_idx] :
    begin
        kwargs = merge(ex[q_tag]["inputs"][q_idx],(id = q_idx,))
        stats = begin @timed quadrature(val, quad_type, ex; kwargs...)::DosData{Int64,Float64} end
        ex[q_tag]["metric"][q_idx]   = (t=stats.time,bytes=stats.bytes,gc=stats.gctime,mem=stats.gcstats)
        v =  Val(quad_type)
        #ex[q_tag]["plot"][q_idx]     = quad_plot( v, quad_type, stats.value, ex, q_idx)
        #ex[q_tag]["result"][q_idx]   = quad_export(v, quad_type, stats.value, ex, q_idx)
        ex[q_tag]["done?"][q_idx]    = true;
        println("did quadrature for "*string((ex[q_tag]["inputs"][q_idx].nbz::Int64)^2)*" points in +"*string(stats.time)*"s")
        #ex[q_tag]["result"][q_idx]
        "hi"
    end
end

#batch quad call
function do_quad(val::Val{:GCE}, quad_type::Symbol, program::Dict{String,Any})::Nothing
    for ex ∈ program["executables"]
        for q_idx ∈ eachindex(ex["quad:"*string(quad_type)]["inputs"])
            do_single_quad(quad_type, ex,q_idx)
        end
    end
    nothing
end

#sinlge quad call
function do_quad(val::Val{:GCE}, quad_type::Symbol, program::Dict{String,Any},     ex_idx::Int64, qscope_idx::Int64)::String
    do_single_quad(quad_type, program["executables"][ex_idx], qscope_idx)
end

function one_dimensional_ex(q_i)
   Dict{String,Any}(
        "base"        => (integrator=:Energy),
        "inputs"      => q_i,
        "metric"      => Array{ NamedTuple{(:t, :bytes, :gc, :mem),Tuple{Float64,Int64,Float64,Base.GC_Diff}},length(size(q_i))}(undef,size(q_i)),
        "metric:plot" => Dict{Symbol,AbstractPlot}(),
        "result"      => Array{String,length(size(q_i))}(undef,size(q_i)),
        "plot"        => Array{Dict{Symbol,AbstractPlot},length(size(q_i))}(undef,size(q_i)),
        "done?"       => fill(false,size(q_i))
    )
end
