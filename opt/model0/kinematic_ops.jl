#data for set of "kinematic" obsevarble operators which are used to compute response functions
struct KinematicOperators{Float_Type <: Real}
    df::Vector{Array{Float_Type,2}}
    dω::Array{Float_Type,2}
    Δ::Array{Array{Complex{Float_Type},2},1}
    re::Array{Array{Complex{Float_Type},2},1}
    rire::Array{Array{Complex{Float_Type},2},2}
end

#Constructor for the Kinematic Operator Bundle
function KinematicOperators(dim::Int64, priorsN::Int64)
    KinematicOperators(
        [zeros(Float64,(dim,dim)) for i=1:priorsN],
        zeros(Float64,(dim,dim)),
        [zeros(ComplexF64,(dim,dim)) for _=1:2],
        [zeros(ComplexF64,(dim,dim)) for _=1:2],
        [zeros(ComplexF64,(dim,dim)) for _=1:2, _=1:2]
    )
end

struct KinematicDensity{Int_Type <: Int, Float_Type <: Real, H_Type <: HamiltonianDensity{Int_Type, Float_Type}}
    priors::Vector{Vector{Float_Type}}
    tbf::H_Type
    k_m::KinematicOperators{Float_Type}
    aux_matrix::NTuple{7,Array{Complex{Float_Type},2}}
    aux_real::NTuple{3,Array{Float_Type,2}}
    eig_aux::Tuple{Vector{Complex{Float_Type}},Vector{Float_Type},Vector{Int_Type},Array{Int_Type,2}}
    E::Eigen{Complex{Float_Type},Float_Type,Array{Complex{Float_Type},2},Array{Float_Type,1}}
end

function kinematic_prealloc(dim::Tuple{Int64,Int64})
    (
        ([zeros(Complex{Float64}, dim) for _=1:7]...,),
        ([zeros(Float64, dim) for _=1:3]...,),
        (zeros(Complex{Float64}, dim[1]), zeros(Float64, dim[1]), zeros(Int64, dim[1]), zeros(Int64,dim)) ,
        Eigen(zeros(Float64, dim[1]),zeros(ComplexF64, dim)),
    )
end

#Constructor with tbi (to be deprecated?)
function KinematicDensity(tbi::TightBindingInfo{Int64,Float64}, priors::Vector{Tuple{Symbol,Float64,Float64,Int64}})
    KinematicDensity((x->[x...]).(collect(Iterators.product([range(p[2],p[3],length=p[4]) for p∈priors]...))[:]), TightBindingDensity(tbi), KinematicOperators(size(tbi.index_rep,1),prod(getindex.(priors,4))), kinematic_prealloc(size(tbi.index_rep))...)
end

function KinematicDensity(tbi::TightBindingInfo{Int64,Float64}, priors::Array{NTuple{N,Float64},N} where N)
    KinematicDensity(((x->[x...]).(priors))[:], TightBindingDensity(tbi), KinematicOperators(size(tbi.index_rep,1),length(priors)), kinematic_prealloc(size(tbi.index_rep))...)
end

function KinematicDensity(tbi::TightBindingInfo{Int64,Float64}, priors::Vector{LinRange{Float64}})
    KinematicDensity(TightBindingDensity(tbi),priors)
end

function KinematicDensity(tbi::TightBindingInfo{Int64,Float64}, priors::Vector{Vector{Float64}})
    KinematicDensity(TightBindingDensity(tbi),priors)
end

#Constructor with Hamiltonian Density

function KinematicDensity(hd::HD where HD <: HamiltonianDensity, priors::Vector{Tuple{Symbol,Float64,Float64,Int64}})
    KinematicDensity((x->[x...]).(collect(Iterators.product([range(p[2],p[3],length=p[4]) for p∈priors]...))[:]), hd, KinematicOperators(size(hd.h_ops.h,1),prod(getindex.(priors,4))), kinematic_prealloc(size(hd.h_ops.h))...)
end

function KinematicDensity(hd::HD where HD <: HamiltonianDensity, priors::Array{NTuple{N,Float64},N} where N)
    KinematicDensity(((x->[x...]).(priors))[:], hd, KinematicOperators(size(hd.h_ops.h,1),length(priors)), kinematic_prealloc(size(hd.h_ops.h))...)
end

function KinematicDensity(hd::HD where HD <: HamiltonianDensity, priors::Vector{LinRange{Float64}})
    KinematicDensity((x->[x...]).(collect(Iterators.product(priors...))[:]), hd, KinematicOperators(size(hd.h_ops.h,1),prod(length.(priors))), kinematic_prealloc(size(hd.h_ops.h))...)
end

function KinematicDensity(hd::HD where HD <: HamiltonianDensity, priors::Vector{Vector{Float64}})
    KinematicDensity(priors, hd, KinematicOperators(size(hd.h_ops.h,1),length(priors)), kinematic_prealloc(size(hd.h_ops.h))...)
end

#Methods to obtain operator sections
function indicator(x::Float64)::Float64
    x > 0.0 ? 0.0 : 1.0
end

@inline function fermi(ϵ::Float64, T::Float64, μ::Number=0.0)::Float64
    z=0.0

    if T==0.0
        z = (ϵ - μ) > 0.0 ? 0.0 : 1.0;
    else
        z = 1.0/(1.0+exp((1.0/(8.617*(1e-5)*T))*(ϵ - μ)));
    end

    return z
end

@inline function rire_add(n::Int64, m::Int64, p::Int64, new_val::Complex{Float64}, va::Array{Complex{Float64},2}, vb::Array{Complex{Float64},2}, ra::Array{Complex{Float64},2}, rb::Array{Complex{Float64},2} )
    if !(p==n||p==m)
    @inbounds @fastmath    new_val +=(va[n,p]*rb[p,m]-rb[n,p]*va[p,m])
    end
end

@inline function km_rire(rireab::Array{Complex{Float64},2}, wab::Array{Complex{Float64},2}, va::Array{Complex{Float64},2}, vb::Array{Complex{Float64},2}, ra::Array{Complex{Float64},2}, rb::Array{Complex{Float64},2}, Δb::Array{Complex{Float64},2}, Δa::Array{Complex{Float64},2}, inv_dω::Array{Float64,2})::Nothing
    dim = 0;
    dim = size(inv_dω,1);
    @inbounds @fastmath for m ∈ 1:dim
        @inbounds @fastmath for n ∈ 1:dim
            new_rire = Complex(0.0)
            @inbounds @fastmath  for p ∈ 1:dim
                rire_add(n,m,p,new_rire,va,vb,ra,rb)
            end
            rireab[n,m] = new_rire
        end
    end

    rireab .+= ( ( (ra .* Δb) .+ (rb .* Δa) ) .+ (im .* wab) )
    rireab .*= (-1 .* inv_dω ) ;

    nothing
end

@inline function kinematic_ops(km::KinematicOperators{Float64}, cfr::HamiltonianOperators{Float64}, priors::Vector{Vector{Float64}}, inv_dω::Array{Float64,2})

    #=
    @inbounds @simd for i=eachindex(inv_dω)
        @fastmath if (inv_dω[i] > 73.50238882763689) | (inv_dω[i] < -73.50238882763689)
            inv_dω[i]=0.0
        end
    end
    =#

    @inbounds km.re[1]   .= -im .* (inv_dω) .* (cfr.v[1]);
    @inbounds km.re[2]   .= -im .* (inv_dω) .* (cfr.v[2]);

    @inbounds km_rire(km.rire[1,1], cfr.a[1,1], cfr.v[1], cfr.v[1], km.re[1], km.re[1], km.Δ[1], km.Δ[1], inv_dω)
    @inbounds km_rire(km.rire[1,2], cfr.a[1,2], cfr.v[1], cfr.v[2], km.re[1], km.re[2], km.Δ[1], km.Δ[2], inv_dω)
    @inbounds km_rire(km.rire[2,1], cfr.a[2,1], cfr.v[2], cfr.v[1], km.re[2], km.re[1], km.Δ[2], km.Δ[1], inv_dω)
    @inbounds km_rire(km.rire[2,2], cfr.a[2,2], cfr.v[2], cfr.v[2], km.re[2], km.re[2], km.Δ[2], km.Δ[2], inv_dω)

    nothing
end


@inline function kinematic_ops(km::KinematicOperators{Float64}, cfr::HamiltonianOperators{Float64}, priors::Vector{Vector{Float64}}, ϵs::Vector{Float64}, r_aux::NTuple{3,Array{Float64,2}})
    #Calculates the Fermi Functions

    #Set some diagonals
    idx=0
    @inbounds @fastmath @simd for i=1:(size(km.dω,1)+1):length(km.dω)
        idx+=1
        km.dω[i]  = ϵs[idx]
        r_aux[1][i] = 0.0
        @inbounds @fastmath for p=eachindex(priors)
            if priors[p][1]==0.0
                if (ϵs[idx] - priors[p][2]) > 0.0
                    km.df[p][i] = 0.0
                else
                    km.df[p][i] = 1.0
                end
            else
                km.df[p][i] = 1.0/(1.0+exp((1.0/(8.617*(1e-5)*priors[p][1]))*(ϵs[idx] - priors[p][2])));
            end
        end
    end

    #Set the off diagonals
    @inbounds @fastmath for i∈ eachindex(ϵs)
         @inbounds @fastmath @simd for j∈ eachindex(ϵs)
            if i != j
                idx = j+(i-1)*size(km.dω,1)
                km.dω[idx]    = ϵs[j]-ϵs[i]
                r_aux[1][idx] = 1.0 / km.dω[j,i]
                km.Δ[1][idx]  = cfr.v[1][j,j] - cfr.v[1][i,i]
                km.Δ[2][idx]  = cfr.v[2][j,j] - cfr.v[2][i,i]
                @inbounds @fastmath for p=eachindex(priors)
                    km.df[p][idx] = km.df[p][j,j] - km.df[p][i,i]
                end
            end
        end
    end

    #Assigns the rest of the operators in the section
    kinematic_ops(km, cfr, priors, r_aux[1] )

    nothing
end

function (k_obs::KinematicDensity)(k::Vector{Float64}) #::Tuple{KinematicOperators{Float64},HamiltonianOperators{Float64}}
    k_obs.tbf(k)
    cartan_frame_transform(k_obs.tbf.h_ops, k_obs.aux_matrix, k_obs.E, k_obs.eig_aux)
    kinematic_ops(k_obs.k_m, k_obs.tbf.h_ops, k_obs.priors, k_obs.E.values, k_obs.aux_real)
    nothing
end
