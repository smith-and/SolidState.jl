#data for set of "kinematic" obsevarble operators which are used to compute response functions
struct KinematicOperators{AR <: AbstractArray, AC <: AbstractArray}
    df::Vector{AR}
    dω::AR
    Δ::Array{AC,1}
    re::Array{AC,1}
    rire::Array{AC,2}
    pz::AC
    pz_0::AC
end

#Constructor for the Kinematic Operator Bundle
function KinematicOperators(dim::Integer, priorsN::Integer,d=2; style::Symbol=:normal)
    if style == :shared
        return KinematicOperators(
            [SharedArray(zeros(Float64,(dim,dim))) for i=1:priorsN],
            SharedArray(zeros(Float64,(dim,dim))),
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d, _=1:d],
            SharedArray(zeros(ComplexF64,(dim,dim))),
            SharedArray(zeros(ComplexF64,(dim,dim)))
        )
    elseif style == :static
        return KinematicOperators(
            [SMatrix{dim,dim}(zeros(Float64,(dim,dim))) for i=1:priorsN],
            SMatrix{dim,dim}(zeros(Float64,(dim,dim))),
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d, _=1:d],
            SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))),
            SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim)))
        )
    else
        return KinematicOperators(
            [(zeros(Float64,(dim,dim))) for i=1:priorsN],
            (zeros(Float64,(dim,dim))),
            [(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [(zeros(ComplexF64,(dim,dim))) for _=1:d, _=1:d],
            (zeros(ComplexF64,(dim,dim))),
            (zeros(ComplexF64,(dim,dim)))
        )
    end
end


struct KinematicDensity{H_Type <: HamiltonianDensity, K_Type <: KinematicOperators, AUXA <: AbstractArray, PA <: AbstractArray}
    hd::H_Type
    k_m::K_Type
    aux_real::AUXA
    priors::PA
end

# Constructor with density
function KinematicDensity(hd::H where H <: HamiltonianDensity{Int64,F}, priors::Vector{Tuple{Symbol,F,F,Int64}}; style::Symbol=:normal) where F <: Number
    inputs = ( hd, KinematicOperators(size(hd.h_ops.h,1),prod(getindex.(priors,4)); style=style), zeros(F, size(hd.h_ops.h)),(x->[x...]).(collect(Iterators.product([range(p[2],p[3],length=p[4]) for p∈priors]...))[:]) )
    KinematicDensity{typeof.(inputs)...}(inputs...)
end

function KinematicDensity(hd::H where H <: HamiltonianDensity{Int64,F}, priors::Vector{LinRange{F}}; style::Symbol=:normal) where F <: AbstractFloat
    inputs = (hd, KinematicOperators(size(hd.h_ops.h,1),prod(length.(priors)); style=style), zeros(F, size(hd.h_ops.h)),(x->[x...]).(collect(Iterators.product(priors...))[:]))
    KinematicDensity{typeof.(inputs)...}(inputs...)
end

function KinematicDensity(hd::H where H <: HamiltonianDensity{Int64,F}, priors::Vector{Vector{F}}; style::Symbol=:normal) where F <: AbstractFloat
    inputs = (hd, KinematicOperators(size(hd.h_ops.h,1),length(priors); style=style), zeros(typeof(priors[1][1]), size(hd.h_ops.h)),priors)
    KinematicDensity{typeof.(inputs)...}(inputs...)
end

function KinematicDensity(hd::HD where HD <: HamiltonianDensity, priors::Array{NTuple{N,Float64},N} where N; style::Symbol=:normal)
    KinematicDensity(hd, KinematicOperators(size(hd.h_ops.h,1),length(priors); style=style), zeros(typeof(0.0),size(hd.h_ops.h)),((x->[x...]).(priors))[:])
end


# Fermi-Dirac Functions
function indicator(x::F)::F where F <: AbstractFloat
    x > 0.0 ? 0.0 : 1.0
end

@inline function fermi(ϵ::F, T::F, μ::Number=0.0)::F where F <: AbstractFloat
    z=0.0

    if T==0.0
        z = (ϵ - μ) > 0.0 ? 0.0 : 1.0;
    else
        z = 1.0/(1.0+exp((1.0/(8.617*(1e-5)*T))*(ϵ - μ)));
    end

    return z
end


@inline function rire_add(n::Int64, m::Int64, p::Int64, new_val::Complex{F}, va::AbstractArray{Complex{F},2}, vb::AbstractArray{Complex{F},2}, ra::AbstractArray{Complex{F},2}, rb::AbstractArray{Complex{F},2} ) where F <: AbstractFloat
    if !(p==n||p==m)
    @inbounds @fastmath    new_val +=(va[n,p]*rb[p,m]-rb[n,p]*va[p,m])
    end
end

@inline function km_rire(rireab::AbstractArray{Complex{F},2}, wab::AbstractArray{Complex{F},2}, va::AbstractArray{Complex{F},2}, vb::AbstractArray{Complex{F},2}, ra::AbstractArray{Complex{F},2}, rb::AbstractArray{Complex{F},2}, Δb::AbstractArray{Complex{F},2}, Δa::AbstractArray{Complex{F},2}, inv_dω::AbstractArray{F,2})::Nothing where F <: AbstractFloat
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

@inline function kinematic_ops(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, inv_dω::AbstractArray{F,2}) where F <: AbstractFloat

    for a ∈ 1:length(cfr.v)
        @inbounds km.re[a]   .= -im .* (inv_dω) .* (cfr.v[a]);
    end

    for a ∈ 1:length(cfr.v)
        for b ∈ 1:length(cfr.v)
            @inbounds km_rire(km.rire[a,b], cfr.a[a,b], cfr.v[a], cfr.v[b], km.re[a], km.re[b], km.Δ[a], km.Δ[b], inv_dω)
        end
    end

    nothing
end


@inline function kinematic_ops(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, ϵs::Vector{F}, r_aux::AbstractArray{F,2}) where F <: AbstractFloat
    #Calculates the Fermi Functions

    #Set some diagonals
    idx=0
    @inbounds @fastmath @simd for i=1:(size(km.dω,1)+1):length(km.dω)
        idx+=1
        km.dω[i]  = ϵs[idx]
        r_aux[i] = 0.0
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
                r_aux[idx] = 1.0 / km.dω[j,i]
                km.Δ[1][idx]  = cfr.v[1][j,j] - cfr.v[1][i,i]
                km.Δ[2][idx]  = cfr.v[2][j,j] - cfr.v[2][i,i]
                @inbounds @fastmath for p=eachindex(priors)
                    km.df[p][idx] = km.df[p][j,j] - km.df[p][i,i]
                end
            end
        end
    end

    #Assigns the rest of the operators in the section
    kinematic_ops(km, cfr, priors, r_aux)

    nothing
end

function (k_obs::KinematicDensity)(k::Vector{F}) where F <: Number #::Tuple{KinematicOperators{AbstractFloat},HamiltonianOperators}
    k_obs.hd(k)
    cartan_transform(k_obs.hd.h_ops)
    kinematic_ops(k_obs.k_m, k_obs.hd.h_ops, k_obs.priors, k_obs.hd.h_ops.E.values, k_obs.aux_real)
    nothing
end

#=
#Constructor for Block priors
function KinematicDensity(tbi::TightBindingInfo{Int64,F}, priors::Array{NTuple{N,F},N} where N) where F <: AbstractFloat
    KinematicDensity(((x->[x...]).(priors))[:], TightBindingDensity(tbi), KinematicOperators(size(tbi.index_rep,1),length(priors)), SharedArray(zeros(F, size(tbi.index_rep))))
end

function KinematicDensity(tbi::TightBindingInfo{Int64,F}, priors::Vector{Tuple{Symbol,F,F,Int64}}) where F <: Number
    inputs = ( (x->[x...]).(collect(Iterators.product([range(p[2],p[3],length=p[4]) for p∈priors]...))[:]), TightBindingDensity(tbi), KinematicOperators(size(tbi.index_rep,1),prod(getindex.(priors,4))), SharedArray(zeros(F, size(tbi.index_rep))))
    signature = (typeof(1),typeof(1.0),typeof(inputs[2]),typeof(inputs[3]))
    KinematicDensity{signature...}(inputs...)
end

function KinematicDensity(tbi::TightBindingInfo{Int64,F}, priors::Vector{LinRange{F}}) where F <: AbstractFloat
    Kinematic(TightBindingDensity(tbi),priors)
end

function KinematicDensity(tbi::TightBindingInfo{Int64,F}, priors::Vector{Vector{F}}) where F <: AbstractFloat
    KinematicDensity(TightBindingDensity(tbi),priors)
end
=#
#=
@inbounds @simd for i=eachindex(inv_dω)
    @fastmath if (inv_dω[i] > cutoff) | (inv_dω[i] < -cutoff)
        inv_dω[i]=0.0
    end
end
=#
