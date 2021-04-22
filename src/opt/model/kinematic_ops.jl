
### We define some ways to collect the priors

function collect_priors(priors::Vector{Tuple{Symbol,F,F,Int64}})
    (x->[x...]).(collect(Iterators.product([range(p[2],p[3],length=p[4]) for p∈priors]...))[:])
end

function collect_priors(priors::Array{NTuple{N,Float64},N} where N)
    ((x->[x...]).(priors))[:]
end

function collect_priors(priors::Vector{LinRange{F}})
    (x->[x...]).(collect(Iterators.product(priors...))[:])
end

function collect_priors(priors::Vector{Vector{F}})
    priors
end

### Now We define some operators

abstract type OperatorSet end
abstract type OperatorDensity end

@generated function (od::OperatorDensity)(hd::HamiltonianDensity) #::Tuple{KinematicOperators{AbstractFloat},HamiltonianOperators}
    name = fieldname(od,1)
    quote
        $(Symbol(name,:(_ops)))(od.$name, hd.h_ops, od.priors, hd.h_ops.E.values, od.aux_real)
        nothing
    end
end


#data for set of "kinematic" obsevarble operators which are used to compute response functions
struct KinematicOperators{AR <: AbstractArray, AC <: AbstractArray} <: OperatorSet
    df::Vector{AR}
    dω::AR
    Δ::Array{AC,1}
    re::Array{AC,1}
end

#Constructor for the Kinematic Operator Bundle
function KinematicOperators(dim::Integer, priorsN::Integer,d=2; style::Symbol=:normal)
    if style == :shared
        return KinematicOperators(
            [SharedArray(zeros(Float64,(dim,dim))) for i=1:priorsN],
            SharedArray(zeros(Float64,(dim,dim))),
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d],
        )
    elseif style == :static
        return KinematicOperators(
            [SMatrix{dim,dim}(zeros(Float64,(dim,dim))) for i=1:priorsN],
            SMatrix{dim,dim}(zeros(Float64,(dim,dim))),
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d],
        )
    else
        return KinematicOperators(
            [(zeros(Float64,(dim,dim))) for i=1:priorsN],
            (zeros(Float64,(dim,dim))),
            [(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [(zeros(ComplexF64,(dim,dim))) for _=1:d],
        )
    end
end


# Fermi-Dirac Functions
function indicator(x::F)::F where F <: AbstractFloat
    x > 0.0 ? 0.0 : 1.0
end

@inline function fermi_bose(ϵ::F, T::F, μ::Number=0.0, ζ = 1)::F where F <: AbstractFloat
    z=0.0

    if T==0.0
        z = (ϵ - μ) > 0.0 ? 0.0 : 1.0;
    else
        z = 1.0/(1.0+ζ*exp((1.0/(8.617*(1e-5)*T))*(ϵ - μ)));
    end

    return z
end

@inline function k_m_ops(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, ϵs::Vector{F}, r_aux::AbstractArray{F,2}) where F <: AbstractFloat
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
    for a ∈ 1:length(cfr.v)
        @inbounds km.re[a]   .= -im .* (r_aux) .* (cfr.v[a]);
    end

    nothing
end

struct KinematicDensity{H_Type <: HamiltonianDensity, O_Type <: OperatorSet, AUXA <: AbstractArray, PA <: AbstractArray} <: OperatorDensity
    k_m::O_Type
    aux_real::AUXA
    priors::PA
end

function KinematicDensity(hd::H where H <: HamiltonianDensity{Int64,F}, priors0; style::Symbol=:normal) where F <: AbstractFloat
    priors = collect_priors(priors0)
    inputs = (KinematicOperators(size(hd.h_ops.h,1),length(priors); style=style), zeros(F, size(hd.h_ops.h)),priors)
    KinematicDensity{typeof.(inputs)...}(inputs...)
end
