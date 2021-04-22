###Second Harmonic Effect

export SHG
struct SHG{TType <: TensorChart} <: ReponseChart
    shg::TType
end

struct SHG_OPERATORS{Re1Type <: AbstractArray, Re2Type <: AbstractArray, Ri1Type <: AbstractArray, Ri2Type <: AbstractArray}
    Re1::Re1Type
    Re2::Re2Type
    Ri1::Ri1Type
    Ri2::Ri2Type
end

# Calculates the [re,ri] commutator matrix elements
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

# Defined in note
@inline function shg_Ri(dωmn::Float64, reanm::A, rebmn::A, recmn::A, Δbmn::A, Δcmn::A, rrbanm::A, rrabnm::A, rrcanm::A, rracnm::A, rrcbmn::A, rrbcmn::A)::Tuple{A,A} where A <: Complex{Float64}
    @fastmath Complex(0.0,-0.5/(dωmn)^2).*(
        reanm*(rebmn*Δcmn + recmn*Δbmn) + dωmn*((rrcanm*rebmn + rrbanm*recmn) - 0.5*(rrabnm*recmn + rracnm*rebmn)),
        2.0*reanm*(dωmn*(rrcbmn + rrbcmn) - 2.0*(rebmn*Δcmn + recmn*Δbmn))
    )
end

@inline function shg_Re(dωml::Float64,dωln::Float64, reanm::C, rebml::C, rebln::C, recml::C, recln::C)::C where C <: Complex{Float64}
    @fastmath return Complex((0.5/(dωln-dωml)),0.0)*reanm*(rebml*recln+recml*rebln)
end

@inline function shg_operators(χ::TensorChart, K::KinematicOperators, cfr::HamiltonianOperators, dim_ℋ::Int)

    for a ∈ 1:length(cfr.v)
        for b ∈ 1:length(cfr.v)
            @inbounds km_rire(K.rire[a,b], cfr.a[a,b], cfr.v[a], cfr.v[b], K.re[a], K.re[b], K.Δ[a], K.Δ[b], inv_dω)
        end
    end

    @fastmath @inbounds (Ri1,Ri2)   = shg_Ri(  K.dω[mn], K.re[a][nm], K.re[b][mn], K.re[c][mn], K.Δ[b][mn], K.Δ[c][mn], K.rire[b,a][nm], K.rire[a,b][nm], K.rire[c,a][nm], K.rire[a,c][nm], K.rire[c,b][mn], K.rire[b,c][mn])

    Re1 = Complex(0.0)
    Re2 = Complex(0.0)
    for l ∈ 1:dim_ℋ
        if (l!=n)&&(l!=m)
            @fastmath ml = m + (l-1)*dim_ℋ
            @fastmath ln = l + (n-1)*dim_ℋ
            @fastmath @inbounds Re2 += shg_Re(K.dω[ml],K.dω[ln],K.re[a][nm],K.re[b][ml],K.re[b][ln],K.re[c][ml],K.re[c][ln])

            @fastmath lm = l + (m-1)*dim_ℋ
            @fastmath nl = n + (l-1)*dim_ℋ
            @fastmath @inbounds Re1 += shg_Re(K.dω[lm],K.dω[mn],K.re[a][nl],K.re[b][lm],K.re[b][mn],K.re[c][lm],K.re[c][mn])
            @fastmath @inbounds Re1 += shg_Re(K.dω[mn],K.dω[nl],K.re[a][lm],K.re[b][mn],K.re[b][nl],K.re[c][mn],K.re[c][nl])
        end
    end

end

# @fastmath @inbounds χ.data[idx] += K.df[ip][mn]*((ops.Re1+ops.Ri1)*(1.0/Complex(K.dω[mn]-ω,-δ))+(ops.Re2+ops.Ri2)*(1.0/Complex(K.dω[mn]-2ω,-2δ)))
function response_summation(χ::TensorChart, T::OperatorDensity, K::KinematicOperators)
    mn = 0 ; nm = 0 ; ln = 0 ; nl = 0 ; ml = 0 ; lm = 0 ;
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                @fastmath nm = n + (m-1)*dim_ℋ
                for (ii,index) ∈ enumerate(χ.indices)
                    for (ib,base) ∈ enumerate(χ.base)
                        for (ip,priors) ∈ enumerate(χ.priors)
                            z = Complex(0.0)
                            for θ ∈ keys(T)
                                z += T[θ][nm]*G(K.df[ip][mn], K.dω[mn], priors..., base...)
                            end
                            @fastmath χ.data[ii + ((ip-1) + (ib-1)*χ.l_p)*χ.l_i] = z
                        end
                    end
                end
            end
        end
    end
end

"""
  - Have a DataCollection where all the type of Ts and Gs that need to be calculated are pooled and uniqued
    - Would require the domains also be described or just identical
        - like I mean that we might want to look at different signals over different frequencyy domains so we would compute the greens functions over different ranges
        - also I am trying to figure out if I should precompute the whole coefficient matrix or greens functions
            - probably whole coefficient matrices but not the greens functions
    - Then those are calculated and shared amongst the different terms
        - good luck
        - Might be mainly useful for the greens function reuse
    - The clear refactor that I am trying to accomplish here is that I realized that I am doing a three band summation twice
        - once in the Kinematic operators
        - once in the evalutaion
    - In addition to solving this, it present the oppurtunity to have more focussed operator evalutation and improve the more simple DataMaps
    - Looking over the previous evaluation code we see that the summation is evaluated as a three band summation nested inside a two band loop
        - probably should stick with that
"""


# 
#
# function shg_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64, ops::OperatorDensity)
#     mn = 0 ; nm = 0 ; ln = 0 ; nl = 0 ; ml = 0 ; lm = 0 ;
#     for n ∈ 1:dim_ℋ
#         for m ∈ 1:dim_ℋ
#             if m!=n
#                 @fastmath mn = m + (n-1)*dim_ℋ
#                 @fastmath nm = n + (m-1)*dim_ℋ
#                 for (ii,(a,b,c)) ∈ enumerate(tc.indices)
#                     for (ib,(ω,)) ∈ enumerate(tc.base)
#                         for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
#                             @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
#                             @fastmath @inbounds tc.data[idx] += K.df[ip][mn]*((ops.Re1+ops.Ri1)*(1.0/Complex(K.dω[mn]-ω,-δ))+(ops.Re2+ops.Ri2)*(1.0/Complex(K.dω[mn]-2ω,-2δ)))
#                         end
#                     end
#                 end
#             end
#         end
#     end
# end
