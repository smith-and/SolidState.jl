
###############################
### Linear Conductivity 
###############################

const ϵ = 55.26349406*(2π)^3/1e6 # 1/(MV⋅m)

ϵ2 = 1/42.106021882608665

export SHG
abstract type SHG <: OpticsChart end

@inline function shg_Ri(dωmn::Float64, reanm::A, rebmn::A, recmn::A, Δbmn::A, Δcmn::A, rrbanm::A, rrabnm::A, rrcanm::A, rracnm::A, rrcbmn::A, rrbcmn::A)::Tuple{A,A} where A <: Complex{Float64}
    @fastmath (
        -im*0.5/(dωmn)^2*(reanm*(rebmn*Δcmn + recmn*Δbmn) + dωmn*((rracnm*rebmn + rrabnm*recmn) - 0.5*(rrbanm*recmn + rrcanm*rebmn))),
        -im/(dωmn)^2*reanm*(dωmn*(rrcbmn + rrbcmn) - 2.0*(rebmn*Δcmn + recmn*Δbmn))
    )
end

@inline function shg_Re(dωml::Float64,dωln::Float64, reanm::C, rebml::C, rebln::C, recml::C, recln::C)::C where C <: Complex{Float64}
    @fastmath return Complex((0.5/(dωln-dωml)),0.0)*reanm*(rebml*recln+recml*rebln)
end

function SHG_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
    mn = 0 ; nm = 0 ; ln = 0 ; nl = 0 ; ml = 0 ; lm = 0 ;
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                @fastmath nm = n + (m-1)*dim_ℋ
                for (ii,(a,b,c)) ∈ enumerate(tc.indices)
                    @fastmath @inbounds (Ri1,Ri2)   = shg_Ri(  K.dω[mn],K.re[a][nm],K.re[b][mn],K.re[c][mn],K.Δ[b][mn],K.Δ[c][mn],K.rire[b,a][nm],K.rire[a,b][nm],K.rire[c,a][nm],K.rire[a,c][nm],K.rire[c,b][mn],K.rire[b,c][mn])
                    Re1 = Complex(0.0)
                    Re2 = Complex(0.0)
                    for l ∈ 1:dim_ℋ
                        if (l!=n)&&(l!=m)#&&(-1e-4 < K.dω[l,n] < 1e-4)&&(-1e-4 < K.dω[l,n] < 1e-4)
                            @fastmath ml = m + (l-1)*dim_ℋ
                            @fastmath ln = l + (n-1)*dim_ℋ
                            @fastmath @inbounds Re2 += -2.0*shg_Re(K.dω[ml],K.dω[ln],K.re[a][nm],K.re[b][ml],K.re[b][ln],K.re[c][ml],K.re[c][ln])

                            @fastmath lm = l + (m-1)*dim_ℋ
                            @fastmath nl = n + (l-1)*dim_ℋ
                            @fastmath @inbounds Re1 += shg_Re(K.dω[lm],K.dω[mn],K.re[a][nl],K.re[b][lm],K.re[b][mn],K.re[c][lm],K.re[c][mn])
                            @fastmath @inbounds Re1 += shg_Re(K.dω[mn],K.dω[nl],K.re[a][lm],K.re[b][mn],K.re[b][nl],K.re[c][mn],K.re[c][nl])
                        end
                    end
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += K.df[ip][mn]*((Re1+Ri1)*(1.0/Complex(K.dω[mn]-ω,-δ))+(Re2+Ri2)*(1.0/Complex(K.dω[mn]-2ω,-2δ)))
                        end
                    end
                end
            end
        end
    end
    tc.data.= tc.data.*(36.474728077745624)*2
    copy(tc.data)
end
