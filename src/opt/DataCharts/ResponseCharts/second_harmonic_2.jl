###Second Harmonic Effect

export SHG2
struct SHG2{TType <: TensorChart} <: ReponseChart
    shg2::TType
end

@inline function shg2_Ri(dωmn::Float64, reanm::A, rebmn::A, recmn::A, Δbmn::A, Δcmn::A, rrbanm::A, rrabnm::A, rrcanm::A, rracnm::A, rrcbmn::A, rrbcmn::A)::A where A <: Complex{Float64}
    @fastmath Complex(0.0,-0.5/(dωmn)^2)*2.0*reanm*(dωmn*(rrcbmn + rrbcmn) - 2.0*(rebmn*Δcmn + recmn*Δbmn))
end

@inline function shg2_Re(dωml::Float64,dωln::Float64, reanm::C, rebml::C, rebln::C, recml::C, recln::C)::C where C <: Complex{Float64}
    @fastmath return Complex((0.5/(dωln-dωml)),0.0)*reanm*(rebml*recln+recml*rebln)
end

function shg2_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
    mn = 0 ; nm = 0 ; ln = 0 ; nl = 0 ; ml = 0 ; lm = 0 ;
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                @fastmath nm = n + (m-1)*dim_ℋ
                for (ii,(a,b,c)) ∈ enumerate(tc.indices)
                    @fastmath @inbounds Ri2   = shg2_Ri(  K.dω[mn],K.re[a][nm],K.re[b][mn],K.re[c][mn],K.Δ[b][mn],K.Δ[c][mn],K.rire[b,a][nm],K.rire[a,b][nm],K.rire[c,a][nm],K.rire[a,c][nm],K.rire[c,b][mn],K.rire[b,c][mn])
                    Re2 = Complex(0.0)
                    for l ∈ 1:dim_ℋ
                        if (l!=n)&&(l!=m)
                            @fastmath ml = m + (l-1)*dim_ℋ
                            @fastmath ln = l + (n-1)*dim_ℋ
                            @fastmath @inbounds Re2 += shg2_Re(K.dω[ml],K.dω[ln],K.re[a][nm],K.re[b][ml],K.re[b][ln],K.re[c][ml],K.re[c][ln])
                        end
                    end
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += K.df[ip][mn]*((Re2+Ri2)*(1.0/Complex(K.dω[mn]-2ω,-2δ)))
                        end
                    end
                end
            end
        end
    end
end
