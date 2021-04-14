##Linear Polarization Tensor

export LP
struct LP{TType <: TensorChart} <: ReponseChart
    lp::TType
end

function lp_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
    R   = Complex(0.0)
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                for (ii,(a,b)) ∈ enumerate(tc.indices)
                    @fastmath @inbounds R = (transpose(H.v[a])[mn]*H.v[b][mn])*(1.0/(K.dω[mn]))
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += im*R*(K.df[ip][mn]/(K.dω[mn]-Complex(ω,δ)))
                        end
                    end
                end
            end
        end
    end
end
