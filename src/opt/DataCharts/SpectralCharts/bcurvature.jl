##Berry Curvature
export BC
struct BC{TType <: TensorChart} <: SpectralChart
    bc::TType
end

function bc_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
    R   = Complex(0.0)
    for m ∈ 1:dim_ℋ
        for (ib,(ω,)) ∈ enumerate(tc.base)
            for (ip, p) ∈ enumerate(tc.priors)
                for (ii,((a,b),n)) ∈ enumerate(tc.indices)
                    if m!=n
                        @fastmath mn = m + (n-1)*dim_ℋ
                        @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                        @fastmath @inbounds tc.data[idx] += im*(transpose(H.v[a])[mn]*H.v[b][mn])*(1.0/(K.dω[mn]^2))
                    end
                end
            end
        end
    end
end
