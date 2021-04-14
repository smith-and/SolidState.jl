##Density of States

export JDOS
struct JDOS{TType <: TensorChart} <: SpectralChart
    jdos::TType
end

function jdos_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
    R   = Complex(0.0)
    for m ∈ 1:dim_ℋ
        for (ip, p) ∈ enumerate(tc.priors)
            for (ib,(ω,)) ∈ enumerate(tc.base)
                for (ii,((a,b),n)) ∈ enumerate(tc.indices)
                    if m!=n
                        @fastmath mn = m + (n-1)*dim_ℋ
                        @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                        @fastmath @inbounds tc.data[idx] += K.df[ip][mn]/(K.dω[mn]-Complex(ω,δ))
                    end
                end
            end
        end
    end
end
