##Density of States
export DOS
struct DOS{TType <: TensorChart} <: SpectralChart
    dos::TType
end

function dos_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
    for (ib,(ω,)) ∈ enumerate(tc.base)
        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
            for (ii,(i,)) ∈ enumerate(tc.indices)
                z=Complex(0.0)
                @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                for n ∈ 1:dim_ℋ
                    @fastmath ind = n + (n-1)*dim_ℋ
                    @fastmath z += 1.0/(K.dω[ind]-Complex(ω,δ))
                end
                tc.data[idx] = z
            end
        end
    end
end
