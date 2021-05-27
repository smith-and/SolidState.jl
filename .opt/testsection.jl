#export TESTCHART
struct TESTCHART{TType <: TensorChart} <: DataChart
    test::TType
end

function test_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
    for (ib,b) ∈ enumerate(tc.base)
        for (ip,p) ∈ enumerate(tc.priors)
            for (ii,i) ∈ enumerate(tc.indices)
                idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                tc.data[idx] = b[1]
            end
        end
    end
end
