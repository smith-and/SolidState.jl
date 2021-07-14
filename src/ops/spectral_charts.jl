"""
   abstract type SpectralChart <: DataChart end
"""
abstract type SpectralChart <: DataChart end

###############################
### Volume Form
###############################

export VOL
struct VOL{TType <: TensorChart}
    vol::TType
end

function vol_evaluation(tc::TensorChart, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
    for (ii,(a,b,c)) ∈ enumerate(tc.indices)
        for (ib,(ω,)) ∈ enumerate(tc.base)
            for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                @fastmath @inbounds tc.data[idx] = Complex(1.0)
            end
        end
    end
    SArray(tc.data)
end

#################################
#### Density of States
#################################

export DOS
struct DOS{TType <: TensorChart} <: SpectralChart
    dos::TType
end

function dos_evaluation(tc::TensorChart, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
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
    copy(tc.data)
end

#################################
#### Joint Density of States
#################################

export JDOS
struct JDOS{TType <: TensorChart} <: SpectralChart
    jdos::TType
end

function jdos_evaluation(tc::TensorChart, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
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
    copy(tc.data)
end

#####################################
#### Off Diagonal Berry Connection
#####################################

export ODBC
struct ODBC{TType <: TensorChart} <: SpectralChart
    odbc::TType
end

function odbc_evaluation(tc::TensorChart, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
    m = Int(dim_ℋ/2)
    n = m+1
    @fastmath mn = m + (n-1)*dim_ℋ
    @fastmath nm = n + (m-1)*dim_ℋ
    for (ib,(ω,)) ∈ enumerate(tc.base)
        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
            for (ii,(a,b,c)) ∈ enumerate(tc.indices)
                @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                tc.data[idx] = K.re[a][nm]*(K.re[b][mn]*K.Δ[c][mn]+K.re[c][mn]*K.Δ[b][mn])
            end
        end
    end
    copy(tc.data)
end

#################################
#### Berry Curvature
#################################
export BC
struct BC{TType <: TensorChart} <: SpectralChart
    bc::TType
end

function bc_evaluation(tc::TensorChart, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
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
    copy(tc.data)
end
