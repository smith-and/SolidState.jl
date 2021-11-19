abstract type SpectralChart <: ChartType end

#################################
#### Bands
#################################

export BANDS
abstract type BANDS <: SpectralChart end

function BANDS_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
    for (ib,k) ∈ enumerate(tc.base)
        for (ip,(μ,)) ∈ enumerate(tc.priors)
            H(k)
            tc.data[:,:,ip,ib] = hcat(eigen(Hermitian(H.h_ops.h))...)
        end
    end
    tc.data
end

#################################
#### Density of States
#################################

export DOS
abstract type DOS <: SpectralChart end

function DOS_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
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
#### Projected Density of States
#################################

export PDOS
abstract type PDOS <: SpectralChart end
function PDOS_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
    K = K0.k_m
    H = K0.hd
    conjugate_projectors!(P,H)
    tc.data .= 0.0
    for (ip, (T,μ,δ)) ∈ enumerate(tc.priors)
        for (ib,(ω,)) ∈ enumerate(tc.base)
            for (ii,(l1,l2)) ∈ enumerate(tc.indices)
                @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                @fastmath G = (1.0 ./ (diag(K.dω).-Complex(ω,δ)))
                @fastmath @inbounds tc.data[idx] =  tr(P.proj[l1]*Diagonal(G)*P.proj[l2])
            end
        end
    end
    copy(tc.data)
end

#################################
#### Joint Density of States
#################################

export JDOS
abstract type JDOS <: SpectralChart end
function JDOS_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
    R   = Complex(0.0)
    for m ∈ 1:dim_ℋ
        for n ∈ 1:dim_ℋ
            for (ip, (T,μ,δ)) ∈ enumerate(tc.priors)
                for (ib,(ω,)) ∈ enumerate(tc.base)
                    for (ii,(a,)) ∈ enumerate(tc.indices)
                        if m!=n
                            @fastmath mn = m + (n-1)*dim_ℋ
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += -K.df[ip][mn]/(K.dω[mn]-Complex(ω,δ))
                        end
                    end
                end
            end
        end
    end
    copy(tc.data)
end

#####################################
#### Projected Density of States
#####################################

export PJDOS
abstract type PJDOS <: SpectralChart end
function PJDOS_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
    K = K0.k_m
    H = K0.hd
    conjugate_projectors!(P,H)
    tc.data .= 0.0
    for (ip, (T,μ,δ)) ∈ enumerate(tc.priors)
        for (ib,(ω,)) ∈ enumerate(tc.base)
            for (ii,(l1,l2)) ∈ enumerate(tc.indices)
                @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                @fastmath G = (-(K.df[ip].* P.OD) ./ (K.dω.-Complex(ω,δ)))
                @fastmath @inbounds tc.data[idx] =  sum(P.proj[l1]*G*P.proj[l2])
            end
        end
    end
    copy(tc.data)
end

#####################################
#### Off Diagonal Berry Connection
#####################################

export ODBC

abstract type ODBC <: SpectralChart end
function ODBC_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
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

abstract type BC <: SpectralChart end
function BC_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
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
