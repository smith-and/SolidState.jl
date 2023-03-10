abstract type OpticsChart  <: ChartType end

###############################
### Linear Permitivity
###############################

export LP
abstract type LP <: OpticsChart end

function LP_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    tc.data .= 0.0
    R   = Complex(0.0)
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                @fastmath nm = n + (m-1)*dim_ℋ
                for (ii,(a,b)) ∈ enumerate(tc.indices)
                    @fastmath @inbounds R = (K.re[a][nm]*K.re[b][mn])
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += 0.729494556*R*(K.df[ip][nm]/(K.dω[mn]-Complex(ω,δ)))
                        end
                    end
                end
            end
        end
    end
    copy(tc.data)
end


###############################
### Linear Permitivity
###############################

export CC2L
abstract type CC2L<: OpticsChart end

function CC2L_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    SolidState.conjugate_projectors!(P,K0.hd)
    K = K0.k_m
    H = K0.hd.h_ops
    tc.data .= 0.0
    R   = Complex(0.0)
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                @fastmath nm = n + (m-1)*dim_ℋ
                for (ii,(a,b)) ∈ enumerate(tc.indices)
                    @fastmath @inbounds R = real((P.proj[1]*H.v[a]*P.proj[1])[nm]*(P.proj[2]*H.v[b]*P.proj[2])[mn])
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += 4*im*ω*R*(K.df[ip][nm]/(K.dω[mn])^2/(K.dω[mn]-Complex(ω,δ)))
                        end
                    end
                end
            end
        end
    end
    copy(tc.data)
end

###############################
###Galvanic Effect
###############################

export GE
abstract type GE <: OpticsChart end

function GE_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    K0(k)
    K = K0.k_m
    H = K0.hd.h_ops
    tc.data .= 0.0
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                for (ii,(a,b,c)) ∈ enumerate(tc.indices)
                    R   = Complex(0.0)
                    for l ∈ 1:dim_ℋ
                        if (l!=n)&&(l!=m)
                            @fastmath nl = n + (l-1)*dim_ℋ
                            @fastmath lm = l + (m-1)*dim_ℋ
                            @fastmath @inbounds R += K.pz[nl]*(H.v[c][lm]*H.v[b][mn]-H.v[b][lm]*H.v[c][mn])/(Complex(K.dω[nl],0.02))
                            @fastmath @inbounds R += K.pz[lm]*(H.v[c][mn]*H.v[b][nl]-H.v[b][mn]*H.v[c][nl])/(Complex(K.dω[lm],0.02))
                        end
                    end
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += R*K.df[ip][mn]*(1.0/Complex(K.dω[mn]-ω,-δ)-1.0/Complex(K.dω[mn]+ω,-δ))
                        end
                    end
                end
            end
        end
    end
    copy(tc.data)
end

###############################
###Second Harmonic Effect
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

###############################
###Second Harmonic Effect Full
###############################

const ϵ = 55.26349406*(2π)^3/1e6 # 1/(MV⋅m)

ϵ2 = 1/42.106021882608665

export SHGF
abstract type SHGF <: OpticsChart end

@inline function shgf_Ri(dωmn::Float64, reanm::A, rebmn::A, recmn::A, Δbmn::A, Δcmn::A, rrbanm::A, rrabnm::A, rrcanm::A, rracnm::A, rrcbmn::A, rrbcmn::A)::Tuple{A,A} where A <: Complex{Float64}
    @fastmath (
        -im*0.5/(dωmn)^2*(reanm*(rebmn*Δcmn + recmn*Δbmn) + dωmn*((rracnm*rebmn + rrabnm*recmn) - 0.5*(rrbanm*recmn + rrcanm*rebmn))),
        -im/(dωmn)^2*reanm*(dωmn*(rrcbmn + rrbcmn) - 2.0*(rebmn*Δcmn + recmn*Δbmn))
    )
end

@inline function shgf_Re(dωml::Float64,dωln::Float64, reanm::C, rebml::C, rebln::C, recml::C, recln::C)::C where C <: Complex{Float64}
    @fastmath return Complex((0.5/(dωln-dωml)),0.0)*reanm*(rebml*recln+recml*rebln)
end

function SHGF_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
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

###############################
### Second Harmonic Effect
###############################

const ϵ = 55.26349406*(2π)^3/1e6 # 1/(MV⋅m)

ϵ2 = 1/42.106021882608665

export LPSHG
abstract type LPSHG <: OpticsChart end

@inline function lpshg_Ri(dωmn::Float64, reanm::A, rebmn::A, recmn::A, Δbmn::A, Δcmn::A, rrbanm::A, rrabnm::A, rrcanm::A, rracnm::A, rrcbmn::A, rrbcmn::A)::Tuple{A,A} where A <: Complex{Float64}
    @fastmath (
        -im*0.5/(dωmn)^2*(reanm*(rebmn*Δcmn + recmn*Δbmn) + dωmn*((rracnm*rebmn + rrabnm*recmn) - 0.5*(rrbanm*recmn + rrcanm*rebmn))),
        -im/(dωmn)^2*reanm*(dωmn*(rrcbmn + rrbcmn) - 2.0*(rebmn*Δcmn + recmn*Δbmn))
    )
end

@inline function lpshg_Re(dωml::Float64,dωln::Float64, reanm::C, rebml::C, rebln::C, recml::C, recln::C)::C where C <: Complex{Float64}
    @fastmath return Complex((0.5/(dωln-dωml)),0.0)*reanm*(rebml*recln+recml*rebln)
end

function LPSHG_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
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
                for (ii,((l1,l2),(a,b,c))) ∈ enumerate(tc.indices)
                    ra = P.proj[l1]*K.re[a]*P.proj[l2]
                    @fastmath @inbounds (Ri1,Ri2)   = lpshg_Ri(  K.dω[mn],ra[nm],K.re[b][mn],K.re[c][mn],K.Δ[b][mn],K.Δ[c][mn],K.rire[b,a][nm],K.rire[a,b][nm],K.rire[c,a][nm],K.rire[a,c][nm],K.rire[c,b][mn],K.rire[b,c][mn])
                    Re1 = Complex(0.0)
                    Re2 = Complex(0.0)
                    for l ∈ 1:dim_ℋ
                        if (l!=n)&&(l!=m)#&&(-1e-4 < K.dω[l,n] < 1e-4)&&(-1e-4 < K.dω[l,n] < 1e-4)
                            @fastmath ml = m + (l-1)*dim_ℋ
                            @fastmath ln = l + (n-1)*dim_ℋ
                            @fastmath @inbounds Re2 += -2.0*lpshg_Re(K.dω[ml],K.dω[ln],ra[nm],K.re[b][ml],K.re[b][ln],K.re[c][ml],K.re[c][ln])

                            @fastmath lm = l + (m-1)*dim_ℋ
                            @fastmath nl = n + (l-1)*dim_ℋ
                            @fastmath @inbounds Re1 += lpshg_Re(K.dω[lm],K.dω[mn],K.re[a][nl],K.re[b][lm],K.re[b][mn],K.re[c][lm],K.re[c][mn])
                            @fastmath @inbounds Re1 += lpshg_Re(K.dω[mn],K.dω[nl],K.re[a][lm],K.re[b][mn],K.re[b][nl],K.re[c][mn],K.re[c][nl])
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

###############################
### Static Limit of SHG
###############################

export StaticSHG
abstract type StaticSHG <: OpticsChart end


function static_shg_direct(re,Δ,w,δω,b,c,mn)
    @fastmath @inbounds -im*(re[b][mn]*Δ[c][mn]+re[c][mn]*Δ[b][mn]+im*(w[b,c][mn]+ w[c,b][mn])./2)./(2 * δω[mn].^3)
end

function static_shg_direct(H,K,(b,c),mn)
    static_shg_direct(K.re,K.Δ,H.a,K.dω,b,c,mn)
end

function static_shg_virtual_coef(dω,ll,nm,ml,ln)
    @fastmath @inbounds dω[ll]/dω[nm]*dω[ml]*dω[ln] - (dω[ml]-dω[ln])/(2dω[nm]^3)
end

function static_shg_virtual(re,dω,(b,c),(m,n),dim)
    z = Complex(0.0)
    nm = n+(m-1)*dim
    for l in dim
        ml = m+(l-1)*dim
        ln = l+(n-1)*dim
        ll = l+(l-1)*dim
        z += static_shg_virtual_coef(dω,ll,nm,ml,ln)*(re[b][ml]*re[c][ln]+re[c][ml]*re[b][ln])/2
    end
    z
end

function StaticSHG_evaluation(tc::DataChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
    tc.data .= 0.0
    mn = 0 ; nm = 0 ; ln = 0 ; nl = 0 ; ml = 0 ; lm = 0 ;
    for n ∈ 1:dim_ℋ
        for m ∈ 1:dim_ℋ
            if m!=n
                @fastmath mn = m + (n-1)*dim_ℋ
                @fastmath nm = n + (m-1)*dim_ℋ
                for (ii,(a,b,c)) ∈ enumerate(tc.indices)
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += K.df[ip][nm]*K.re[a][nm]*(static_shg_direct(H,K,(b,c),mn)+static_shg_virtual(K.re,K.dω,(b,c),(m,n),dim_ℋ))
                        end
                    end
                end
            end
        end
    end
    copy(tc.data)
end


###############################
###Second Harmonic Effect Better
###############################

const ϵ = 55.26349406*(2π)^3/1e6 # 1/(MV⋅m)

abstract type SHG0 <: OpticsChart end

@inline function shg_Ri0(dωmn::Float64, reanm::A, rebmn::A, recmn::A, Δbmn::A, Δcmn::A, rrbanm::A, rrabnm::A, rrcanm::A, rracnm::A, rrcbmn::A, rrbcmn::A)::Tuple{A,A} where A <: Complex{Float64}
    @fastmath Complex(0.0,-0.5/(dωmn)^2).*(
        reanm*(rebmn*Δcmn + recmn*Δbmn) + dωmn*((rrcanm*rebmn + rrbanm*recmn) - 0.5*(rrabnm*recmn + rracnm*rebmn)),
        2.0*reanm*(dωmn*(rrcbmn + rrbcmn) - 2.0*(rebmn*Δcmn + recmn*Δbmn))
    )
end

@inline function shg_Re0(dωml::Float64,dωln::Float64, reanm::C, rebml::C, rebln::C, recml::C, recln::C)::C where C <: Complex{Float64}
    @fastmath return Complex((0.5/(dωln-dωml)),0.0)*reanm*(rebml*recln+recml*rebln)
end

function SHG0_evaluation(tc::DataChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
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
                        if (l!=n)&&(l!=m)
                            @fastmath ml = m + (l-1)*dim_ℋ
                            @fastmath ln = l + (n-1)*dim_ℋ
                            @fastmath @inbounds Re2 += shg_Re(K.dω[ml],K.dω[ln],K.re[a][nm],K.re[b][ml],K.re[b][ln],K.re[c][ml],K.re[c][ln])


                            @fastmath @inbounds Re1 += shg_Re(K.dω[lm],K.dω[mn],K.re[a][nl],K.re[b][lm],K.re[b][mn],K.re[c][lm],K.re[c][mn])
                            @fastmath @inbounds Re1 += shg_Re(K.dω[mn],K.dω[nl],K.re[a][lm],K.re[b][mn],K.re[b][nl],K.re[c][mn],K.re[c][nl])
                        end
                    end
                    for (ib,(ω,)) ∈ enumerate(tc.base)
                        for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                            @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                            @fastmath @inbounds tc.data[idx] += K.df[ip][mn]/ϵ*((Re1+Ri1)*(1.0/Complex(K.dω[mn]-ω,-δ))+(Re2+Ri2)*(1.0/Complex(K.dω[mn]-2ω,-2δ)))
                        end
                    end
                end
            end
        end
    end
    # tc.data.= tc.data./ϵ
    copy(tc.data)
end
