

@inline function rire_add(n::Int64, m::Int64, p::Int64, new_val::Complex{F}, va::AbstractArray{Complex{F},2}, vb::AbstractArray{Complex{F},2}, ra::AbstractArray{Complex{F},2}, rb::AbstractArray{Complex{F},2} ) where F <: AbstractFloat
    if !(p==n||p==m)
    @inbounds @fastmath    new_val +=(va[n,p]*rb[p,m]-rb[n,p]*va[p,m])
    end
end

@inline function km_rire(rireab::AbstractArray{Complex{F},2}, wab::AbstractArray{Complex{F},2}, va::AbstractArray{Complex{F},2}, vb::AbstractArray{Complex{F},2}, ra::AbstractArray{Complex{F},2}, rb::AbstractArray{Complex{F},2}, Δb::AbstractArray{Complex{F},2}, Δa::AbstractArray{Complex{F},2}, inv_dω::AbstractArray{F,2})::Nothing where F <: AbstractFloat
    dim = 0;
    dim = size(inv_dω,1);
    @inbounds @fastmath for m ∈ 1:dim
        @inbounds @fastmath for n ∈ 1:dim
            new_rire = Complex(0.0)
            @inbounds @fastmath  for p ∈ 1:dim
                rire_add(n,m,p,new_rire,va,vb,ra,rb)
            end
            rireab[n,m] = new_rire
        end
    end

    rireab .+= ( ( (ra .* Δb) .+ (rb .* Δa) ) .+ (im .* wab) )
    rireab .*= (-1 .* inv_dω ) ;

    nothing
end

@inline function kinematic_ops(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, inv_dω::AbstractArray{F,2}) where F <: AbstractFloat

    for a ∈ 1:length(cfr.v)
        @inbounds km.re[a]   .= -im .* (inv_dω) .* (cfr.v[a]);
    end

    for a ∈ 1:length(cfr.v)
        for b ∈ 1:length(cfr.v)
            @inbounds km_rire(km.rire[a,b], cfr.a[a,b], cfr.v[a], cfr.v[b], km.re[a], km.re[b], km.Δ[a], km.Δ[b], inv_dω)
        end
    end

    nothing
end


@inline function kinematic_ops(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, ϵs::Vector{F}, r_aux::AbstractArray{F,2}) where F <: AbstractFloat
    #Calculates the Fermi Functions

    #Set some diagonals
    idx=0
    @inbounds @fastmath @simd for i=1:(size(km.dω,1)+1):length(km.dω)
        idx+=1
        km.dω[i]  = ϵs[idx]
        r_aux[i] = 0.0
        @inbounds @fastmath for p=eachindex(priors)
            if priors[p][1]==0.0
                if (ϵs[idx] - priors[p][2]) > 0.0
                    km.df[p][i] = 0.0
                else
                    km.df[p][i] = 1.0
                end
            else
                km.df[p][i] = 1.0/(1.0+exp((1.0/(8.617*(1e-5)*priors[p][1]))*(ϵs[idx] - priors[p][2])));
            end
        end
    end

    #Set the off diagonals
    @inbounds @fastmath for i∈ eachindex(ϵs)
         @inbounds @fastmath @simd for j∈ eachindex(ϵs)
            if i != j
                idx = j+(i-1)*size(km.dω,1)
                km.dω[idx]    = ϵs[j]-ϵs[i]
                r_aux[idx] = 1.0 / km.dω[j,i]
                km.Δ[1][idx]  = cfr.v[1][j,j] - cfr.v[1][i,i]
                km.Δ[2][idx]  = cfr.v[2][j,j] - cfr.v[2][i,i]
                @inbounds @fastmath for p=eachindex(priors)
                    km.df[p][idx] = km.df[p][j,j] - km.df[p][i,i]
                end
            end
        end
    end

    #Assigns the rest of the operators in the section
    kinematic_ops(km, cfr, priors, r_aux)

    nothing
end

function (k_obs::KinematicDensity)(k::Vector{F}) where F <: Number #::Tuple{KinematicOperators{AbstractFloat},HamiltonianOperators}
    k_obs.hd(k)
    cartan_transform(k_obs.hd.h_ops)
    kinematic_ops(k_obs.k_m, k_obs.hd.h_ops, k_obs.priors, k_obs.hd.h_ops.E.values, k_obs.aux_real)
    nothing
end


###Second Harmonic Effect

export SHG
struct SHG{TType <: TensorChart} <: ReponseChart
    shg::TType
end

@inline function shg_Ri(dωmn::Float64, reanm::A, rebmn::A, recmn::A, Δbmn::A, Δcmn::A, rrbanm::A, rrabnm::A, rrcanm::A, rracnm::A, rrcbmn::A, rrbcmn::A)::Tuple{A,A} where A <: Complex{Float64}
    @fastmath Complex(0.0,-0.5/(dωmn)^2).*(
        reanm*(rebmn*Δcmn + recmn*Δbmn) + dωmn*((rrcanm*rebmn + rrbanm*recmn) - 0.5*(rrabnm*recmn + rracnm*rebmn)),
        2.0*reanm*(dωmn*(rrcbmn + rrbcmn) - 2.0*(rebmn*Δcmn + recmn*Δbmn))
    )
end

@inline function shg_Re(dωml::Float64,dωln::Float64, reanm::C, rebml::C, rebln::C, recml::C, recln::C)::C where C <: Complex{Float64}
    @fastmath return Complex((0.5/(dωln-dωml)),0.0)*reanm*(rebml*recln+recml*rebln)
end

function shg_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
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
end
