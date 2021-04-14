###Galvanic Effect

export GE
struct GE{TType <: TensorChart} <: ReponseChart
    ge::TType
end

function ge_evaluation(tc::TensorChart,K::KinematicOperators, H::HamiltonianOperators, dim_ℋ::Int64)
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
end

#=
#Irreducible Tensors
@inline function cpge_IT(dω::Array{Float64,2}, v::Array{Array{Complex{Float64},2},1}, pz::Array{Complex{Float64},2}, (a,b,c)::NTuple{3,Int64}, (l,m,n)::NTuple{3,Int64}, dim::Int64)::Complex{Float64}
    @inbounds @fastmath pz[n,l]*(v[b][l,m]*v[c][m,n] - v[c][l,m]*v[b][m,n])*(1.0/(dω[n,l]+Complex(0.0,0.02))) + pz[l,m]*(v[b][m,n]*v[c][n,l] - v[c][m,n]*v[b][n,l])*(1.0/(dω[l,m]+Complex(0.0,0.02)))
end

#Tensor Coefficient Matrix Element
@inline function cpge_rank_3_prep(km::KinematicOperators{Float64}, tr::HamiltonianOperators{Float64}, (a,b,c)::NTuple{3,Int64}, dim::Int64, (m,n,l)::NTuple{3,Int64})::Tuple{Complex{Float64}}
    if  (l==n||m==l||(-D_CUT < km.dω[n,l] < D_CUT)||(-D_CUT < km.dω[l,m] < D_CUT))
        return (Complex(0.0),)
    else
        return (Complex(0.0,-1.0)*cpge_IT(km.dω, tr.v, tr.pz, (a,b,c), (l,m,n), dim),)
    end
end

#Green's Function Element
@inline function cpge_greens_function(df::Float64, ε::Float64, dω::Float64, (ω,δ)::NTuple{2,Float64})::Tuple{Complex{Float64}}
    (Gnωmn(df,dω,ω,δ, 1)-Gnωmn(df,dω,-ω,δ, 1),)
end


=#
