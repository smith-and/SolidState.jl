# Layer Density of States


#Tensor Coefficient Matrix Element
@inline function ldos_tensor_coefficient(km::KinematicOperators{Float64}, tr::HamiltonianOperators{Float64}, (a,)::Tuple{Int64}, (m,n)::NTuple{2,Int64}, dim_ℋ::Int64)::Tuple{Complex{Float64}}
    (tr.l_pjs[a][m,n],)
end

#Green's Function Element
@inline function ldos_greens_function(df::Float64, ε::Float64, δω::Float64,  (ω,δ)::NTuple{2,Float64})::Tuple{Complex{Float64}}
    (Gnωmn(1.0, ε, ω, -δ, 1),)
end

@inline function ldos_screen(dω::Array{Float64,2}, df::Vector{Array{Float64,2}}, (i1,i2)::Tuple{Int64,Int64})::Bool
    true#false
end

@inline function question(v::Val{:rank_2_prep},v2::Val{:ldos})
    true
end
