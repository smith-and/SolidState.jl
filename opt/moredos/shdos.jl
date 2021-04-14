#Joint Density of States
#Tensor Coefficient Matrix Element
@inline function shdos_rank_2_prep(km::KinematicOperators{Float64}, tr::HamiltonianOperators{Float64}, (a,)::NTuple{1,Int64}, dim::Int64, (m,n)::NTuple{2,Int64})::NTuple{2,Complex{Float64}}
    return (Complex(1.0),Complex(1.0))
end

#Green's Function Element
@inline function shdos_greens_function(df::Float64, ε::Float64, δω::Float64, (ω,δ)::NTuple{2,Float64})::NTuple{2,Complex{Float64}}
    (Gnωmn(df, δω, -ω, -δ, 1), Gnωmn(-df/2.0, δω, -ω, -δ, 2))
end

@inline function shdos_screen(dω::Array{Float64,2}, df::Vector{Array{Float64,2}}, (i1,i2)::Tuple{Int64,Int64})::Bool
    (i1==i2)
end

@inline function question(v::Val{:rank_2_prep},v2::Val{:shdos})
    true
end
