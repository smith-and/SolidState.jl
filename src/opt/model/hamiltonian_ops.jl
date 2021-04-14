abstract type HamiltonianDensity{Int_Type <: Integer, Float_Type <: Real} end

struct HamiltonianOperators{A <: AbstractArray, Float_Type <: AbstractFloat, Int_Type <: Integer}
    h::A
    v::Array{A,1}
    a::Array{A,2}
    aux_matrix::A
    eig_aux::Tuple{Vector{Complex{Float_Type}},Vector{Float_Type},Vector{Int_Type},Array{Int_Type,2}}
    E::Eigen{Complex{Float_Type},Float_Type,Array{Complex{Float_Type},2},Array{Float_Type,1}}
end

#Easy constructor for section data
function HamiltonianOperators(dim, d=2; style::Symbol=:normal)::HamiltonianOperators

    if style == :shared
        return HamiltonianOperators(SharedArray(zeros(ComplexF64,(dim,dim))),
                            [SharedArray(zeros(ComplexF64,(dim,dim)))   for _=1:d ],
                            [SharedArray(zeros(ComplexF64,(dim,dim)))   for _=1:d, _=1:d],
                            SharedArray(zeros(ComplexF64, (dim,dim))),
                            (zeros(ComplexF64, dim), zeros(Float64, dim), zeros(Int64, dim), zeros(Int64,(dim,dim))) ,
                            Eigen(zeros(Float64, dim),zeros(ComplexF64, (dim,dim)))
                            )
    elseif style == :static
        return HamiltonianOperators(SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))),
                            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim)))   for _=1:d ],
                            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim)))   for _=1:d, _=1:d],
                            SMatrix{dim,dim}(zeros(ComplexF64, (dim,dim))),
                            (SVector{dim}(zeros(ComplexF64, dim)), SVector{dim}(zeros(Float64, dim)), SVector{dim}(zeros(Int64, dim)), SMatrix{dim,dim}(zeros(Int64,(dim,dim)))) ,
                            Eigen(SVector{dim}(zeros(Float64, dim)),SMatrix{dim,dim}(zeros(ComplexF64, (dim,dim))))
                            )
    else
        return HamiltonianOperators((zeros(ComplexF64,(dim,dim))),
                            [(zeros(ComplexF64,(dim,dim)))   for _=1:d ],
                            [(zeros(ComplexF64,(dim,dim)))   for _=1:d, _=1:d],
                            (zeros(ComplexF64, (dim,dim))),
                            (zeros(ComplexF64, dim), zeros(Float64, dim), zeros(Int64, dim), zeros(Int64,(dim,dim))) ,
                            Eigen(zeros(Float64, dim),zeros(ComplexF64, (dim,dim)))
                            )
    end

end

struct HamiltonianView{Float_Type <: Real}
    h::Array{      SubArray{Complex{Float_Type},0,Array{Complex{Float_Type},1},Tuple{Int64},true},2}
    v::Array{Array{SubArray{Complex{Float_Type},0,Array{Complex{Float_Type},1},Tuple{Int64},true},2},1}
    a::Array{Array{SubArray{Complex{Float_Type},0,Array{Complex{Float_Type},2},Tuple{Int64,Int64},true},2},2}
end

@inline function HamiltonianView(h_ops::H)::H where H <: HamiltonianOperators

    h = Array{SubArray{Complex{T},0,Array{Complex{T},1},Tuple{Int64},true},2}(undef,size(h_ops.h))

    v = Array{Array{SubArray{Complex{T},0,Array{Complex{T},1},Tuple{Int64},true},2},1}(undef,2)
    for i ∈ 1:length(v)
        v[i]   = Array{SubArray{Complex{T},0,Array{Complex{T},1},Tuple{Int64},true},2}(undef,size(h_ops.h))
    end

    a = Array{Array{SubArray{Complex{T},0,Array{Complex{T},2},Tuple{Int64,Int64},true},2},2}(undef,(2,2))
    for i ∈ 1:length(a)
        a[i] = Array{SubArray{Complex{T},0,Array{Complex{T},2},Tuple{Int64,Int64},true},2}(undef,size(h_ops.h))
    end

    return HamiltonianView(h,v,a)
end


# Random Hamiltonian for testing

function RandomHamiltonianOperators(dim, d=2)
    HamiltonianOperators(SharedArray(rand(ComplexF64,(dim,dim))),
                        [SharedArray(rand(ComplexF64,(dim,dim)))   for _=1:d ],
                        [SharedArray(rand(ComplexF64,(dim,dim)))   for _=1:d, _=1:d],
                         SharedArray(zeros(ComplexF64,(dim,dim))),
                        (zeros(ComplexF64, dim), zeros(Float64, dim), zeros(Int64, dim), zeros(Int64,(dim,dim))) ,
                        Eigen(zeros(Float64, dim),zeros(ComplexF64, (dim,dim)))
                        )
end


struct RandomHamiltonian{H <: HamiltonianOperators, Int_Type <: Int, Float_Type <: Real} <: HamiltonianDensity{Int_Type,Float_Type}
    h_ops::H
end

function RandomHamiltonian(dim,d=2)
    H0=RandomHamiltonianOperators(dim,d)
    RandomHamiltonian{typeof(H0),typeof(1),typeof(1.0)}(RandomHamiltonianOperators(dim,d))
end

function (rh::RandomHamiltonian)(k::Vector{F}) where F <: Number
    dim = size(rh.h_ops.h,1)
    rh.h_ops.h .= SharedArray(rand(ComplexF64,(dim,dim)))
    for i∈ length(rh.h_ops.v)
        rh.h_ops.v[i] .= SharedArray(rand(ComplexF64,(dim,dim)))
    end
    for i∈ length(rh.h_ops.a)
        rh.h_ops.a[i] .= SharedArray(rand(ComplexF64,(dim,dim)))
    end
end
