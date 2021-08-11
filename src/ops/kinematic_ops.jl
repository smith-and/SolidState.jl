########################################
#### Hamiltonian Operators
########################################

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

########################################
#### Cartan Frame Diagonalization
########################################

# In-place version of eigen
@inline function cartan!(A::AbstractArray{Complex{F},2}, E::Eigen{Complex{F},F,Array{Complex{F},2},Array{F,1}}, aux::Tuple{Vector{Complex{F}},Vector{F},Vector{Int64},Array{Int64,2}}) where F <: AbstractFloat

    LinearAlgebra.chkstride1(A)
    n = Int64(LinearAlgebra.checksquare(A))
    lda = Int64(max(1,LinearAlgebra.stride(A,2)))
    m = Ref{Int64}()

    w = E.values
    Z = E.vectors

    isuppz = aux[4]
    work   = aux[1]
    lwork  = Int64(-1)
    rwork  = aux[2]
    lrwork = Int64(-1)
    iwork  = aux[3]
    liwork = Int64(-1)
    info   = Ref{Int64}()

    for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
        ccall((:zheevr_, LinearAlgebra.BLAS.libblas), Cvoid,
              (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{Int64},
               Ptr{ComplexF64}, Ref{Int64}, Ref{ComplexF64}, Ref{ComplexF64},
               Ref{Int64}, Ref{Int64}, Ref{ComplexF64}, Ptr{Int64},
               Ptr{AbstractFloat}, Ptr{ComplexF64}, Ref{Int64}, Ptr{Int64},
               Ptr{ComplexF64}, Ref{Int64}, Ptr{AbstractFloat}, Ref{Int64},
               Ptr{Int64}, Ref{Int64}, Ptr{Int64}),
              'V', 'A', 'U', n,
              A, lda, 0.0, 0.0,
              1, n, -1.0, m,
              w, Z, n, isuppz,
              work, lwork, rwork, lrwork,
              iwork, liwork, info)
        #LinearAlgebra.LAPACK.chklapackerror(Int32(info[]))
        if i == 1
            lwork = Int64(real(work[1]))
            resize!(work, lwork)
            lrwork = Int64(rwork[1])
            resize!(rwork, lrwork)
            liwork = iwork[1]
            resize!(iwork, liwork)
        end
    end
    E.values  .= w
    E.vectors .= Z
    LinearAlgebra.sorteig!(E.values, E.vectors)

    nothing
end

@inline function action_Ad_H(U::AbstractArray{Complex{F},2},M::A where A <: Union{AbstractArray{SharedArray{Complex{F},2}},AbstractArray{Array{Complex{F},2}}},Y::AbstractArray{Complex{F},2}, i1::Int64)::Nothing where F <: AbstractFloat
    @fastmath @inbounds mul!(Y,Hermitian(M[i1]),U,true,false);
    @fastmath @inbounds mul!(M[i1],(U)',Y, true, false)
    nothing
end

@inline function cartan_frame(F::Eigen{Complex{F0},F0,Array{Complex{F0},2},Array{F0,1}},rep::H0, y::AbstractArray{Complex{F0},2})::Nothing where {F0 <: AbstractFloat, H0 <: HamiltonianOperators}
    #storing eigenvalue differences and eigenvalues in hamiltonian
    #action_Ad_H(F.vectors,rep.h,y)
    rep.h .= Complex(0.0)
    dim =size(rep.h,1)
    for i=1:dim
        @fastmath @inbounds rep.h[i+(i-1)*dim] = F.values[i]
        nothing
    end

    for i = 1:length(rep.v)
        action_Ad_H(F.vectors,rep.v,y,i);
        nothing
    end

    for i = 1:length(rep.a)
        action_Ad_H(F.vectors,rep.a,y,i);
        nothing
    end

    nothing

end
#
# @inline function cartan_transform(h_bundle::HamiltonianOperators)
#     try
#         cartan!(h_bundle.h, h_bundle.E, h_bundle.eig_aux)
#         cartan_frame(h_bundle.E, h_bundle, h_bundle.aux_matrix)
#     catch
#         let E = eigen!(Hermitian(h_bundle.h))
#             h_bundle.E.values .= E.values
#             h_bundle.E.vectors .= E.vectors
#         end
#
#         cartan_frame(h_bundle.E, h_bundle, h_bundle.aux_matrix)
#     end
# end

@inline function cartan_transform(h_bundle::HamiltonianOperators)
    let E = eigen!(Hermitian(h_bundle.h))
        h_bundle.E.values .= E.values
        h_bundle.E.vectors .= E.vectors
    end

    cartan_frame(h_bundle.E, h_bundle, h_bundle.aux_matrix)
end

########################################
#### Kinematic Operators
########################################

#data for set of "kinematic" obsevarble operators which are used to compute response functions
struct KinematicOperators{AR <: AbstractArray, AC <: AbstractArray}
    df::Vector{AR}
    dω::AR
    Δ::Array{AC,1}
    re::Array{AC,1}
    rire::Array{AC,2}
    pz::AC
    pz_0::AC
end

#Constructor for the Kinematic Operator Bundle
function KinematicOperators(dim::Integer, priorsN::Integer,d=2; style::Symbol=:normal)
    if style == :shared
        return KinematicOperators(
            [SharedArray(zeros(Float64,(dim,dim))) for i=1:priorsN],
            SharedArray(zeros(Float64,(dim,dim))),
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SharedArray(zeros(ComplexF64,(dim,dim))) for _=1:d, _=1:d],
            SharedArray(zeros(ComplexF64,(dim,dim))),
            SharedArray(zeros(ComplexF64,(dim,dim)))
        )
    elseif style == :static
        return KinematicOperators(
            [SMatrix{dim,dim}(zeros(Float64,(dim,dim))) for i=1:priorsN],
            SMatrix{dim,dim}(zeros(Float64,(dim,dim))),
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))) for _=1:d, _=1:d],
            SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim))),
            SMatrix{dim,dim}(zeros(ComplexF64,(dim,dim)))
        )
    else
        return KinematicOperators(
            [(zeros(Float64,(dim,dim))) for i=1:priorsN],
            (zeros(Float64,(dim,dim))),
            [(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [(zeros(ComplexF64,(dim,dim))) for _=1:d],
            [(zeros(ComplexF64,(dim,dim))) for _=1:d, _=1:d],
            (zeros(ComplexF64,(dim,dim))),
            (zeros(ComplexF64,(dim,dim)))
        )
    end
end


struct KinematicDensity{H_Type <: HamiltonianDensity, K_Type <: KinematicOperators, AUXA <: AbstractArray, PA <: AbstractArray}
    hd::H_Type
    k_m::K_Type
    aux_real::AUXA
    priors::PA
end

# Constructor with density
function KinematicDensity(hd::H where H <: HamiltonianDensity{Int64,F}, priors::Vector{Tuple{Symbol,F,F,Int64}}; style::Symbol=:normal) where F <: Number
    inputs = ( hd, KinematicOperators(size(hd.h_ops.h,1),prod(getindex.(priors,4)); style=style), zeros(Complex{F}, size(hd.h_ops.h)),(x->[x...]).(collect(Iterators.product([range(p[2],p[3],length=p[4]) for p∈priors]...))[:]) )
    KinematicDensity{typeof.(inputs)...}(inputs...)
end

function KinematicDensity(hd::H where H <: HamiltonianDensity{Int64,F}, priors::Vector{LinRange{F}}; style::Symbol=:normal) where F <: AbstractFloat
    inputs = (hd, KinematicOperators(size(hd.h_ops.h,1),prod(length.(priors)); style=style), zeros(Complex{F}, size(hd.h_ops.h)),(x->[x...]).(collect(Iterators.product(priors...))[:]))
    KinematicDensity{typeof.(inputs)...}(inputs...)
end

function KinematicDensity(hd::H where H <: HamiltonianDensity{Int64,F}, priors::Vector{Vector{F}}; style::Symbol=:normal) where F <: AbstractFloat
    inputs = (hd, KinematicOperators(size(hd.h_ops.h,1),length(priors); style=style), zeros(eltype(hd.h_ops.h), size(hd.h_ops.h)),priors)
    KinematicDensity{typeof.(inputs)...}(inputs...)
end

function KinematicDensity(hd::HD where HD <: HamiltonianDensity, priors::Array{NTuple{N,Float64},N} where N; style::Symbol=:normal)
    KinematicDensity(hd, KinematicOperators(size(hd.h_ops.h,1),length(priors); style=style), zeros(Complex{typeof(0.0)},size(hd.h_ops.h)),((x->[x...]).(priors))[:])
end

########################################
#### [ri,re] matrix elements
########################################

@inline function rire_add(n::Int64, m::Int64, p::Int64, new_val::Complex{F}, va::AbstractArray{Complex{F},2}, vb::AbstractArray{Complex{F},2}, ra::AbstractArray{Complex{F},2}, rb::AbstractArray{Complex{F},2} ) where F <: AbstractFloat
    if !(p==n||p==m)
        @inbounds @fastmath    new_val += (va[n,p]*rb[p,m]-rb[n,p]*va[p,m])
    end
end

@inline function km_rire(rireab::AbstractArray{Complex{F},2}, wab::AbstractArray{Complex{F},2}, va::AbstractArray{Complex{F},2}, vb::AbstractArray{Complex{F},2}, ra::AbstractArray{Complex{F},2}, rb::AbstractArray{Complex{F},2}, Δb::AbstractArray{Complex{F},2}, Δa::AbstractArray{Complex{F},2}, inv_dω::AbstractArray{Complex{F},2})::Nothing where F <: AbstractFloat
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

@inline function rire_ops(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, inv_dω::AbstractArray{Complex{F},2}) where F <: AbstractFloat
    for a ∈ 1:length(cfr.v)
        for b ∈ 1:length(cfr.v)
            @inbounds km_rire(km.rire[a,b], cfr.a[a,b], cfr.v[a], cfr.v[b], km.re[a], km.re[b], km.Δ[a], km.Δ[b], inv_dω)
        end
    end

    nothing
end

@inline function kinematic_ops_fast(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, ϵs::Vector{F}, r_aux::AbstractArray{Complex{F},2}) where F <: AbstractFloat
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
                km.re[1][idx] = Complex(0.0,-(r_aux[idx]).re) * (cfr.v[1][idx]) ;
                km.re[2][idx] = Complex(0.0,-(r_aux[idx]).re) * (cfr.v[2][idx]) ;
                @inbounds @fastmath for p=eachindex(priors)
                    km.df[p][idx] = km.df[p][j,j] - km.df[p][i,i]
                end
            end
        end
    end

    # Assigns the rest of the operators in the section
    # rire_ops(km, cfr, priors, r_aux)

    nothing
end

@inline function kinematic_ops(km::KinematicOperators, cfr::HamiltonianOperators, priors::Vector{Vector{F}}, ϵs::Vector{F}, r_aux::AbstractArray{Complex{F},2}) where F <: AbstractFloat
    kinematic_ops_fast(km,cfr,priors,ϵs,r_aux)
    # Assigns the rest of the operators in the section
    rire_ops(km, cfr, priors, r_aux)
    nothing
end

function (k_obs::KinematicDensity)(k::Vector{F}) where F <: Number #::Tuple{KinematicOperators{AbstractFloat},HamiltonianOperators}
    k_obs.hd(k)
    cartan_transform(k_obs.hd.h_ops)
    kinematic_ops(k_obs.k_m, k_obs.hd.h_ops, k_obs.priors, k_obs.hd.h_ops.E.values, k_obs.aux_real)

    nothing
end

function evaluate_km(k_obs::KinematicDensity, k::Vector{F}) where F <: Number #::Tuple{KinematicOperators{AbstractFloat},HamiltonianOperators}
    k_obs.hd(k)
    cartan_transform(k_obs.hd.h_ops)
    kinematic_ops_fast(k_obs.k_m, k_obs.hd.h_ops, k_obs.priors, k_obs.hd.h_ops.E.values, k_obs.aux_real)

    nothing
end

########################################
#### Projection Operators
########################################
abstract type ProjectionOperator end

function init_sqoi(asd,names)
    asdb = asd|>SolidState.ASDBasics
    sqoi = Dict{Symbol, Any}(:asdb=>asdb,:names=>names)
    map(names) do name
        sqoi[name] = getfield.(asdb["sqo"],name)
        sqoi[Symbol(:unq,name)] = sqoi[name]|>unique
    end
    sqoi
end

########################################
#### Atom and Layer Projector
########################################

struct ProjectorSL{AC <: AbstractArray, BV <: AbstractVector, DA <: AbstractDict} <: ProjectionOperator
    layer::AC
    atom::AC
    layer0::BV
    atom0::BV
    sqoi::DA
end

function ProjectorSL(dim, Nlayer, Natom, sqoi=Dict{Symbol,Any}())
    ProjectorSL(
        zeros(typeof(0.0),(dim,Nlayer)),
        zeros(typeof(0.0),(dim,Natom)),
        [(BitArray(undef,dim) .= 0 )   for _=1:Natom],
        [(BitArray(undef,dim) .= 0 )   for _=1:Nlayer],
        sqoi
    )
end

function init_proj!(p::ProjectorSL,name)
    map(enumerate(p.sqoi[Symbol(:unq,name)])) do (i,a)
        inds = (a.==p.sqoi[name])
        getfield(p,Symbol(name,0))[i] .= inds
        # proj = getfield(p,name)[i]
        # proj[diagind(proj)[inds]] .= Complex(1.0)
        # proj
    end
end

function ProjectorSL(asd::Dict{String,Any})
    sqoi = init_sqoi(asd,[:layer,:atom])
    p = ProjectorSL(length(sqoi[:asdb]["sqo"]),length(sqoi[Symbol(:unq,:layer)]),length(sqoi[Symbol(:unq,:atom)]),sqoi)

    init_proj!.(Ref(p),[:atom,:layer])
    p
end

function eigenvec_weight_sl(evec,inds)
    sum(abs2.(evec[inds]))
end

function calc_weights(hd,p::ProjectorSL,name)

    dim = size(getfield(p,name),1)
    Np = size(getfield(p,name),2)

    for i in 1:dim
        for j in 1:Np
            getfield(p,name)[i + (j-1)*dim] = eigenvec_weight_sl(hd.h_ops.E.vectors[:,i], getfield(p,Symbol(name,0))[j])
        end
    end
    getfield(p,name)
end

########################################
#### Layer Parity Projector
########################################

struct ProjectorLP{AC <: AbstractArray, BV <: AbstractArray, DA <: AbstractDict} <: ProjectionOperator
    layer::AC
    layer0::BV
    sqoi::DA
end

function ProjectorLP(dim, Nlayer, sqoi=Dict{Symbol,Any}())
    ProjectorLP(
        zeros(typeof(0.0),(dim,Nlayer)),
        [zeros(typeof(Complex(0.0)),(dim,dim)) for _=1:Nlayer],
        sqoi
    )
end

function init_proj!(p::ProjectorLP,name)
    map(enumerate(p.sqoi[Symbol(:unq,name)])) do (i,a)
        inds = (a.==p.sqoi[name])
        getfield(p,Symbol(name,0))[i] .= inds
        # proj = getfield(p,name)[i]
        # proj[diagind(proj)[inds]] .= Complex(1.0)
        # proj
    end
end

function ProjectorLP(asd::Dict{String,Any})
    names = [:layer]
    sqoi = init_sqoi(asd,names)

    p = ProjectorLP(length(sqoi[:asdb]["sqo"]),length.(getindex.(Ref(sqoi),Symbol.(:unq,names))...),sqoi)

    init_proj!.(Ref(p),names)

    p
end

function eigenvec_weight_lp(evec,op)
    evec'*op*evec
end

function calc_weights(hd,p::ProjectorLP,name)

    dim = size(getfield(p,name),1)
    Np = size(getfield(p,name),2)

    for i in 1:dim
        for j in 1:Np
            getfield(p,name)[i + (j-1)*dim] = eigenvec_weight_lp(hd.h_ops.E.vectors[:,i], getfield(p,Symbol(name,0))[j])
        end
    end
    getfield(p,name)
end
