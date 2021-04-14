abstract type HamiltonianDensity{Int_Type <: Integer, Float_Type <: Real} end


struct HamiltonianOperators{Float_Type <: Real}
    h::Array{Complex{Float_Type},2}
    v::Array{Array{Complex{Float_Type},2},1}
    a::Array{Array{Complex{Float_Type},2},2}
    pz::Array{Complex{Float_Type},2}
    l_pjs::Array{Array{Complex{Float_Type},2},1}
    sl_pjs::Array{Array{Complex{Float_Type},2},1}
    structure::Array{Float_Type,2}
    vds::Array{Array{Complex{Float_Type},1},1}
    pz_0::Array{Complex{Float_Type},2}
    l_pjs_0::Array{Array{Complex{Float_Type},2},1}
    sl_pjs_0::Array{Array{Complex{Float_Type},2},1}
end

#=

struct HamiltonianPoperators{Float_Type <: Real}
    h::Array{Complex{Float_Type},2}
    v::Array{Array{Complex{Float_Type},2},1}
    a::Array{Array{Complex{Float_Type},2},2}
end

struct ProjectionOperators{Float_Type <: Real, A_type <: AbstractArray}
    pz::A_type{Complex{Float_Type},2}
    l_pjs::Array{A_type{Complex{Float_Type},2},1}
    sl_pjs::Array{A_type{Complex{Float_Type},2},1}
    structure::Array{Float_Type,2}
    vds::Array{Array{Complex{Float_Type},1},1}
    pz_0::Array{Complex{Float_Type},2}
    l_pjs_0::Array{Array{Complex{Float_Type},2},1}
    sl_pjs_0::Array{Array{Complex{Float_Type},2},1}
end

=#

#Easy constructor for section data
function HamiltonianOperators(dim,projN,slv_prj_N)
    HamiltonianOperators(Array(zeros(ComplexF64,(dim,dim))),
                        [Array(zeros(ComplexF64,(dim,dim))),   zeros(ComplexF64,(dim,dim))],
                        [Array(zeros(ComplexF64,(dim,dim)))    for _=1:2, _=1:2],
                        Array(zeros(ComplexF64,(dim,dim))),
                        [Array(zeros(ComplexF64,(dim,dim)))    for _=1:projN],
                        [Array(zeros(ComplexF64,(dim,dim)))   for _=1:slv_prj_N],
                        Array(zeros(Float64,(dim,dim))),
                        [Array(zeros(ComplexF64,(dim)))        for _=1:dim],
                        Array(zeros(ComplexF64, (dim,dim))),
                        [Array(zeros(ComplexF64,(dim,dim)))        for _=1:projN],
                        [Array(zeros(ComplexF64,(dim,dim)))        for _=1:slv_prj_N]
                        )
end

#Easy constructor for section data
function HamiltonianOperators(dim::Int64,pz::Array{Complex{Float64},2},pjs::Vector{Array{Complex{Float64},2}},slv_prjs::Vector{Array{Complex{Float64},2}})
    HamiltonianOperators(Array(zeros(ComplexF64,(dim,dim))),
                        [Array(zeros(ComplexF64,(dim,dim))), Array(zeros(ComplexF64,(dim,dim)))],
                        [Array(zeros(ComplexF64,(dim,dim))) for _=1:2, _=1:2],
                        Array(pz),
                        Array.(pjs),
                        Array.(slv_prjs),
                        Array(zeros(Float64,(dim,dim))),
                        [Array(zeros(ComplexF64,(dim)))     for _=1:dim],
                        Array(pz),
                        Array.(pjs),
                        Array.(slv_prjs)
                        )
end


struct HamiltonianView{Float_Type <: Real}
    h::Array{      SubArray{Complex{Float_Type},0,Array{Complex{Float_Type},1},Tuple{Int64},true},2}
    v::Array{Array{SubArray{Complex{Float_Type},0,Array{Complex{Float_Type},1},Tuple{Int64},true},2},1}
    a::Array{Array{SubArray{Complex{Float_Type},0,Array{Complex{Float_Type},2},Tuple{Int64,Int64},true},2},2}
end


@inline function HamiltonianView(h_ops::HamiltonianOperators{T})::HamiltonianView{T} where T <: Real

    h = Array{SubArray{Complex{T},0,Array{Complex{T},1},Tuple{Int64},true},2}(undef,size(h_ops.h))

    v = Array{Array{SubArray{Complex{T},0,Array{Complex{T},1},Tuple{Int64},true},2},1}(undef,2)
    v[1]   = Array{SubArray{Complex{T},0,Array{Complex{T},1},Tuple{Int64},true},2}(undef,size(h_ops.h))
    v[2]   = Array{SubArray{Complex{T},0,Array{Complex{T},1},Tuple{Int64},true},2}(undef,size(h_ops.h))

    a = Array{Array{SubArray{Complex{T},0,Array{Complex{T},2},Tuple{Int64,Int64},true},2},2}(undef,(2,2))
    a[1,1] = Array{SubArray{Complex{T},0,Array{Complex{T},2},Tuple{Int64,Int64},true},2}(undef,size(h_ops.h))
    a[1,2] = Array{SubArray{Complex{T},0,Array{Complex{T},2},Tuple{Int64,Int64},true},2}(undef,size(h_ops.h))
    a[2,1] = Array{SubArray{Complex{T},0,Array{Complex{T},2},Tuple{Int64,Int64},true},2}(undef,size(h_ops.h))
    a[2,2] = Array{SubArray{Complex{T},0,Array{Complex{T},2},Tuple{Int64,Int64},true},2}(undef,size(h_ops.h))
    return HamiltonianView(h,v,a)
end

#Methods to change the frame of the operators to the eigenbasis
@inline function action_Ad_H(U::Array{Complex{Float64},2},M::Array{Complex{Float64},2},Y::Array{Complex{Float64},2})::Nothing
    mul!(Y,Hermitian(M),U,true,false)
    mul!(M,(U)',Y, true, false)
    nothing
end

@inline function action_Ad_H(U::Array{Complex{Float64},2},M::AbstractArray{Array{Complex{Float64},2}},Y::Array{Complex{Float64},2}, i1::Int64, i2::Int64)::Nothing
    @inbounds mul!(Y,Hermitian(M[i1,i2]),U,true,false);
    @inbounds mul!(M[i1,i2],(U)',Y, true, false)
    nothing
end

@inline function action_Ad_H(U::Array{Complex{Float64},2},M::AbstractArray{Array{Complex{Float64},2}},Y::Array{Complex{Float64},2}, i1::Int64)::Nothing
    @inbounds mul!(Y,Hermitian(M[i1]),U,true,false);
    @inbounds mul!(M[i1],(U)',Y, true, false)
    nothing
end


@inline function action_Ad_H(T::Array{Complex{Float64},2},U::Array{Complex{Float64},2},M::Array{Complex{Float64},2},Y::Array{Complex{Float64},2})::Nothing
    mul!(Y,Hermitian(M),U, true,false)
    mul!(T,(U)',Y, true, false)
    nothing
end

@inline function cartan_frame(F::Eigen{Complex{Float64},Float64,Array{Complex{Float64},2},Array{Float64,1}},rep::HamiltonianOperators{Float64},y::NTuple{7,Array{ComplexF64,2}})::Nothing
    #storing eigenvalue differences and eigenvalues in hamiltonian
    #action_Ad_H(F.vectors,rep.h,y)
    rep.h .= Complex(0.0)
    dim =size(rep.h,1)
    @inbounds @simd for i=1:dim
        rep.h[i+(i-1)*dim] = F.values[i]
    end

    action_Ad_H(rep.pz, F.vectors, rep.pz_0, y[1]);


    action_Ad_H(F.vectors,rep.v,y[2],1);
    action_Ad_H(F.vectors,rep.v,y[3],2);
    action_Ad_H(F.vectors,rep.a,y[4],1,1);
    action_Ad_H(F.vectors,rep.a,y[5],1,2);

    action_Ad_H(F.vectors,rep.a,y[6],2,1);

    action_Ad_H(F.vectors,rep.a,y[7],2,2);


    @inbounds for i=1:length(rep.l_pjs)
        action_Ad_H(rep.l_pjs[i],  F.vectors, rep.l_pjs_0[i],   y[i])
    end

    #Storing the structure factors
    rep.structure .= abs2.(F.vectors)
    nothing
end

@inline function cartan!(A::Array{Complex{Float64},2}, E::Eigen{Complex{Float64},Float64,Array{Complex{Float64},2},Array{Float64,1}}, aux::Tuple{Vector{Complex{Float64}},Vector{Float64},Vector{Int64},Array{Int64,2}})

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
               Ptr{Float64}, Ptr{ComplexF64}, Ref{Int64}, Ptr{Int64},
               Ptr{ComplexF64}, Ref{Int64}, Ptr{Float64}, Ref{Int64},
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

@inline function cartan_frame_transform(h_bundle::HamiltonianOperators{Float64}, aux_matrix::NTuple{7,Array{Complex{Float64},2}})
    cartan_frame(eigen!(Hermitian(h_bundle.h)), h_bundle, aux_matrix)
end

@inline function cartan_frame_transform(h_bundle::HamiltonianOperators{Float64}, aux_matrix::NTuple{7,Array{Complex{Float64},2}}, E::Eigen{Complex{Float64},Float64,Array{Complex{Float64},2},Array{Float64,1}}, eig_aux::Tuple{Vector{Complex{Float64}},Vector{Float64},Vector{Int64},Array{Int64,2}})
    #cartan_frame(eigen!(Hermitian(h_bundle.h)), h_bundle, aux_matrix)
    cartan!(h_bundle.h, E, eig_aux)
    cartan_frame(E, h_bundle, aux_matrix)
end
