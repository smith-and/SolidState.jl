#Methods to change the frame of the operators to the eigenbasis

#= Projection Operator Basis Changes
    action_Ad_H(rep.pz, F.vectors, rep.pz_0, y[1]);
@inbounds for i=1:length(rep.l_pjs)
    action_Ad_H(rep.l_pjs[i],  F.vectors, rep.l_pjs_0[i],   y[i])
end

#Storing the structure factors
rep.structure .= abs2.(F.vectors)
nothing
=#

#In-place version of eigen
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

@inline function cartan_transform(h_bundle::HamiltonianOperators)
    try
        cartan!(h_bundle.h, h_bundle.E, h_bundle.eig_aux)
        cartan_frame(h_bundle.E, h_bundle, h_bundle.aux_matrix)
    catch
<<<<<<< HEAD
        eigen(Hermitian(h_bundle.h)) |> identity() do E
=======
        let E = eigen!(Hermitian(h_bundle.h))
>>>>>>> development
            h_bundle.E.values .= E.values
            h_bundle.E.vectors .= E.vectors
        end

        cartan_frame(h_bundle.E, h_bundle, h_bundle.aux_matrix)
    end
end
