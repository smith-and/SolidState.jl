# A functional layer which allows the specific data charts define maps for data
"""
    DataMap{ChartType <: DataChart, KType <: KinematicDensity, LType <: AbstractArray}
"""
struct DataMap{ChartType <: DataChart, OType <: OperatorDensity, KType <: KinematicDensity, LType <: AbstractArray}
    d::Int
    Λ::LType
    K::KType
    chart::ChartType
    ops::OType
end

"""
"""
function range_scope(rng::Vector{Tuple{Symbol,N,N,Int64}} where N <: Number)::(Array{NTuple{N1,N}, N1} where {N1, N <: Number})
    collect(Iterators.product(LinRange.(getindex.(rng,2),getindex.(rng,3),getindex.(rng,4) )...))
end

"""
"""
function range_scope(type::Symbol,dim::Int64)::(Array{NTuple{N1,N}, N1} where {N1, N <: Number})
    if type==:dims
        return collect(Iterators.product(1:dim))
    end
end

"""
    bandwidth_range(asd, hd, Nb)
"""
function bandwidth_range(asd, hd, Nb)
    pts = getindex.(Ref(getindex(asd|>ASDGeometry,"bz_hs")),["Q1","M1","K1","Γ"])
    emin = 0.0 ; emax = 0.0;
    for k∈pts
        hd([0.0,0.0,0.0])
        emin = min(hd.h_ops.h|>Hermitian|>eigmin,emin)
        emax = max(hd.h_ops.h|>Hermitian|>eigmax,emax)
    end
    range_scope([(:ω,1.1*emin,1.1*emax,Nb)])
end

"""
    dim_kron(dim::Int, inds::AbstractVector)
dim_kron(10,[(1,2),(1,1),(2,2)])
"""
function dim_kron(dim::Int, inds::AbstractVector)
    Iterators.product(inds,collect(1:dim))|>collect
end
"""
    dim_kron(dim::Int, inds::NTuple{N,T} where {N, T <: AbstractRange})

dim_kron(10,(1:2,1:2))
"""
function dim_kron(dim::Int, inds::NTuple{N,T} where {N, T <: AbstractRange})
    Iterators.product(Iterators.product(inds...),1:dim)|>collect
end

"""
    input_process(asd,hd,indices,prange,brange)
"""
function input_process(asd,hd,indices,prange,brange)
    priors = range_scope(prange)
    if typeof(brange)==Symbol
        base = range_scope(brange,dim_h)
    elseif typeof(brange)==Tuple{Symbol,Symbol,Int}
        base = bandwidth_range(asd,hd,brange[3])
    else
        base = range_scope(brange)
    end

    newindices = indices
    if indices[1]==:dimkron

        inds = indices[2]
        dim = size(hd.h_ops.h,1)

        newindices = dim_kron(dim,inds)
    end

    (newindices,priors,base)
end

"""
    domain_map(Λ0,asd)
"""
function domain_map(Λ0,asd)
    if Λ0 == :auto
        return (asd|>ASDGeometry)["dxtal"][1][1:3,1:2]
    elseif Λ0 == :box
        return norm((asd|>ASDGeometry)["dxtal"][1][1:3,1])*[1.0 0.0 ; 0.0 1.0 ; 0.0 0.0]
    else
        return Λ0
    end
end

"""
    DataMap(dtype::(Type{T} where T <: DataChart), asd::Dict{String,Any}, indices0, prange, brange; Λ0=:auto, style=:normal)
"""
function DataMap(dtype::(Type{T} where T <: DataChart), asd::Dict{String,Any}, indices0, prange, brange; Λ0=:auto, style=:normal)
    Λ = domain_map(Λ0,asd)
    hd = asd|>TightBindingDensity
    (indices,priors,base) = input_process(asd,hd,indices0,prange,brange)
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,priors);
    chart = dtype(TensorChart(indices,priors,base,Complex{Float64};style=style))

    DataMap(dim_h,Λ,K,chart)
end

"""
    DataMap(dtype::(Type{T} where T <: DataChart), asd::Dict{String,Any}, hd::HamiltonianDensity, indices0, prange, brange; Λ0=:auto, style=:normal)
"""
function DataMap(dtype::(Type{T} where T <: DataChart), asd::Dict{String,Any}, hd::HamiltonianDensity, indices0, prange, brange; Λ0=:auto, style=:normal)
    Λ = domain_map(Λ0,asd)
    (indices,priors,base) = input_process(asd,hd,indices0,prange,brange)
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,priors)
    chart = dtype(TensorChart(indices,priors,base,Complex{Float64};style=style))

    DataMap(dim_h,Λ,K,chart)
end

"""
    DataMap(model::Symbol, dtype::(Type{T} where T <: DataChart), indices0, prange, brange; cachedir, mn, Λ0=:auto, style=:normal)
"""
function DataMap(model::Symbol, dtype::(Type{T} where T <: DataChart), indices0, prange, brange; cachedir, mn, Λ0=:auto, style=:normal)
    asd = BSON.load("$cachedir/$model/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/$model/hd-$(mn[1])-$(mn[2]).bson")
    Λ = domain_map(Λ0,asd)
    (indices,priors,base) = input_process(asd,hd,indices0,prange,brange)
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,priors)
    chart = dtype(TensorChart(indices,priors,base,Complex{Float64};style=style))

    DataMap(dim_h,Λ,K,chart)
end

"""
    @generated function evaluate_map(dm::DataMap, k::AbstractVector, section::T where T <: DataChart)
"""
@generated function evaluate_map(dm::DataMap, k::AbstractVector, section::T where T <: DataChart)
    name = fieldname(section,1)
    chart = fieldtypes(fieldtype(section,1))[end-3]
    if chart <: MArray
        quote
            dm.hd(k)
            dm.K(dm.hd)
            section.$name.data .= 0.0
            $(Symbol(name,:(_operators)))(section.$name, dm.K.k_m, dm.hd.h_ops, dm.d)
            $(Symbol(name,:(_evaluation)))(section.$name, dm.K.k_m, dm.K.hd.h_ops, dm.d)
            SArray(section.$name.data)
        end
    else
        quote
            dm.hd(k)
            dm.K(dm.hd)
            section.$name.data .= 0.0
            $(Symbol(name,:(_operators)))(section.$name, dm.K.k_m, dm.K.hd.h_ops, dm.d)
            $(Symbol(name,:(_evaluation)))(section.$name, dm.K.k_m, dm.K.hd.h_ops, dm.d, dm.ops)
            copy(section.$name.data)
        end
    end
end

"""
    (dm::DataMap)(n::AbstractVector)

passes  dm.Λ*n to the DataMap
"""
function (dm::DataMap)(n::AbstractVector)
    evaluate_map(dm, dm.Λ*n, dm.chart)
end

"""
(dm::DataMap)(n::NTuple{N,F} where {N,F <: Number})

passes eachcol(dm.Λ).*n|>sum to the DataMap
"""
function (dm::DataMap)(n::NTuple{N,F} where {N,F <: Number})
    evaluate_map(dm, eachcol(dm.Λ).*n|>sum, dm.chart)
end
