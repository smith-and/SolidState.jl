"""
   abstract type ChartType end
"""
abstract type ChartType end

"""
const RangeType = Union{
   Tuple{Symbol,NTuple{N,T}} where {N, T <: AbstractRange}, Tuple{Symbol,AV} where AV <: AbstractVector,
   Vector{Tuple{Symbol,N,N,I}},
   Tuple{Symbol,Symbol,I},
   Vector{NTuple{NT,I}}} where {N <: Number, I <: Integer, NT
}

"""
const RangeType = Union{
         Tuple{Symbol,NTuple{N,T}} where {N, T <: AbstractRange}, Tuple{Symbol,AV} where AV <: AbstractVector,
         Vector{Tuple{Symbol,N,N,I}},
         Tuple{Symbol,Symbol,I},
         Vector{NTuple{NT,I}}} where {N <: Number, I <: Integer, NT
      }

"""
struct ChartInfo{T <: ChartType, RT1 <: RangeType, RT2 <: RangeType, RT3 <: RangeType, CAT <: Union{Tuple{Int,Int}}}
    type::Type{T}
    indices::RT1
    priors::RT2
    base::RT3
    mn::CAT
end
"""
struct ChartInfo{T <: ChartType, RT1 <: RangeType, RT2 <: RangeType, RT3 <: RangeType, CAT <: Union{Tuple{Int,Int}}}
    type::Type{T}
    indices::RT1
    priors::RT2
    base::RT3
    mn::CAT
end

"""
function ChartInfo()
   ChartInfo(ChartType,[(1,)],[(1,)],[(1,)],(0,0))
end
"""
function ChartInfo()
   ChartInfo(ChartType,[(1,)],[(1,)],[(1,)],(0,0))
end

"""
   DataChart{IndexT <: AbstractArray, PriorT <: AbstractArray, BaseT <: AbstractArray, DataT <: AbstractArray}

A generic structure to hold Tensor Data
"""
struct DataChart{CT <: ChartType, IndexT <: AbstractArray, PriorT <: AbstractArray, BaseT <: AbstractArray, DataT <: AbstractArray}
    name::Type{CT}
    indices::IndexT
    priors::PriorT
    base::BaseT
    data::DataT
    l_i::Int
    l_p::Int
    l_b::Int
end

"""
   DataChart(Indices::AbstractArray, Priors::AbstractArray, Base::AbstractArray, DType::DataType; style=:normal)
"""
function DataChart(chart_type::(Type{T} where T <: ChartType), Indices::AbstractArray, Priors::AbstractArray, Base::AbstractArray, DType::DataType; style=:normal)
   if      style == :static
      DataChart(chart_type,Indices, Priors, Base, @MArray(zeros(DType,(size(Indices)...,size(Priors)...,size(Base)...)...)) ,length(Indices),length(Priors),length(Base))
   elseif  style ==:shared
      DataChart(chart_type,SharedArray(Indices), SharedArray(Priors), SharedArray(Base), SharedArray(zeros(DType,(size(Indices)...,size(Priors)...,size(Base)...))),length(Indices),length(Priors),length(Base))
   else
      DataChart(chart_type,Indices, Priors, Base, (zeros(DType,(size(Indices)...,size(Priors)...,size(Base)...)...)),length(Indices),length(Priors),length(Base))
   end
end

#############################
#### DataMap
#############################

struct DataMap{ChartType <: DataChart, KType <: KinematicDensity, PT<: Projector, LType <: AbstractArray, IType <: NamedTuple}
    d::Int
    Ω::Float64
    Λ::LType
    K::KType
    chart::ChartType
    P::PT
    input::IType
end

function range_scope(rng::Vector{Tuple{Symbol,N,N,Int64}} where N <: Number)::(Array{NTuple{N1,N}, N1} where {N1, N <: Number})
    collect(Iterators.product(LinRange.(getindex.(rng,2),getindex.(rng,3),getindex.(rng,4) )...))
end

function range_scope(type::Symbol,dim::Int64)::(Array{NTuple{N1,N}, N1} where {N1, N <: Number})
    if type==:dims
        return collect(Iterators.product(1:dim))
    end
end

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

function dim_kron(dim::Int, inds::AbstractVector)
    Iterators.product(inds,collect(1:dim))|>collect
end

function dim_kron(dim::Int, inds::NTuple{N,T} where {N, T <: AbstractRange})
    Iterators.product(Iterators.product(inds...),1:dim)|>collect
end

function input_process(asd,hd,indices,prange,brange)
    priors = range_scope(prange)
    if typeof(brange)==Symbol
        base = range_scope(brange,dim_h)
    elseif typeof(brange)==Tuple{Symbol,Symbol,Int}
        base = bandwidth_range(asd,hd,brange[3])
    elseif indices[1]==:bands||indices[1]==:bandstructure
        base = SolidState.Main.path_points(asd,brange...)
    elseif brange[1] == :dim
        base = collect(1:size(hd.h_ops.h,1))
    elseif brange[1] == :one
        base = collect(1:1)
    else
        base = range_scope(brange)
    end

    newindices = indices
    if indices[1]==:dimkron
        inds = indices[2]
        dim = size(hd.h_ops.h,1)

        newindices = dim_kron(dim,inds)
    elseif indices[1]==:bandstructure
        dim = size(hd.h_ops.h,1)
        newindices = ones(dim,dim+1)
    elseif indices[1]==:bands
        dim = size(hd.h_ops.h,1)
        newindices = ones(dim)
    end

    (newindices,priors,base)
end

function domain_map(Λ0,asd)
    if Λ0 == :auto
        return (asd|>ASDGeometry)["dxtal"][1][1:3,1:2]
    elseif Λ0 == :box
        return norm((asd|>ASDGeometry)["dxtal"][1][1:3,1])*[1.0 0.0 ; 0.0 1.0 ; 0.0 0.0]
    else
        return Λ0
    end
end

function DataMap( asd0::Function, dtype::(Type{T} where T <: ChartType), indices0, prange, brange; pnames = Symbol[], Λ0=:auto, style=:normal)
    asd = asd0()
    Λ = domain_map(Λ0,asd)
    Ω = asd["blv"]|>det

    hd = asd|>TightBindingDensity
    (indices,priors,base) = input_process(asd,hd,indices0,prange,brange)
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,priors);
    chart = DataChart(dtype,indices,priors,base,Complex{Float64};style=style)
    P = Projector(asd,pnames)

    inputs = (
        asd = asd0,
        mn = (1,1),
        indices= indices0,
        priors = prange,
        base   = brange,
        Λ      = Λ0,
        projectors = pnames,
        style = style,
    )

    println(inputs)

    DataMap(dim_h,Ω,Λ,K,chart,P,inputs)
end

function DataMap( asd0::Function, hd::H where H <: HamiltonianDensity, dtype::(Type{T} where T <: ChartType), indices0, prange, brange; pnames = [:layer],  Λ0=:auto, style=:normal)
    asd = asd0()
    Λ = domain_map(Λ0,asd)
    Ω = asd["blv"]|>det

    (indices,priors,base) = input_process(asd,hd,indices0,prange,brange)
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,priors)
    chart = DataChart(dtype,indices,priors,base,Complex{Float64};style=style)
    P = Projector(asd,pnames)

    inputs = (
        asd = asd0,
        mn = (1,1),
        indices= indices0,
        priors = prange,
        base   = brange,
        Λ      = Λ0,
        projectors = pnames,
        style = style,
    )

    DataMap(dim_h,Ω,Λ,K,chart,P,inputs)
end

"""
    DataMap(model::Symbol, dtype::(Type{T} where T <: ChartType), indices0, prange, brange; cachedir, mn, Λ0=:auto, style=:normal)

Uses cached model
"""
function DataMap(asd, mn::Tuple, dtype::(Type{T} where T <: ChartType), indices0, prange, brange; cachedir=ENV["cachedir"],  pnames = [:layer], Λ0=:auto, style=:normal)
    asd0 = BSON.load("$cachedir/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")
    Λ = domain_map(Λ0,asd0)
    Ω = asd0["blv"]|>det

    (indices,priors,base) = input_process(asd0,hd,indices0,prange,brange)
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,priors)
    chart = DataChart(dtype,indices,priors,base,Complex{Float64};style=style)
    P = Projector(asd0,pnames)

    inputs = (
        asd    = asd,
        mn     = mn,
        indices= indices0,
        priors = prange,
        base   = brange,
        Λ      = Λ0,
        projectors = pnames,
        style = style,
    )

    DataMap(dim_h,Ω,Λ,K,chart,P,inputs)
end

"""
    DataMap(; asd, dtype::(Type{T} where T <: ChartType), indices0, prange, brange, cachedir, mn, Λ0=:auto, style=:normal)

Uses cached model
"""
function DataMap(; asd, mn, dtype::(Type{T} where T <: ChartType), indices,priors,base, cachedir=ENV["cachedir"],  pnames = [:layer], Λ0=:auto, style=:normal, kargs...)
    asd0 = BSON.load("$cachedir/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")
    Λ = domain_map(Λ0,asd0)
    Ω = asd0["blv"]|>det
    (indices0,prange,brange) = input_process(asd0,hd,indices,priors,base)
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,prange)
    chart = DataChart(dtype,indices0,prange,brange, Complex{Float64};style=style)
    P = Projector(asd0,pnames)

    inputs = (
        asd    = asd,
        mn     = mn,
        indices= indices,
        priors = priors,
        base   = base,
        Λ      = Λ0,
        projectors = pnames,
        style = style,
    )

    DataMap(dim_h,Ω,Λ,K,chart,P,inputs)
end

"""
    DataMap(pmodel,pchart; cachedir=ENV["cachedir"], Λ0=:auto, style=:normal, kargs...)

"""
function DataMap(pmodel,pchart; cachedir=ENV["cachedir"],  pnames = [:layer], Λ0=:auto, style=:normal, kargs...)
    asd0 = SolidState.make_model(pmodel)
    hd   = SolidState.make_tb_model(pmodel)
    Λ = domain_map(Λ0,asd0)
    Ω = asd0["blv"]|>det
    (indices0,prange,brange) = input_process(asd0,hd,pchart[:indices],pchart[:priors],pchart[:base])
    dim_h   = size(hd.h_ops.h,1)
    K       = KinematicDensity(hd,prange)
    chart = DataChart(pchart[:dtype],indices0,prange,brange,Complex{Float64};style=style)
    P = Projector(asd0,pnames)

    inputs = (
        asd    = asd,
        mn     = mn,
        indices= indices0,
        priors = prange,
        base   = brange,
        Λ      = Λ0,
        projectors = pnames,
        style = style,
    )

    DataMap(dim_h,Ω,Λ,K,chart,P,inputs)
end

@generated function evaluate_map(dm::DataMap, k::AbstractVector)
    dtype = ((dm|>fieldtypes)[5]|>fieldtypes)[1]
    dname = "$dtype"[6:end-1]
    quote
        $(Symbol(dname,:(_evaluation)))(dm.chart, dm.P, dm.K, k, dm.d)
    end
end

function (dm::DataMap)(n::AbstractVector)
    evaluate_map(dm, dm.Λ*n)
end

function (dm::DataMap)()
    evaluate_map(dm, [0.0,0.0,0.0])
end


"""
    VOL
"""
abstract type VOL <: ChartType end

function VOL_evaluation(tc::DataChart, P::Projector, K0::KinematicDensity, k::AbstractVector , dim_ℋ::Int64)
    evaluate_km(K0,k)
    K = K0.k_m
    H = K0.hd
    tc.data .= 0.0
    for (ii,(a,b,c)) ∈ enumerate(tc.indices)
        for (ib,(ω,)) ∈ enumerate(tc.base)
            for (ip,(T,μ,δ)) ∈ enumerate(tc.priors)
                @fastmath idx = ii + ((ip-1) + (ib-1)*tc.l_p)*tc.l_i
                @fastmath @inbounds tc.data[idx] = Complex(1.0)
            end
        end
    end
    SArray(tc.data)
end

##################################
#### DataIntegral
##################################

"""
    DataIntegral{DM <: DataMap, DA <: AbstractArray, AV <: AbstractVector, IA <: AbstractArray, EA <: AbstractArray, CIT <: ChartInfo}
"""
struct DataIntegral{DM <: DataMap, DA <: AbstractArray, AV <: AbstractVector, IA <: AbstractArray, EA <: AbstractArray, CIT <: ChartInfo}
    dm::DM
    a::AV
    b::AV
    times::EA
    err::EA
    evals::IA
    bytes::IA
    data::DA
    chartinfo::CIT
    cachedir::String
end

#Construction from map
"""
    DataIntegral(dm::DM, a::AV=[0.0,0.0], b::AV=[1.0,1.0];cachedir::String="none", ranges::ChartInfo=ChartInfo())::DataIntegral where {AV <: AbstractVector, DM <: DataMap, PT <: Union{AbstractWorkerPool,Symbol}}
"""
function DataIntegral(dm::DM, a::AV=[0.0,0.0], b::AV=[1.0,1.0];cachedir=ENV["cachedir"], ranges::ChartInfo=ChartInfo())::DataIntegral where {AV <: AbstractVector, DM <: DataMap, PT <: Union{AbstractWorkerPool,Symbol}}
    DataIntegral(dm,a,b,typeof(1.0)[],typeof(1.0)[],typeof(1)[],typeof(1)[],typeof(dm.chart.data)[],ranges,cachedir)
end

#Construction with map description
"""
    DataIntegral(datatype::Type{T} where T <: ChartType,asd::Dict{String,Any},hd::HamiltonianDensity,indices,priors,base; a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector, DM <: DataMap}
"""
function DataIntegral(asd::Dict{String,Any},hd::HamiltonianDensity,datatype::ChartType,indices,priors,base; a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector, DM <: DataMap}
    dm  = DataMap(asd,hd,datatype,indices,priors,base);
    DataIntegral(dm,a,b, ranges=ChartInfo(datatype,indices,priors,base,(0,0)),cachedir=pwd())
end

"""
    DataIntegral(asd, mn, datatype::Type{T} where T <: ChartType,indices,priors, base, cachedir=ENV["cachedir"], a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector}

Construction with cached objects & description
"""
function DataIntegral(asd, mn, datatype,indices,priors, base, cachedir=ENV["cachedir"], a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector}
    ASD = BSON.load("$cachedir/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")
    dm = DataMap(()->ASD,hd,datatype,indices,priors,base);
    DataIntegral(dm,a,b; ranges=ChartInfo(datatype,indices,priors,base,mn),cachedir="$cachedir/$asd")
end

"""
    DataIntegral(; asd, dtype::Type{T} where T <: ChartType,indices,priors, base, cachedir, mn, a::AV=[0.0,0.0], b::AV =[1.0,1.0], kargs...)::DataIntegral  where {AV <: AbstractVector}

Construction with cached objects & description
"""
function DataIntegral(; asd, mn, dtype::ChartType, indices, priors, base, a::AV=[0.0,0.0], b::AV =[1.0,1.0], cachedir=ENV["cachedir"], kargs...)::DataIntegral  where {AV <: AbstractVector}
    asd0 = BSON.load("$cachedir/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/$asd/hd-$(mn[1])-$(mn[2]).bson")
    dm = DataMap(asd0,hd,dtype,indices,priors,base);
    DataIntegral(dm,a,b; ranges=ChartInfo(dtype,indices,priors,base,mn),cachedir="$cachedir/$asd")
end

"""
    DataIntegral(ci::ChartInfo, cachedir::String; a::AbstractVector=[0.0,0.0], b::AbstractVector =[1.0,1.0])::DataIntegral
"""
function DataIntegral(ci::ChartInfo, cachedir::String=ENV["cachedir"]; a::AbstractVector=[0.0,0.0], b::AbstractVector =[1.0,1.0])::DataIntegral
    asd = BSON.load("$cachedir/asd-$(ci.mn[1])-$(ci.mn[2]).bson")
    hd  = data_import("$cachedir/hd-$(ci.mn[1])-$(ci.mn[2]).bson")
    dm = DataMap(asd,hd,datatype,indices,priors,base);
    DataIntegral(dm,a,b; ranges=ci,cachedir=cachedir)
end

"""
    integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::Symbol = :none; rtol=1e-12, atol=1e-12, evals::Int=typemax(Int), kwargs...)::Float64

Core Call to HCubature
"""
function integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::Symbol; rtol=1e-12, atol=1e-12, evals::Int=typemax(Int), kwargs...)::Float64
    dm.chart.data .*= 0.0
    integral,int_error = hcubature(dm, a, b; norm=norm, rtol=rtol, atol=atol, maxevals=evals, initdiv=1)::Tuple{typeof(dm.chart.data),typeof(1.0)}
    dm.chart.data .= integral./(dm.Ω/(2π)^size(dm.Λ,1))

    int_error/(dm.chart.data|>length)
end

"""
function strip_divide(a,b,nw)

 For parallel Integration, to divide hyper k-cube into (k-1) slabs. cube is defined between to points a & b. nw specifies the number of workers, i.e. strips
"""
function strip_divide(a,b,nw)
    as = [[a[1],b[2]*(i-1)/nw] for i ∈ 1:nw]
    bs = [[b[1],b[2]*i/nw] for i∈1:nw]
    return as,bs
end


"""
    cointegrate(ci::ChartInfo,evals,cachedir,a,b)

"""
function cointegrate(ci::ChartInfo,evals,cachedir,a,b)
    mn = ci.mn
    asd = BSON.load("$cachedir/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/hd-$(mn[1])-$(mn[2]).bson")
    dm  = DataMap(()->asd,hd,ci.type,ci.indices,ci.priors,ci.base);
    di  = DataIntegral(dm,a,b)
    di(evals)

    (di.data, di.err)
end

"""
    integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::AbstractWorkerPool; rtol=1e-12, atol=1e-12, evals::AbstractRange, ranges::ChartInfo=ChartInfo(), cachedir::String)
"""
function integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::AbstractWorkerPool; rtol=1e-12, atol=1e-12, evals::AbstractRange, ranges::ChartInfo=ChartInfo(), cachedir::String=ENV["cachedir"])
    #Divide integration region
    nw = length(pool)
    as,bs = strip_divide(a,b,nw)
    #Hold global eval count fixed
    Nevals = Int.(round.(evals./nw))
    futures = Vector{Future}(undef, nw)

    for i=1:nw
        futures[i] = @spawnat(pool,cointegrate(ranges,Nevals,cachedir,as[i],bs[i]))
    end

    #initialize
    err = zeros(typeof(0.0),length(Nevals))
    data = Vector{typeof(dm.chart.data)}(undef,length(Nevals))
    for i=1:length(Nevals)
        data[i] = zero(dm.chart.data)
    end
    for (i,fut) ∈ enumerate(futures)
        wdata,werr = fetch(fut)
        #::Tuple{Vector{typeof(dm.chart.data)},Vector{Float64}}
        err .+= werr.^2
        for i ∈ 1:length(wdata)
            data[i] .+= wdata[i]
        end
    end
    err .= sqrt.(err)./(dm.chart.data|>length)

    di = DataIntegral(dm,a,b,ranges=ranges,cachedir=ENV["cachedir"])

    push!.(Ref(di.err), err)
    push!.(Ref(di.evals) , Int.(collect(Nevals).*nw))
    push!.(Ref(di.data) , data)

    di
end
#
# dm::DM
# a::AV
# b::AV
# times::EA
# err::EA
# evals::IA
# bytes::IA
# data::DA
# chartinfo::CIT
# cachedir::String

function (di::DataIntegral)(evals::Union{AbstractRange,Vector{Int64}}, pool::Symbol = :none, atol = 1e-20, rtol = 1e-20)
    for eval ∈ evals
        di(eval, pool, atol, rtol);
        nothing
    end
    di
end

function (di::DataIntegral)(evals::Int, pool::Symbol = :none, atol = 1e-20, rtol = 1e-20)::DataIntegral
    #integration line
    stats = @timed(integrate(di.dm, di.a, di.b, pool; evals = evals, rtol=rtol, atol=atol, ranges=di.chartinfo, cachedir=di.cachedir))

    #Store integration informaiton
    push!(di.evals, evals)
    push!(di.err,   stats.value)
    push!(di.times, stats.time/3600)
    push!(di.bytes, stats.bytes)
    push!(di.data,  copy(di.dm.chart.data))

    #print to stdout
    #println("did $(fieldnames(di.dm.chart|>typeof)[1]) integral w/ $(evals) pts in $(round(stats.time,sigdigits=3))s or $(round(stats.time/3600,sigdigits=3)) SUs");flush(stdout)
    #flush(stdout)

    di
end

# """
#     (di::DataIntegral)(evals::Union{AbstractRange,Int}, pool::AbstractWorkerPool, atol = 1e-20, rtol = 1e-20)::DataIntegral
# """
function (di::DataIntegral)(evals::Int64, pool::AbstractWorkerPool, atol = 1e-20, rtol = 1e-20)::DataIntegral
    integrate(di.dm, di.a, di.b, pool; evals = evals:evals, rtol=rtol, atol=atol, ranges=di.chartinfo, cachedir=di.cachedir)
end


##############################
#### DataSection
##############################

"""
    DataSection{DM <: DataMap, TC <: DataChart}

A structure to hold sampling information of DataMaps over the sampling base space
"""
struct DataSection{DM <: DataMap, TC <: DataChart}
    datadim::Int
    dm::DM
    chart::TC
end

"""
    combine_bases(oldbase,dombase)
"""
function combine_bases(oldbase,dombase)
    dimbase = length(oldbase)
    new_type = (oldbase[1]...,dombase[1]...)|>typeof

    new_base = Array{new_type,length(size(oldbase))+length(size(dombase))}(undef,(size(oldbase)...,size(dombase)...))
    for (i,b) ∈ enumerate(oldbase)
        for (j,d) ∈ enumerate(dombase)
            new_base[i+(j-1)*dimbase]=(b...,d...)
        end
    end
    new_base
end

"""
    prior_base_combine(oldpriors,oldbase)
"""
function prior_base_combine(oldpriors,oldbase)
    dimpriors  = length(oldpriors)
    new_type = (oldpriors[1]...,oldbase[1]...)|>typeof

    new_priors = Array{new_type,length(size(oldpriors))+length(size(oldbase))}(undef,(size(oldpriors)...,size(oldbase)...))
    for (i,p) ∈ enumerate(oldpriors)
        for (j,b) ∈ enumerate(oldbase)
            new_priors[i+(j-1)*dimpriors]=(p...,b...)
        end
    end
    new_priors
end

"""
    DataSection(dm::DataMap, dom::Vector{Tuple{Symbol,N,N,Int64}} where N <: Number)

A Data Section constructor with a range input
"""
function DataSection(dm::DataMap, dom::Vector{Tuple{Symbol,N,N,Int64}} where N <: Number)
    chart = dm.chart
    priors  = prior_base_combine(chart.priors,chart.base)
    dombase = range_scope(dom)
    DataSection(chart.l_i*chart.l_p*chart.l_b, dm, DataChart(dm.chart.name,chart.indices,priors,dombase,eltype(chart.data)))
end

"""
    DataSection(dm::DataMap, dom::AbstractArray)

DataSection constructor
"""
function DataSection(dm::DataMap, dom::AbstractArray)
    chart = dm.chart
    priors  = prior_base_combine(chart.priors,chart.base)

    DataSection(chart.l_i*chart.l_p*chart.l_b, dm,DataChart(dm.chart.name,chart.indices,priors,dom,eltype(chart.data)))
end

"""
    DataSection(;asd,mn,dtype,N,indices,priors,base,cachedir,offset, kargs...)

Dictionary based constructor
"""
function DataSection(;offset, N, kargs...)
    dom = [(:n1,-0.5+offset[1],0.5+offset[1],N),(:n2,-0.5+offset[2],0.5+offset[2],N)]
    dm = DataMap(; kargs...)

    DataSection(dm,dom)
end


# """
#     (ds::DataSection)()
# A call to sample the data map of the defined sampling space and thus form a section
# """
function (ds::DataSection)()
    for (i,b) ∈ enumerate(ds.chart.base)
        ds.chart.data[(1+ds.datadim*(i-1)):ds.datadim*i] .= ds.dm(b)[:]
    end
end



"""
    class_weights(w::Tuple{Symbol,Symbol})

Returns symmetrized square triangular integral mesh weights. The weight corresponds to the average over global triangularization schemes and give corners 1.5, edges 3, and bulk 6.
"""
function class_weights(w::Tuple{Symbol,Symbol})
    if     w == (:bulk,:bulk)
        return 6.0
    elseif w == (:bulk,:min)
        return 3.0
    elseif w == (:bulk,:max)
        return 3.0
    elseif w == (:min,:bulk)
        return 3.0
    elseif w == (:min,:min)
        return 1.5
    elseif w == (:min,:max)
        return 1.5
    elseif w == (:max,:bulk)
        return 3.0
    elseif w == (:max,:min)
        return 1.5
    elseif w == (:max,:max)
        return 1.5
    end
end


"""
    weight_classify(x::NTuple{N, Float64} where N ,mins::AbstractArray,maxes::AbstractArray,i::Int)

Return a tuple classifying the location of each coordinate in the grid
"""
function weight_classify(x::NTuple{N, Float64} where N ,mins::AbstractArray,maxes::AbstractArray,i::Int)
    abs(x[i]-mins[i]) < 1e-10 ? :min :
        abs(x[i]-maxes[i]) < 1e-10 ? :max :
            :bulk
end


"""
    edge_test(x::NTuple{N, Float64} where N ,mins::AbstractArray,maxes::AbstractArray)

"""
function edge_test(x,mins,maxes)
    ([weight_classify(x,mins,maxes,i) for i ∈ 1:length(x)]...,)
end

"""
    square_integration_weights(base)

Returns the weights associated to the base of a data section
"""
function square_integration_weights(base)

    axs = [sort(union(getindex.(base,i)[:])) for i ∈ 1:length(base[1])]
    mins  = [min(axs[i]...) for i ∈ 1:length(axs) ]
    maxes = [max(axs[i]...) for i ∈ 1:length(axs) ]

    edge_weights = edge_test.(base,Ref(mins),Ref(maxes))

    class_weights.(edge_weights)
end

"""
    base_ranges(ds::DataSection, i_b::Int)
"""
function base_ranges(ds::DataSection, i_b::Int)
    (1:(ds.chart.l_i*ds.chart.l_p)).+(i_b-1)*ds.chart.l_i*ds.chart.l_p
end

"""
    integral_data(ds::DataSection)

Returns the integral in an AbstractArray type that corresponds to that of the data in DataMap of ds.
"""
function integral_data(ds::DataSection)
    base_index_ranges = collect.(base_ranges.(Ref(ds),1:ds.chart.l_b));
    base_weights = square_integration_weights(ds.chart.base);

    ds_integral = zero(ds.chart.data[base_index_ranges[1]])
    for (i,brng) ∈ enumerate(base_index_ranges)
        ds_integral .+= ds.chart.data[brng]*base_weights[i]
    end

    di_copy = zero(ds.dm.chart.data)
    di_copy[:] .= ds_integral ./((sqrt(length(base_weights))-1)^2*6*(ds.dm.Ω/(2π)^size(ds.dm.Λ,1)))

    di_copy
end

"""
    integrate(ds::DataSection)

Returns DataIntegral with identical DataMap and the output of integrate(ds) pushed into the first data position
"""
function integrate(ds::DataSection)
    di = DataIntegral(ds.dm)
    push!(di.err,0.0)
    push!(di.times,0.0)
    push!(di.evals,size(ds.chart.base,1))
    push!(di.data,integral_data(ds))
    di
end

"""
    integrate!(di::DataIntegral, ds::DataSection)

In place push of data section integral to the passed DataIntegral
"""
function integrate!(di::DataIntegral, ds::DataSection)
    push!(di.err,0.0)
    push!(di.times,0.0)
    push!(di.evals,size(ds.chart.base,1))
    push!(di.data,integral_data(ds))
    di
end

###########################################
#### Importing and Exporting DataStructures
###########################################
"""
    data_export(datadir,strct)
"""
function data_export(datadir,strct)
    dict = OrderedDict{Symbol,Any}()
    push!(dict,:type=>(strct|>typeof))
    dictdata = OrderedDict{Symbol,Any}()
    for name ∈ strct|>typeof|>fieldnames
        push!(dictdata,name=>getfield(strct,name))
    end
    push!(dict,:data=>dictdata)
    bson(datadir,dict)
    return datadir
end

"""
    data_import(datadir)
"""
function data_import(datadir)
    dict = BSON.load(datadir,@__MODULE__)
    dict[:type]((dict[:data]|>values|>collect)...)
end
#
# function data_export(datadir,strct)
#     bson(datadir,Dict(:data=>strct))
# end
#
# function data_import(datadir)
#     BSON.load(datadir,@__MODULE__)[:data]
# end
