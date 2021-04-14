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
function DataIntegral(dm::DM, a::AV=[0.0,0.0], b::AV=[1.0,1.0];cachedir::String="none", ranges::ChartInfo=ChartInfo())::DataIntegral where {AV <: AbstractVector, DM <: DataMap, PT <: Union{AbstractWorkerPool,Symbol}}
    DataIntegral(dm,a,b,typeof(1.0)[],typeof(1.0)[],typeof(1)[],typeof(1)[],typeof(getfield(dm.chart,1).data)[],ranges,cachedir)
end

#Construction with map description
"""
    DataIntegral(datatype::Type{T} where T <: DataChart,asd::Dict{String,Any},hd::HamiltonianDensity,indices,priors,base; a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector, DM <: DataMap}
"""
function DataIntegral(datatype::Type{T} where T <: DataChart,asd::Dict{String,Any},hd::HamiltonianDensity,indices,priors,base; a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector, DM <: DataMap}
    dm  = DataMap(datatype,asd,hd,indices,priors,base);
    DataIntegral(dm,a,b, ranges=ChartInfo(datatype,indices,priors,base,(0,0)),cachedir=pwd())
end

"""
    DataIntegral(datatype::Type{T} where T <: DataChart,indices,priors, base, cachedir, mn; a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector}
Construction with cached objects & description
"""
function DataIntegral(datatype::Type{T} where T <: DataChart,indices,priors, base, cachedir, mn; a::AV=[0.0,0.0], b::AV =[1.0,1.0])::DataIntegral  where {AV <: AbstractVector}
    asd = BSON.load("$cachedir/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedir/hd-$(mn[1])-$(mn[2]).bson")
    dm = DataMap(datatype,asd,hd,indices,priors,base);
    DataIntegral(dm,a,b; ranges=ChartInfo(datatype,indices,priors,base,mn),cachedir)
end

"""
    DataIntegral(ci::ChartInfo, cachedir::String; a::AbstractVector=[0.0,0.0], b::AbstractVector =[1.0,1.0])::DataIntegral
"""
function DataIntegral(ci::ChartInfo, cachedir::String; a::AbstractVector=[0.0,0.0], b::AbstractVector =[1.0,1.0])::DataIntegral
    asd = BSON.load("$cachedir/asd-$(ci.mn[1])-$(ci.mn[2]).bson")
    hd  = data_import("$cachedir/hd-$(ci.mn[1])-$(ci.mn[2]).bson")
    dm = DataMap(datatype,asd,hd,indices,priors,base);
    DataIntegral(dm,a,b; ranges=ci,cachedir)
end


#Pool Routines
include(".integrate.jl")

function (di::DataIntegral)(evals::Union{AbstractRange,Vector{Int64}}, pool::Symbol = :none, atol = 1e-20, rtol = 1e-20)
    for eval ∈ evals
        di(eval, pool, atol, rtol);
        nothing
    end
    di
end

# """
#     (di::DataIntegral)(evals::Int, pool::Symbol = :none, atol = 1e-20, rtol = 1e-20)::DataIntegral
#
# """
function (di::DataIntegral)(evals::Int, pool::Symbol = :none, atol = 1e-20, rtol = 1e-20)::DataIntegral
    #integration line
    stats = @timed(integrate(di.dm, di.a, di.b, pool; evals = evals, rtol=rtol, atol=atol, ranges=di.chartinfo, cachedirr=di.cachedir))

    #Store integration informaiton
    push!(di.evals, evals)
    push!(di.err,   stats.value)
    push!(di.times, stats.time/3600)
    push!(di.bytes, stats.bytes)
    push!(di.data,  copy(getfield(di.dm.chart,1).data))

    #print to stdout
    #println("did $(fieldnames(di.dm.chart|>typeof)[1]) integral w/ $(evals) pts in $(round(stats.time,sigdigits=3))s or $(round(stats.time/3600,sigdigits=3)) SUs");flush(stdout)
    #flush(stdout)

    di
end

"""
    di(evals::Union{AbstractRange,Int}, pool::AbstractWorkerPool, atol = 1e-20, rtol = 1e-20)::DataIntegral
"""
function (di::DataIntegral)(evals::Union{AbstractRange,Int}, pool::AbstractWorkerPool, atol = 1e-20, rtol = 1e-20)::DataIntegral
    integrate(di.dm, di.a, di.b, pool; evals = evals, rtol=rtol, atol=atol, ranges=di.chartinfo, cachedirr=di.cachedir)
end
