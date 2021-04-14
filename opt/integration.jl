
function integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector; rtol=1e-12, atol=1e-12, maxevals=typemax(Int))::Float64
    getfield(dm.chart,1).data .*= 0.0
    integral,int_error = hcubature(dm, a, b; norm=norm, rtol=rtol, atol=atol, maxevals=maxevals, initdiv=1)
    getfield(dm.chart,1).data .= integral
    int_error/(getfield(dm.chart,1).data|>length)
end

function integrate(datatype::Type{T} where T <: DataChart,asd::Dict{String,Any},hd::HamiltonianDensity,indices,priors,base,evals; tag="", rootdir=pwd(),a=[0.0,0.0],b=[1.0,1.0])
    dm  = DataMap(datatype,asd,hd,indices,priors,base);
    di  = DataIntegral(dm,a,b);
    di(evals);
    data_export("$rootdir/$(fieldname(datatype,1))$tag.bson",di);
end


function strip_divide(a,b,nw)
    as = [[a[1],b[2]*(i-1)/nw] for i ∈ 1:nw]
    bs = [[b[1],b[2]*i/nw] for i∈1:nw]
    return as,bs
end

function cointegrate(dm0::DM, a::AbstractVector, b::AbstractVector; rtol=1e-12, atol=1e-12, maxevals=typemax(Int))::Tuple{DM,Float64} where DM <: DataMap
    dm = deepcopy(dm0)
    integral,int_error = hcubature(dm, a, b; norm=norm, rtol=rtol, atol=atol, maxevals=maxevals, initdiv=1)
    getfield(dm.chart,1).data .= integral
    (dm,int_error)

end

function integrate(dm::DM, a::AbstractVector, b::AbstractVector, worker_pool::Union{WorkerPool,CachingPool}; rtol=1e-12, atol=1e-12, maxevals=typemax(Int))::Float64 where DM <: DataMap
    getfield(dm.chart,1).data .*= 0.0
    nw = length(worker_pool)
    as,bs = strip_divide(a,b,nw)
    futures = Vector{Future}(undef, nw)
    # Send strip integrals out
    for i=1:nw
        futures[i] = remotecall(cointegrate,worker_pool,dm,as[i],bs[i];maxevals=maxevals)
    end

    # Fetch the worker data
    err = 0.0
    for (i,data) ∈ enumerate(futures)
        wdata,werr = fetch(data)::Tuple{DM,Float64}
        err += werr^2
        getfield(dm.chart,1).data .+= getfield(wdata.chart,1).data
    end
    err =sqrt(err)
    err/(getfield(dm.chart,1).data|>length)
end
