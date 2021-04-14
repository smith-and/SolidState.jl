"""
    integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::Symbol = :none; rtol=1e-12, atol=1e-12, evals::Int=typemax(Int), kwargs...)::Float64

Core Call to HCubature
"""
function integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::Symbol = :none; rtol=1e-12, atol=1e-12, evals::Int=typemax(Int), kwargs...)::Float64
    getfield(dm.chart,1).data .*= 0.0
    integral,int_error = hcubature(dm, a, b; norm=norm, rtol=rtol, atol=atol, maxevals=evals, initdiv=1)::Tuple{typeof(getfield(dm.chart,1).data),typeof(1.0)}
    getfield(dm.chart,1).data .= integral

    int_error/(getfield(dm.chart,1).data|>length)
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
    cointegrate(ci::ChartInfo,evals,cachedirr,a,b)

"""
function cointegrate(ci::ChartInfo,evals,cachedirr,a,b)
    mn = ci.mn
    asd = BSON.load("$cachedirr/asd-$(mn[1])-$(mn[2]).bson")
    hd  = data_import("$cachedirr/hd-$(mn[1])-$(mn[2]).bson")
    dm  = DataMap(ci.type,asd,hd,ci.indices,ci.priors,ci.base);
    di  = DataIntegral(dm,a,b)
    di(evals)

    (di.data, di.err)
end

"""
    integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::AbstractWorkerPool; rtol=1e-12, atol=1e-12, evals::AbstractRange, ranges::ChartInfo=ChartInfo(), cachedirr::String)
"""
function integrate(dm::DM where DM <: DataMap, a::AbstractVector, b::AbstractVector, pool::AbstractWorkerPool; rtol=1e-12, atol=1e-12, evals::AbstractRange, ranges::ChartInfo=ChartInfo(), cachedirr::String)
    #Divide integration region
    nw = length(pool)
    as,bs = strip_divide(a,b,nw)
    #Hold global eval count fixed
    Nevals = Int.(round.(evals./nw))
    futures = Vector{Future}(undef, nw)
    for i=1:nw
        futures[i] = @spawnat(pool,cointegrate(ranges,Nevals,cachedirr,as[i],bs[i]))
    end

    #initialize
    err = zeros(typeof(0.0),length(Nevals))
    data = Vector{typeof(getfield(dm.chart,1).data)}(undef,length(Nevals))
    for i=1:length(Nevals)
        data[i] = zero(getfield(dm.chart,1).data)
    end

    for (i,fut) ∈ enumerate(futures)
        wdata,werr = fetch(fut)
        #::Tuple{Vector{typeof(getfield(dm.chart,1).data)},Vector{Float64}}
        err .+= werr.^2
        for i ∈ 1:length(wdata)
            data[i] .+= wdata[i]
        end
    end
    err .= sqrt.(err)./(getfield(dm.chart,1).data|>length)

    DataIntegral(dm,a,b,zeros(eltype(err),length(Nevals)),err,Int.(collect(Nevals).*nw),zeros(Int,length(Nevals)),data,ranges,cachedirr)
end


### DataSection Integration

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

    di_copy = zero(getfield(ds.dm.chart,1).data)
    di_copy[:] .= ds_integral

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
