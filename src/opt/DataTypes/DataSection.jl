
"""
    DataSection{DM <: DataMap, TC <: TensorChart}

A structure to hold sampling information of DataMaps over the sampling base space
"""
struct DataSection{DM <: DataMap, TC <: TensorChart}
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
"""
function DataSection(dm::DataMap, dom::Vector{Tuple{Symbol,N,N,Int64}} where N <: Number)
    chart = getfield(dm.chart,1)
    priors  = prior_base_combine(chart.priors,chart.base)
    dombase = range_scope(dom)
    DataSection(chart.l_i*chart.l_p*chart.l_b, dm, TensorChart(chart.indices,priors,dombase,eltype(chart.data)))
end

"""
    DataSection(dm::DataMap, dom::AbstractArray)

DataMap constructor
"""
function DataSection(dm::DataMap, dom::AbstractArray)
    chart = getfield(dm.chart,1)
    priors  = prior_base_combine(chart.priors,chart.base)

    DataSection(chart.l_i*chart.l_p*chart.l_b, dm,TensorChart(chart.indices,priors,dom,eltype(chart.data)))
end

"""
    DataSection(;asd,mn,dtype,N,indices,priors,base,cachedir,offset, kargs...)

Dictionary based constructor
"""
function DataSection(;asd,mn,dtype,N,indices,priors,base,cachedir,offset, kargs...)
    dom = [(:n1,-0.5+offset[1],0.5+offset[1],N),(:n2,-0.5+offset[2],0.5+offset[2],N)]
    dm = DataMap(asd, dtype, indices, priors, base; cachedir=cachedir, mn=mn)

    DataSection(dm,dom)
end


"""
    (ds::DataSection)()
A call to sample the data map of the defined sampling space and thus form a section
"""
function (ds::DataSection)()
    for (i,b) ∈ enumerate(ds.chart.base)
        ds.chart.data[(1+ds.datadim*(i-1)):ds.datadim*i] .= ds.dm(b)[:]
    end
end
