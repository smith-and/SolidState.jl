#DataStructures
include("DataTypes/DataChart.jl")
include("DataTypes/DataMap.jl")
include("DataTypes/DataSection.jl")
# include("DataTypes/DataIntegral.jl")

# Importing and Exporting DataStructures
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
    dict = BSON.load(datadir)
    dict[:type]((dict[:data]|>values|>collect)...)
end
