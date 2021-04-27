using OrderedCollections,LinearAlgebra,SolidState,BSON
using SolidState: SHG, ChartInfo, cθ

# datatype=SHG
# indices =[(2,2,2)]
# base    =[(:ω,0.0,10.0,400)]
function mem_test(ASD::Symbol, comargs::Vector{Tuple{Int64,Int64}}, datatype, indices, base; cachedir=cachedir)
    priors  =[(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)]
    a       =[0.0,0.0]
    b       =[1.0,1.0]

    com_dict = OrderedDict{Int,OrderedDict{Symbol,Any}}()
    for (i,mn) ∈ enumerate(comargs)
        di = DataIntegral(
                DataMap(
                    datatype,
                    BSON.load("$cachedir/$ASD/asd-$(mn[1])-$(mn[2]).bson"),
                    data_import("$cachedir/$ASD/hd-$(mn[1])-$(mn[2]).bson"),
                    indices,priors,base
                ),
                a,b;
                ranges=ChartInfo(datatype,indices,priors,base,mn),cachedir=cachedir*"/$ASD"
        )

        di.dm(rand(2))

        top_out = read(`top -bn1 -p $(getpid())`, String)
        res_size = split(split(top_out,  "\n")[end-1])[6]

        sizes = OrderedDict{Symbol,Any}(
            :mn     => mn,
            :angle  => cθ(mn...)*180/π,
            :dim    => di.dm.d,
            #:asd    => Base.summarysize(asd)/1e9,
            #:hd     => Base.summarysize(hd)/1e9,
            :chart  => Base.summarysize(getfield(di,1).chart)/1e9,
            :di     => Base.summarysize(di)/1e9,
            :top    => res_size,
        )

        println("")
        println("sizes for $mn: ");flush(stdout)
        for p ∈ sizes
            println("$(p.first): $(p.second)")
        end

        push!(com_dict, i=>sizes)
    end

    bson("$(@__DIR__)/$ASD-mem_test.bson", com_dict)
end

function mem_test(asd, comargs::Tuple{Symbol,Int,Int})

    cachedir = "$(@__DIR__)/../.cache"

    SolidState.make_models(asd, comargs..., cachedir=cachedir)

    mem_test(asd, SolidState.twist_series(comargs...), cachedir=cachedir)
end
