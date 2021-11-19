
########################################################################
#### Projectors
########################################################################
function sqo_info(asd)
    asdb = asd|>SolidState.ASDBasics
    names = asdb["sqo"][1]|>typeof|>fieldnames
    sqoi = Dict{Symbol, Any}()
    map(names) do name
        sqoi[name] = getfield.(asdb["sqo"],name)
    end
    sqoi
end

function sqo_unq(sqoi)
    unqsqoi = Dict{Symbol, Any}()
    map(sqoi|>keys|>collect) do name
        unqsqoi[name] = sqoi[name]|>unique
    end
    unqsqoi
end

function create_projectors(asd,name)
    sqoi = sqo_info(asd)
    sqou = sqo_unq(sqoi)

    map(sqou[name]) do idx
        Matrix(Diagonal(Complex.(Float64.(Ref(idx).==sqoi[name]))))
    end
end

struct Projector{LV <: AbstractArray, BV <: AbstractArray, CV <: AbstractArray, AA <: AbstractArray}
    name::Symbol
    labels::LV
    proj::BV
    ref::CV
    OD::AA
end

export Projector
function Projector(asd,feature::Symbol)
    sqoi = sqo_info(asd)
    dim  = length(collect(values(sqoi))[1])
    Projector(
        feature,
        sqo_unq(sqoi)[feature],
        create_projectors(asd,feature),
        create_projectors(asd,feature),
        ones(dim,dim)-Diagonal(ones(dim))
    )
end

function product_projector(P1,P2)
    name = Symbol(P1.name,"_",P2.name)

    labels = vcat(map(P1.labels) do p1
        map(P2.labels) do p2
            p1,p2
        end
    end...)

    projs = vcat(map(P1.proj) do p1
        map(P2.proj) do p2
            p1*p2
        end
    end...)

    Projector(name,labels,projs,projs)
end

function Projector(asd,names::Vector{Symbol})
    P = Projector(asd,names[1])
    for (i,name) in enumerate(names[2:end])
        P = product_projector(P,Projector(asd,name))
    end
    P
end

function conjugate_projectors!(P,hd)
    U = hd.h_ops.E.vectors
    for (i,π0) in enumerate(P.ref)
        P.proj[i] .= U*π0*(U')
    end
end
