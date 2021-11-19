########################################
#### Projection Operators
########################################
struct Projector{AC <: AbstractArray, BV <: AbstractArray}
    proj0::AC
    proj::BV
end

function Projector(dim, Nfeature)
    Projector(
        zeros(typeof(0.0),(dim,Nfeature)),
        [zeros(typeof(Complex(0.0)),(dim,dim)) for _=1:Nfeature],
    )
end

#### Layer Projectors
function init_sqoi(asd,names)
    asdb = asd|>SolidState.ASDBasics
    sqoi = Dict{Symbol, Any}(:asdb=>asdb,:names=>names)
    map(names) do name
        sqoi[name] = getfield.(asdb["sqo"],name)
        sqoi[Symbol(:unq,name)] = sqoi[name]|>unique
    end
    sqoi
end

function init_layer_proj!(p::Projector,name)
    map(enumerate(p.sqoi[Symbol(:unq,name)])) do (i,a)
        inds = (a.==p.sqoi[name])
        getfield(p,Symbol(name,0))[i] .= inds
        # proj = getfield(p,name)[i]
        # proj[diagind(proj)[inds]] .= Complex(1.0)
        # proj
    end
end

function layer_projector(asd::Dict{String,Any})
    names = [:layer]
    sqoi = init_sqoi(asd,names)
    Nlayer = length.(getindex.(Ref(sqoi),Symbol.(:unq,names))...)|>length
    p = ProjectorLP(length(sqoi[:asdb]["sqo"]),Nlayer,sqoi)

    init_proj!.(Ref(p),names[1])

    p
end

#### Sublatttice Projectors

function init_proj!(p::Projector,name)
    map(enumerate(p.sqoi[Symbol(:unq,name)])) do (i,a)
        inds = (a.==p.sqoi[name])
        getfield(p,Symbol(name,0))[i] .= inds
        # proj = getfield(p,name)[i]
        # proj[diagind(proj)[inds]] .= Complex(1.0)
        # proj
    end
end

function sublattice_projector(asd::Dict{String,Any})
    sqoi = init_sqoi(asd,[:layer,:atom])
    p = ProjectorSL(length(sqoi[:asdb]["sqo"]),length(sqoi[Symbol(:unq,:layer)]),length(sqoi[Symbol(:unq,:atom)]),sqoi)

    init_proj!.(Ref(p),[:atom,:layer])
    p
end
