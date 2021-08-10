#######################################
#### Generalize Slater Koster Functions
#######################################

#Generalize Slater Koster Hopping Functions
#radial functions
function rbinomial(n::Int,k::Int)
    n<0 ? (-1)^k*gamma(k-n)/(gamma(k+1)*gamma(-n)) : gamma(n+1)/(gamma(k+1)*gamma(n-k+1))
end

function amid(l::Int,m::Int)
    (l,m)==(0,0) ? 1 : sum(2*s+1 for s=0:max((l-1),0))+m+l+1
end

function constreg(a,d,v,e)
    function(r)
        norm(r)>1e-10 ? v*a : e
    end
end

function constregfunctions(energyscales)
    refcouplingfs = map((a,d,v,e)->constreg(a,d,v,e),energyscales.as,energyscales.ds,energyscales.vs,energyscales.es);
    function(l1, l2, m1, m2,r)
        refcouplingfs[amid(l1,m1),amid(l2,m2)](r)
    end
end


function expreg(a,d,v,e)
    function(r)
        norm(r)>1e-10 ? v*exp(-(norm(r)-a)/d) : e
    end
end

function expregfunctions(energyscales)
    refcouplingfs = map((a,d,v,e)->expreg(a,d,v,e),energyscales.as,energyscales.ds,energyscales.vs,energyscales.es);
    function(l1, l2, m1, m2,r)
        refcouplingfs[amid(l1,m1),amid(l2,m2)](r)
    end
end

function covreg(a,d,v,e)
    function(r)
        norm(r)>1e-10 ? v/(1+(norm(r)-a)/d)^3 : e
    end
end

function covregfunctions(energyscales)
    refcouplingfs = map((a,d,v,e)->covreg(a,d,v,e),energyscales.as,energyscales.ds,energyscales.vs,energyscales.es);
    function(l1, l2, m1, m2,r)
        refcouplingfs[amid(l1,m1),amid(l2,m2)](r)
    end
end;

#spin weighted spherical harmonics
function infTest(x)
    abs(x)==Inf ? 0 : x
end

function swhNorm(l,m,s)
    sqrt((2l+1)/(4π)*gamma(l+m+1)*gamma(l-m+1)/gamma(l+s+1)/gamma(l-s+1))
end

function swhPolar(l,m,s,θ)
    (abs(m)>l)||(abs(s)>l) ? 0 :
        begin
            x=0
            for r=0:1:(l-s)
                x+=binomial(l-s,r)*binomial(l+s,r+s-m)*(-1)^r*infTest((sin(θ/2))^(2l-(2r+s-m)))*(cos(θ/2))^(2r+s-m)
            end
            x
        end
end

function sYlm(l,m,s)
    function(θ,ϕ)
        ((abs(m)>l)||(abs(s)>l) ? 0 : (exp(complex(0,m*ϕ))*swhNorm(l,m,s)*(-1)^(l-s+m)*swhPolar(l,m,s,θ)))
    end
end
#Generalize Slater Koster Functions
function gskAngular(m1,m2,θ,ϕ,l,lp,m,mp)
    4π/sqrt((2*l+1)*(2*lp+1))*(-1)^(m+m1)* exp(complex(0,(m1-m2)*ϕ))*sYlm(l,-m,m1)(θ,ϕ)*sYlm(lp,mp,-m2)(θ,ϕ)
end

function gsk(radials)
    function(l,m,lp,mp)
        function(r3d)
            x=0; for  m1=-l:1:l, m2=(-lp):1:lp
                x+=radials(l,lp,m1,m2,norm(r3d))*gskAngular(m1,m2,atan(sqrt(r3d[1]^2+r3d[2]^2),r3d[3]),atan(r3d[2],r3d[1]),l,lp,m,mp)
            end
            x
        end
    end
end

#=

function cJ(J,m1,m2,l,lp,m,mp)
    (abs(mp-m)>J)|| (abs(m2-m1)>J) ? 0 : (-1)^(m+m1)*sqrt(4*π*(2*J+1))*wigner3j(l,lp,J,-m,mp,mp-m)*wigner3j(l,lp,J,m1,-m2,m2-m1)
end

function gskAngular0(m1,m2,θ,ϕ,l,lp,m,mp)
    x=0.0
    for  J=abs(l-lp):1:(l+lp)
        x+=cJ(J,m1,m2,l,lp,m,mp)*sYlm(J,mp-m,m1-m2)(θ,ϕ)
    end
    x*exp(complex(0,(m1-m2)*ϕ))
end

=#

#Local Hopping parameterization from ASD
function ASDGSK(ASD::Dict{String,Any})::Dict{NTuple{N,String} where N,Any}
    if ASD["regulator"]=="cov"
        return Dict(x=>gsk(covregfunctions(ASD["scales"][x])) for x=keys(ASD["scales"]))
    elseif  ASD["regulator"]=="exp"
        return Dict(x=>gsk(expregfunctions(ASD["scales"][x])) for x=keys(ASD["scales"]))
    elseif ASD["regulator"] == "const"
        return Dict(x=>gsk(constregfunctions(ASD["scales"][x])) for x=keys(ASD["scales"]))
    end
end

#######################################
#### Tightbinding Model
#######################################

unqDict(unq)=Dict{UInt64,Any}(hash(unq[i])=>i for i=1:length(unq))
unqht(collection)=unique(collection) |> x->(x,unqDict(x))

#overall calculation of the information in the Hamiltonian Graph
struct TightBindingInfo{Int_Type <: Integer, Float_Type  <: Real} #<: HamiltonianDensity
    index_rep::Array{Int_Type, 2}
    weights::Vector{Vector{Complex{Float_Type}}}
    covers::Vector{Vector{Vector{Float_Type}}}
    pz::Array{Complex{Float_Type},2}
    pjs::Vector{Array{Complex{Float_Type},2}}
    slv_pjs::Vector{Array{Complex{Float_Type},2}}
end

#asos
aso(x)= (atom=(x.atom),spin=(x.spin),orbital=(x.orbital))
aso(x,y)= (atom=(x.atom,y.atom),spin=[x.spin...,y.spin...],orbital=[x.orbital...,y.orbital...])
asoreg(asdb)=function(x)asdb["gsk"][x.atom](x.orbital...)end
asounqhash(unq)= broadcast(aso, unq, reshape(unq,(1,length(unq)))) |> unqht;

function asoInformation(asdb::Dict{String,Any})::(Tuple{Vector{Function},Function})
    #list of atomic, orbital, spin information for 1 body operators
    asos=aso.(asdb["sqo"])
    asosunq=unique(asos)
    asosIndx=indexin(asos,asosunq);

    #forming aso index tables for 2-body operators
    unqaso, asoht = asounqhash(asosunq);
    asoidxmap = Array{Any,2}(undef,(length(asosunq),length(asosunq)));
    asoidxmap.=hash.(broadcast(aso, asosunq, reshape(asosunq,(1,length(asosunq)))))
    for i=1:length(asosunq), j=1:length(asosunq) asoidxmap[i,j]=asoht[asoidxmap[i,j]] end
    asoIndexer(i,j)= asoidxmap[asosIndx[i],asosIndx[j]];

    #generalized slater koster interactions for two body interactions
    regf=asoreg(asdb);
    regs=regf.(unqaso);

    (regs,asoIndexer)
end

#sublattice info
function cutpush(cut,c0,w)
    if cut(w)
        push!(c0,w)
    end
end

function xcover(asdb,vr)
    lbase=asdb["lbase"];
    cut=asdb["cutfunction"];
    w=Array{Float64,1}(undef,(3))
    function coverer(v::Nothing)
        Array{Float64,1}[]
    end
    function coverer(v::Array{Float64,1})::Array{Array{Float64,1},1}
        c=Array{Float64,1}[]
        #(w.= lbase[i].+v;w)
        for i=1:length(lbase) cutpush(cut,c, lbase[i].+v)  end
        return c #vr.(c)
    end

    return coverer
end

function sublatticeInformation(asdb::Dict{String,Any})::(Tuple{Vector{Vector{Vector{Float64}}},Function})
    dim=length(asdb["sqo"]);
    #vector rounding function
    vr=vrd(18);

    #determining sublattice positions of operators in the basis
    τs= Array{Array{Float64,1},1}(undef,(dim)); for i=1:dim τs[i] = asdb["sqo"][i].pos end;
    τunq=unique(τs);#unique set (useful when there are many orbitals per site)
    τIdxs=indexin(τs,τunq)#index of basis sublattice in the unique set

    #determining the sublattice differences
    τdim=length(τunq);
    dτs=Array{Any,2}(nothing,(τdim,τdim));
    #for i=1:τdim, j=i:τdim dτs[i,j]=vr(τunq[i]-τunq[j]) end
    for i=1:τdim, j=i:τdim dτs[i,j]=(τunq[i]-τunq[j]) end

    #construction of the basis indexer for unique sublattice differences
    dτunq, dτmap=unqht(dτs);
    dτIdxs=Array{Any,2}(nothing,(τdim,τdim));
    for i=1:τdim, j=i:τdim  dτIdxs[i,j]=dτmap[hash(dτs[i,j])] end
    function dτIndexer(i::Int64,j::Int64)::Int64
        dτIdxs[τIdxs[i],τIdxs[j]]
    end

    #lattice cover and cut function
    xc=xcover(asdb,vr);
    #calculating the lattice cover of the unique sublattice differences
    covers=xc.(dτunq);

    (covers,dτIndexer)
end

#product indexing
function productIndexing(elementIndexer::Function,dim::Int64)
    #forming product index table
    elTable =Array{Any,2}(nothing,(dim,dim));
    for i=1:dim, j=i:dim elTable[i,j]= elementIndexer(i,j)  end;

    #index table & list for the unique products
    eunq, emap=unqht(elTable);
    efind=(k->emap[hash(k)]);

    index_table = zeros(Int64,(dim,dim));
    for i=1:dim, j=i:dim index_table[i,j] = efind(elTable[i,j])::Int64  end;

    (index_table,eunq)
end

#Generalized Slater-Koster weight & cover map over index space.

struct weightcover{Float_Type <: Real , Func_Type <: Function}
    covers::Vector{Vector{Vector{Float_Type}}}
    regs::Vector{Func_Type}
end

function (wc::weightcover)(x::Nothing)::Tuple{Array{Complex{Float64},1},Vector{Vector{Float64}}}
    ([Complex(0.0)],[[0.0,0.0,0.0]])::Tuple{Array{Complex{Float64},1},Vector{Vector{Float64}}}
end

function (wc::weightcover)(x::Tuple{Int64,Int64})::Tuple{Array{Complex{Float64},1},Vector{Vector{Float64}}}
    ((wc.regs[x[1]].(wc.covers[x[2]]))::Array{Complex{Float64},1},wc.covers[x[2]]::Vector{Vector{Float64}})
end

export form_model
function form_model(asd,wha)
    merge(asd,wha)
end

function form_model(asd,wha,p)
    try
        ASD = asd(;p...)
        form_model(ASD,wha(ASD;p...))
    catch
        ASD = asd()
        form_model(ASD,wha(ASD;p...))
    end
end

function modify_asd(asd,key,val)
    if :layer==key
        LayerASD(asd,val)
    elseif :shift==key
        LayerShiftASD(asd,val)
    elseif :twist==key
        CommensurateASD(asd,val)
    else
        asd
    end
end

function form_model(p::OrderedDict)
    asd = form_model(p[:asd],p[:wha],p)
    for key in keys(p)
        asd = modify_asd(asd,key,p[key])
    end
    asd
end


function tb_info(asd)

    asdb=ASDBasics(asd);
    push!(asdb, "gsk"=>ASDGSK(asd))
    #atomic site orbital information
    (regs,asoIndexer)=asoInformation(asdb)

    #sublattice difference information
    (covers,dτIndexer)=sublatticeInformation(asdb)
    #constructing the unique product element and an index table for the two body representation
    elementIndexer(i,j) = covers[dτIndexer(i,j)]==[] ? nothing : (asoIndexer(i,j),dτIndexer(i,j))
    (elTable,eunq)=productIndexing(elementIndexer,length(asdb["sqo"]));

    index_table = fill(0,size(elTable))
    index_table .= elTable;
    #function casting covers as cover->(weight,cover) where the weights at gsk(aso) for the element
    wc=weightcover(covers,regs);
    #forming functions over momentum label for the unqiue elements k->[unique matrix elements...]
    weights_covers = wc.(eunq);
    weights_unq = (x->x[1]).(weights_covers);
    covers_unq = (x->x[2]).(weights_covers);
    #subspace projectors

    layer_indexes = getindex.(asd["sites"],asd["sk"]["Layer"])
    l_pjs= [zeros(Complex{Float64},size(elTable)) for _=1:max(layer_indexes...) ]
    for i=1:length(asdb["sqo"])
        l_pjs[asdb["sqo"][i].layer][i,i] = Complex(1.0,0.0)
    end

    #layer polarization
    z_center= unique(layer_indexes) |> x->sum(x)/length(x)
    pz = Array{Complex{Float64},2}(undef, size(index_table));
    pz .= Diagonal([asdb["sqo"][i].layer - z_center for i=1:length(asdb["sqo"])])
    #pz .= Diagonal([i <= length(asd["sites"])/2 ? 1.0 : -1.0 for i=1:length(asd["sites"])]);

    #vector rounding function
    vr=vrd(18);
    dim=2;
    h_dim=length(asdb["sqo"]);
    #determining sublattice positions of operators in the basis
    τs= Array{Array{Float64,1},1}(undef,(h_dim)); for i=1:h_dim τs[i] = asdb["sqo"][i].pos end;
    τunq=unique(τs);#unique set (useful when there are many orbitals per site)
    τIdxs=indexin(τs,τunq)#index of basis sublattice in the unique set
    τ_dim = length(τunq)

    slv_partition = [ findall(idx->(idx==τ_idx), τIdxs) .|> x -> (x,x)  for τ_idx=1:τ_dim]
    slv_pjs = [zeros(ComplexF64,(h_dim,h_dim)) for _=1:1]

    (index_table, weights_unq, covers_unq, pz, l_pjs, slv_pjs)
end

function tb_info(asd,wha)
    tb_info(merge(asd,wha))
end

function tb_info(asd, wha, p)
    tb_info(form_model(asd,wha,p))
end

function TightBindingInfo(p::OrderedDict)
    TightBindingInfo(tb_info(form_model(p))...)
end


function TightBindingInfo(asd)
    TightBindingInfo(tb_info(asd)...)
end

function TightBindingInfo(asd,wha)
    TightBindingInfo(tb_info(asd,wha)...)
end

function TightBindingInfo(asd,wha,p)
    TightBindingInfo(tb_info(asd,wha,p)...)
end

#type for information of a single matrix element for Bundle and Tangent Bundle and Stability Bundle.
struct TBElement{Float_Type <: Real}
    h::Vector{Complex{Float_Type}}
    v::Vector{Complex{Float_Type}}
    a::Array{Complex{Float_Type},2}
end


#Construct Zero element
const TBE_Zero  = ([Complex(0.0)],[Complex(0.0),Complex(0.0)], [Complex(0.0) Complex(0.0) ; Complex(0.0) Complex(0.0)])
function TBElement()
    TBElement(TBE_Zero...)
end

#Method for TBElement type to zero internal data
const TBE_Zeros = ((Complex(0.0)),(Complex(0.0),Complex(0.0)), (Complex(0.0), Complex(0.0), Complex(0.0), Complex(0.0)))

function (tbe::TBElement{Float64})()::Nothing
    tbe.h[1]   = Complex(0.0)
    tbe.v[1]   = Complex(0.0)
    tbe.v[2]   = Complex(0.0)
    tbe.a[1,1] = Complex(0.0)
    tbe.a[2,2] = Complex(0.0)
    tbe.a[1,2] = Complex(0.0)
    tbe.a[2,1] = Complex(0.0)
    nothing
end

#Type which handles a TBElement being mapped over the base manifold
# and also has information about the indices of which matrix elements
#are being calculated by the object
struct TBEFunction{Int_Type <: Int, Float_Type <: Real}
    element::TBElement{Float_Type};
    weights::Vector{Complex{Float_Type}}
    cover::Vector{Vector{Float_Type}}
    index_domain::Vector{Int_Type}
    nw::Int64
end

function  TBEFunction(weights::Vector{Complex{Float64}}, cover::Vector{Vector{Float64}}, index_domain::Vector{Int64})
    TBEFunction(TBElement(),weights,cover, index_domain, length(weights))
end

@inline function int_f_3(tbe::TBElement, expΦw::ComplexF64, expΦwc1::ComplexF64, expΦwc2::ComplexF64, expΦwc12::ComplexF64, expΦwc11::ComplexF64, expΦwc22::ComplexF64)::Nothing
    @inbounds @fastmath begin
        tbe.h[1]   +=    expΦw
        tbe.v[1]   += im*expΦwc1
        tbe.v[2]   += im*expΦwc2
        tbe.a[1,1] +=   -expΦwc11
        tbe.a[2,2] +=   -expΦwc22
        tbe.a[1,2] +=   -expΦwc12
        tbe.a[2,1] +=   -expΦwc12
    end
    nothing
end

@inline function int_f_2(tbe::TBElement, expΦw::ComplexF64, expΦwc1::ComplexF64, expΦwc2::ComplexF64, c1::Float64, c2::Float64)::Nothing
    @inbounds @fastmath int_f_3(tbe, expΦw, expΦwc1, expΦwc2, expΦwc2*c1, expΦwc1*c1, expΦwc2*c2)
end

@inline function int_f_1(tbe::TBElement, expΦw::ComplexF64, c1::Float64, c2::Float64)::Nothing
    @fastmath int_f_2(tbe, expΦw, expΦw*c1, expΦw*c2, c1, c2)
end


@inline function int_f_0(tbe::TBElement, w::ComplexF64, k::Vector{Float64}, c1::Float64, c2::Float64)::Nothing
    @inbounds @fastmath int_f_1(tbe, exp(complex(0.0,k[1]*c1+k[2]*c2))*w, c1, c2)
end

@inline function (tbe::TBElement)(w::ComplexF64, c::Vector{Float64}, k::Vector{Float64})::Nothing
    @inbounds int_f_0(tbe, w, k, c[1], c[2])
end


@inline function (tbe_f::TBEFunction)(k::Vector{T} where T <: Real )::TBElement
    tbe_f.element();
    @inbounds for i=1:tbe_f.nw
        tbe_f.element(tbe_f.weights[i],tbe_f.cover[i],k)::Nothing
    end
    return tbe_f.element
end

#Function to create Vector{TBEFunction} for a given TightBindingInfo object.
function element_functions(tbi::TightBindingInfo)
    tbef_unq = Vector{TBEFunction{Int64,Float64}}(undef,size(tbi.weights))
    index_dom(i::Int64)::Vector{Int64} = findall(x-> x===i,tbi.index_rep[:]) ;
    for i=1:length(tbi.weights)
        tbef_unq[i] = TBEFunction(tbi.weights[i],tbi.covers[i],index_dom(i))
    end
    return tbef_unq
end

#Type to facilitate the Bundle map from a base manifold into the fiber associated with HamiltonianOperators
export TightBindingDensity
struct TightBindingDensity{H <: HamiltonianOperators, Int_Type <: Int, Float_Type <: Real} <: HamiltonianDensity{Int_Type,Float_Type}
    el_f::Vector{TBEFunction{Int_Type,Float_Type}}
    h_ops::H
end

#constructor for TightBindingDensity tight-binding function from the tight binding information, i.e. some static graph data for a model
function TightBindingDensity(tbi::TightBindingInfo; style=:normal)
    TightBindingDensity(element_functions(tbi), HamiltonianOperators(size(tbi.index_rep,1),style = style))
end

#constructor for TightBindingDensity tight-binding function from the tight binding information, i.e. some static graph data for a model
function TightBindingDensity(p::OrderedDict; style=:normal)
    TightBindingDensity(TightBindingInfo(p),style = style)
end

function TightBindingDensity(asd::Dict{String,Any}; style=:normal)
    TightBindingDensity(TightBindingInfo(asd),style = style)
end

function TightBindingDensity(asd::Dict{String,Any}, wha::Dict{String,Any}; style=:normal)
    TightBindingDensity(TightBindingInfo(asd,wha),style = style)
end

function TightBindingDensity(asd::Function, wha::Function, p::OrderedDict; style=:normal)
    TightBindingDensity(TightBindingInfo(asd,wha,p),style = style)
end

#function assignment
@inline @inbounds function rep_assign_core(idx::Int64, rep::HamiltonianOperators,vals::TBElement{Float64})::Nothing
    rep.h[idx]      = vals.h[1];
    rep.v[1][idx]   = vals.v[1];
    rep.v[2][idx]   = vals.v[2];
    rep.a[1,1][idx] = vals.a[1,1];
    rep.a[2,1][idx] = vals.a[2,1];
    rep.a[1,2][idx] = vals.a[1,2];
    rep.a[2,2][idx] = vals.a[2,2];
    nothing
end

@inline function index_domain_assign(index_domain, h_ops, tbe::TBElement)
    @inbounds @simd for index ∈ index_domain
        rep_assign_core(index, h_ops, tbe)
    end
end

#method on the TightBindingDensity to obtain the fibers of a Principal G Bundle, its Tangent Bundles and its Stability Bundle
@inline function (tbf::TightBindingDensity)(k::Vector{Float64})::Nothing
    @inbounds @simd for tbe_f ∈ tbf.el_f
        index_domain_assign(tbe_f.index_domain, tbf.h_ops, tbe_f(k))
    end
end
