#asos
aso(x)= (atom=(x.atom),spin=(x.spin),orbital=(x.orbital))
aso(x,y)= (atom=(x.atom,y.atom),spin=[x.spin...,y.spin...],orbital=[x.orbital...,y.orbital...])
asoreg(asdb)=function(x) asdb["gsk"][x.atom](x.orbital...) end
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


#function assignment
#=
@inline @inbounds function rep_assign_core(idx::Int64, rep::HamiltonianOperators{Float64},vals::TBElement{Float64})::Nothing
    rep.h[idx]      = vals.h[1];
    rep.v[1][idx]   = vals.v[1];
    rep.v[2][idx]   = vals.v[2];
    rep.a[1,1][idx] = vals.a[1,1];
    rep.a[2,1][idx] = vals.a[2,1];
    rep.a[1,2][idx] = vals.a[1,2];
    rep.a[2,2][idx] = vals.a[2,2];
    nothing
end
=#

function tb_info(asd)

    asdb=ASDBasics(asd);
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

    #=projector_diags = [[ j*(length(asd["sites"])/get(asd,"bd",size(elTable,1))) < i <= (j+1)*(length(asd["sites"])/get(asd,"bd",size(elTable,1))) ? 1.0 : 0.0 for i=1:length(asd["sites"])] for j=0:(get(asd,"bd",size(elTable,1))-1)]
    for j=1:length(l_pjs)
        for i=1:length(projector_diags[1])
            l_pjs[j][i,i] = projector_diags[j][i]
        end
    end=#
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
    #=for i=1:τ_dim
        for j=1:length(slv_partition[i])
            setindex!(slv_pjs[i], 1.0, slv_partition[i][j]...)
        end
    end=#

#=
    @inbounds for i=eachindex(tbf.el_f)
        @inbounds for j=eachindex(tbf.el_f[i].index_domain)
            rep_assign_core(tbf.el_f[i].index_domain[j], tbf.h_ops, tbf.el_f[i](k))
        end
    end
=#


    (index_table, weights_unq, covers_unq, pz, l_pjs, slv_pjs)
end
