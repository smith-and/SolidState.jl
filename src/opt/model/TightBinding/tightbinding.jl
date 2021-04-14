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

include("tightbinding_utilities.jl")

function TightBindingInfo(asd)
    TightBindingInfo(tb_info(asd)...)
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
function TightBindingDensity(asd::Dict{String,Any}; style=:normal)
    TightBindingDensity(TightBindingInfo(asd),style = style)
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
