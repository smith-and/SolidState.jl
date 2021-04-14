"""
   TensorChart{IndexT <: AbstractArray, PriorT <: AbstractArray, BaseT <: AbstractArray, DataT <: AbstractArray}

A generic structure to hold Tensor Data
"""
struct TensorChart{IndexT <: AbstractArray, PriorT <: AbstractArray, BaseT <: AbstractArray, DataT <: AbstractArray}
   indices::IndexT
   priors::PriorT
   base::BaseT
   data::DataT
   l_i::Int
   l_p::Int
   l_b::Int
end

"""
   TensorChart(Indices::AbstractArray, Priors::AbstractArray, Base::AbstractArray, DType::DataType; style=:normal)
"""
function TensorChart(Indices::AbstractArray, Priors::AbstractArray, Base::AbstractArray, DType::DataType; style=:normal)
   if      style == :static
      TensorChart(Indices, Priors, Base, @MArray(zeros(DType,(size(Indices)...,size(Priors)...,size(Base)...)...)) ,length(Indices),length(Priors),length(Base))
   elseif  style ==:shared
      TensorChart(SharedArray(Indices), SharedArray(Priors), SharedArray(Base), SharedArray(zeros(DType,(size(Indices)...,size(Priors)...,size(Base)...))),length(Indices),length(Priors),length(Base))
   else
      TensorChart(Indices, Priors, Base, (zeros(DType,(size(Indices)...,size(Priors)...,size(Base)...)...)),length(Indices),length(Priors),length(Base))
   end
end

### Application Specific Wrappings for TensorCharts
"""
   abstract type DataChart end
"""
export DataChart
abstract type DataChart end

"""
   abstract type SpectralChart <: DataChart end
"""
abstract type SpectralChart <: DataChart end
include.(readdir(string(@__DIR__)*"/../DataCharts/SpectralCharts",join=true))

"""
   abstract type ReponseChart  <: DataChart end
"""
abstract type ReponseChart  <: DataChart end
include.(readdir(string(@__DIR__)*"/../DataCharts/ResponseCharts",join=true))

"""
const RangeType = Union{
   Tuple{Symbol,NTuple{N,T}} where {N, T <: AbstractRange}, Tuple{Symbol,AV} where AV <: AbstractVector,
   Vector{Tuple{Symbol,N,N,I}},
   Tuple{Symbol,Symbol,I},
   Vector{NTuple{NT,I}}} where {N <: Number, I <: Integer, NT
}

"""
const RangeType = Union{
         Tuple{Symbol,NTuple{N,T}} where {N, T <: AbstractRange}, Tuple{Symbol,AV} where AV <: AbstractVector,
         Vector{Tuple{Symbol,N,N,I}},
         Tuple{Symbol,Symbol,I},
         Vector{NTuple{NT,I}}} where {N <: Number, I <: Integer, NT
      }

"""
struct ChartInfo{T <: DataChart, RT1 <: RangeType, RT2 <: RangeType, RT3 <: RangeType, CAT <: Union{Tuple{Int,Int}}}
    type::Type{T}
    indices::RT1
    priors::RT2
    base::RT3
    mn::CAT
end
"""
struct ChartInfo{T <: DataChart, RT1 <: RangeType, RT2 <: RangeType, RT3 <: RangeType, CAT <: Union{Tuple{Int,Int}}}
    type::Type{T}
    indices::RT1
    priors::RT2
    base::RT3
    mn::CAT
end

"""
function ChartInfo()
   ChartInfo(DataChart,[(1,)],[(1,)],[(1,)],(0,0))
end
"""
function ChartInfo()
   ChartInfo(DataChart,[(1,)],[(1,)],[(1,)],(0,0))
end
