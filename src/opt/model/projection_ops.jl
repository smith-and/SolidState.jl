
struct ProjectionOperators{AR <: AbstractArray, AC <: AbstractArray}
    l_pjs::Array{AC,1}
    sl_pjs::Array{AC,1}
    structure::AR
    vds::Array{AC,1}
    l_pjs_0::Array{AC,1}
    sl_pjs_0::Array{AC,1}
end

function ProjectionOperators(dim, projN, slv_prj_N)
    ProjectionOperators(
        [SharedArray(zeros(ComplexF64,(dim,dim)))   for _=1:projN],
        [SharedArray(zeros(ComplexF64,(dim,dim)))   for _=1:slv_prj_N],
        SharedArray(zeros(Float64,(dim,dim))),
        [SharedArray(zeros(ComplexF64,(dim)))       for _=1:dim],
        [SharedArray(zeros(ComplexF64,(dim,dim)))   for _=1:projN],
        [SharedArray(zeros(ComplexF64,(dim,dim)))   for _=1:slv_prj_N]
    )
end

function ProjectionOperators(dim, pz::Array{Complex{Float64},2},pjs::Vector{Array{Complex{Float64},2}},slv_prjs::Vector{Array{Complex{Float64},2}})
    ProjectionOperators(
        SharedArray(pz),
        SharedArray.(pjs),
        SharedArray.(slv_prjs),
        SharedArray(zeros(Float64,(dim,dim))),
        [SharedArray(zeros(ComplexF64,(dim)))     for _=1:dim],
        SharedArray.(pjs),
        SharedArray.(slv_prjs)
    )
end
