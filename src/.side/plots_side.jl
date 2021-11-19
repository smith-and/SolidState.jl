
##################################################
#### Symmetric Unit Cell
##################################################

const SO3Defining = [ [ 0 0 0 ;
  0 0 -1 ;
  0 1 0 ], [  0 0 1 ;
   0 0 0 ;
  -1 0 0 ], [ 0 -1 0 ;
  1 0 0 ;
  0 0 0 ],
]

axisangle(θ,n) = exp(θ.*sum(n.*SO3Defining))

function group_generate(G0)
        G = deepcopy(G0)
        for g∈G0, h∈G
                if (Ref(h*g) .- G).|>norm.|>(x->x>1e-8)|>all
                        push!(G,h*g)
                end
        end
        if length(G)==length(G0)
                return G0
        else
                group_generate(G)
        end
end

function parallel_project(bt0,direction)
    bt = deepcopy(bt0)
    filter!(bt.vertices) do branch
        bond = branch.a .- branch.b
        abs(acos(dot(direction,bond)/(norm(bond)*norm(direction)))) < 1e-5
    end
    bt
end

function branch_parallel_Q(branch,direction)
    bond = branch.a .- branch.b
    abs(acos(dot(direction,bond)/(norm(bond)*norm(direction)))) < 1e-5
end

function branch_vertical_Q(branch)
    (abs(branch.a[1] -  branch.b[1]) < 1e-5) && (abs(branch.a[2]-branch.b[2]) < 1e-5)
end

function plane_bonds(asd)
    asdg = asd|>SolidState.ASDGeometry
    boundaryQ = SolidState.BoundaryEdgeQ((asdg["uc_c"]...,asdg["uc_c"][1],asdg["uc_c"][2]))
    bt = GraphTree(typeof(EdgeLine(rand(3),rand(3)))[],VertexPoint[],typeof(EdgeLine(rand(3),rand(3)))[],VertexPoint[],zeros(3,3))
    for (i,site1) in enumerate(getindex.(asd["sites"],3))
        for site2 in getindex.(asd["sites"],3)
            if boundaryQ((site1+site2)/2) && (abs(site1[3]-site2[3]) < 1e-5)
                push!(bt.edges,EdgeLine(copy(site1),copy(site2)))
            end
        end
    end
    bt
end

function vertical_bonds(asd)
    asdg = asd|>SolidState.ASDGeometry
    bt = GraphTree(typeof(EdgeLine(rand(3),rand(3)))[],VertexPoint[],typeof(EdgeLine(rand(3),rand(3)))[],VertexPoint[],zeros(3,3))
    for (i,site1) in enumerate(getindex.(asd["sites"],3))
        for site2 in getindex.(asd["sites"],3)
            push!(bt.edges,EdgeLine(copy(site1),copy(site2)))
        end
    end
    btz = parallel_project(bt,[0.0,0.0,1.0])

    btz
end
,
