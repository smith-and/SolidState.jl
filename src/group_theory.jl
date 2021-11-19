
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

function C3VectorRep()
    e = axisangle(0,[0,0,1]).|>round;
    C3=axisangle(-2π/3,[0,0,1]);
    G0 = [e,C3];
    G=group_generate(G0);
end

function D3VectorRep()
    e = axisangle(0,[0,0,1]).|>round;
    C3=axisangle(-2π/3,[0,0,1]);
    C2=axisangle(π,[0,-1,0]).|>round;
    G0 = [e,C3,C2];
    G=group_generate(G0);
end


function site_representation(asd,G1)
    Ω = asd["blv"]|>transpose
    uc = getindex.(asd["sites"],3)
    Gτ = [zeros(length(uc),length(uc)) for _=1:length(G1)]
    for (i,g) ∈ enumerate(G1)
            for (j,τ1) ∈ enumerate(uc)
                    for (k,τ2) ∈ enumerate(uc)
                            #Permutation Entry
                            if norm(τ1-g*τ2) < 1e-4
                                    Gτ[i][k,j] = 1.0
                            #Boundary State
                            elseif inv(Ω)*(τ1-g*τ2)|> (x->(norm(round.(x).-x) < 1e-4))
                                    Gτ[i][k,j]  = 1.0
                            #Nothing
                            else
                                    Gτ[i][k,j] = 0.0
                            end
                    end
            end
    end
    Gτ
end

########################################################
###### 3D Symmetrized Unit Cell Plots

function symmetrize_site!(asd,G)
    asd["sites"] = vcat(map(asd["sites"]) do site
        map(G) do g
                newsite = deepcopy(site)
                newsite[3] .= g*site[3]
                newsite
        end
    end...)|> unique

    asd["sites"] = unique(x->round.(x[3],digits=4),asd["sites"])
    nothing
end
