##Carbon Graphene Systems
export G
export GAA
export GBernal

function Graphene_Scales(a,c,δ,vppσ,vppπ,ϵCC)
    Dict(
        ("C","C")=>(
            as=[0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0;0 0 0 a/sqrt(3)],
            ds=[δ δ δ δ; δ δ δ δ; δ δ δ δ; δ δ δ δ],
            vs=[0 0 0 0; 0 vppπ 0 0; 0 0 vppσ 0; 0 0 0 vppπ],
            es=[0 0 0 0; 0 ϵCC 0 0 ; 0 0 ϵCC  0; 0 0 0 ϵCC])
    )
end

function G()
    a=0.246
    c=0.335
    δ= 0.184*a
    vppπ=-2.7
    vppσ=0.48
    ϵCC=0.0
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"=>[
            ("C", 1, [0.0, a/sqrt(3),  0.0], (0//1,[1]), (1,[2]), :circle, :black),
            ("C", 1, [0.0, -a/sqrt(3), 0.0], (0//1,[1]), (1,[2]), :circle, :black),
        ],
        "lbase"=>7,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>Graphene_Scales(a,c,δ,vppσ,vppπ,ϵCC),
        "filling"   => 0.5
    )
end

export G3p
function G3p()
    a=0.246
    c=0.335
    δ= 0.184*a
    vppπ=-2.7
    vppσ=0.48
    ϵCC=0.0
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"=>[
            ("C", 1, [0.0, a/sqrt(3),  0.0], (0//1,[1]), (1,[1,2,3]), :circle, :black),
            ("C", 1, [0.0, -a/sqrt(3), 0.0], (0//1,[1]), (1,[1,2,3]), :circle, :black),
        ],
        "lbase"=>7,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>Graphene_Scales(a,c,δ,vppσ,vppπ,ϵCC),
        "filling"   => 0.5
    )
end

function GAA()
    a=0.246
    c=0.335
    δ= 0.184*a
    vppπ=-2.7
    vppσ=0.48
    ϵCC=0
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("C", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ("C", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ("C", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ("C", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ],
        "lbase"=>7,"cutoff"=>4.00001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>Graphene_Scales(a,c,δ,vppσ,vppπ,ϵCC),
        "filling"   => 0.5
    )
end

function GBernal()
    a=0.246
    c=0.335
    δ= 0.184*a
    vppπ=-2.7
    vppσ=.48
    ϵCC=0
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("C", 1, [0.0, 0.0,        -c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ("C", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ("C", 2, [0.0, 0.0,         c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ("C", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle, :black),
            ],
        "lbase"=>7,"cutoff"=>4.00001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>Graphene_Scales(a,c,δ,vppσ,vppπ,ϵCC),
        "filling"   => 0.5
    )
end
