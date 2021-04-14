##Boron Nitride Systems
export BN
export BNAA
export BNAB
export BNAC

function BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN)
    Dict(
        ("N","N")=>(
            as=[0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0;0 0 0 a/sqrt(3)],
            ds=[δ δ δ δ; δ δ δ δ; δ δ δ δ; δ δ δ δ],
            vs=[0 0 0 0; 0 vppπ 0 0; 0 0 vppσ 0; 0 0 0 vppπ],
            es=[0 0 0 0; 0 ϵNN 0 0 ; 0 0 ϵNN  0; 0 0 0 ϵNN]),
        ("B","B")=>(
            as=[0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0;0 0 0 a/sqrt(3)],
            ds=[δ δ δ δ; δ δ δ δ; δ δ δ δ; δ δ δ δ],
            vs=[0 0 0 0; 0 vppπ 0 0; 0 0 vppσ 0; 0 0 0 vppπ],
            es=[0 0 0 0; 0 ϵBB 0 0 ; 0 0 ϵBB  0; 0 0 0 ϵBB]),
        ("N","B")=>(
            as=[0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0;0 0 0 a/sqrt(3)],
            ds=[δ δ δ δ; δ δ δ δ; δ δ δ δ; δ δ δ δ],
            vs=[0 0 0 0; 0 vppπ 0 0; 0 0 vppσ 0; 0 0 0 vppπ],
            es=[0 0 0 0; 0 0 0 0   ; 0 0 0  0  ; 0 0 0 0]),
        ("B","N")=>(
            as=[0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0;0 0 0 a/sqrt(3)],
            ds=[δ δ δ δ; δ δ δ δ; δ δ δ δ; δ δ δ δ],
            vs=[0 0 0 0; 0 vppπ 0 0; 0 0 vppσ 0; 0 0 0 vppπ],
            es=[0 0 0 0; 0 0 0 0   ; 0 0 0  0  ; 0 0 0 0])
    )
end

function BN()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.37
    vppσ=0.6
    ϵBB=3.4/2
    ϵNN=-5.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"=>[
            ("B", 1, [0.0, a/sqrt(3), 0.0],  (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 1, [0.0, -a/sqrt(3), 0.0], (0//1,[1]), (1,[2]), :circle, :orange)
        ],
        "lbase"=>3,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"=>0.5
    )
end

function BNAB()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.37
    vppσ=0.6
    ϵBB=3.4/2
    ϵNN=-5.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :orange),
            ("B", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle, :orange),
            ],
        "lbase"=>3,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"=>0.5
    )

end

export BNAB2p3
function BNAB2p3()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.7
    vppσ=0.3
    ϵBB=4.4/2
    ϵNN=-4.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[1,2,3]), :circle,  :blue),
            ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[1,2,3]), :circle, :orange),
            ("B", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[1,2,3]), :circle,  :blue),
            ("N", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[1,2,3]), :circle, :orange),
            ],
        "lbase"=>6,"cutoff"=>10.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"=>0.5/3
    )
end

export BNABZero
function BNABZero()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.37
    vppσ=0.6
    ϵBB=3.4/2
    ϵNN=-5.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("B", 1, [0.0, a/sqrt(3),  -100*c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 1, [0.0, -a/sqrt(3), -100*c/2], (0//1,[1]), (1,[2]), :circle, :orange),
            ("B", 2, [0.0, a/sqrt(3),   100*c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 2, [0.0, -a/sqrt(3),  100*c/2], (0//1,[1]), (1,[2]), :circle, :orange),
            ],
        "lbase"=>3,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"=>0.5
    )
end

export BNABZero2p3
function BNABZero2p3()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.37
    vppσ=0.6
    ϵBB=3.4/2
    ϵNN=-5.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("B1", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[1,2,3]), :circle,  :blue),
            ("N1", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[1,2,3]), :circle, :orange),
            ("B2", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[1,2,3]), :circle,  :blue),
            ("N2", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[1,2,3]), :circle, :orange),
            ],
        "lbase"=>3,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales_2L(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"=>0.5
    )
end


function BNAA()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.37
    vppσ=0.6
    ϵBB=3.4/2
    ϵNN=-5.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :orange),
            ("B", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle, :orange),
            ],
        "lbase"=>3,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"=>0.5
    )
end

function BNAC()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.37
    vppσ=0.6
    ϵBB=3.4/2
    ϵNN=-5.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 2*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :diamond, :orange),
            ("B", 2, [0.0, -a/sqrt(3)-.1*a,  c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
            ("N", 2, [0.0, a/sqrt(3)-.1*a,   c/2], (0//1,[1]), (1,[2]), :diamond, :orange),
            ],
        "lbase"=>3,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN)
    )
end


#2+2 Layering
export BNABBA
function BNABBA()
    a=.2512
    c=.323
    δ= 0.184*a
    vppπ=-2.37
    vppσ=0.6
    ϵBB=3.4/2
    ϵNN=-5.4/2
    Dict(
        "blv"=>[a 0 0 ; a/2 sqrt(3)/2*a 0 ; 0 0 4*c ] ,
        "sk"=>Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites" => [
            ("B", 1, [0.0, -a/sqrt(3),  -3c/2], (0//1,[1]), (1,[2]), :circle, :blue),
            ("N", 1, [0.0, a/sqrt(3), -3c/2], (0//1,[1]), (1,[2]), :circle, :orange),
            ("B", 1, [0.0, a/sqrt(3), -c/2],  (0//1,[1]), (1,[2]), :circle, :blue),
            ("N", 1, [0.0, -a/sqrt(3),  -c/2],  (0//1,[1]), (1,[2]), :circle, :orange),
            ("B", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle, :blue),
            ("N", 2, [0.0, -a/sqrt(3),  c/2],   (0//1,[1]), (1,[2]), :circle, :orange),
            ("B", 2, [0.0, -a/sqrt(3),  3c/2], (0//1,[1]), (1,[2]), :circle, :blue),
            ("N", 2, [0.0, a/sqrt(3),  3c/2],  (0//1,[1]), (1,[2]), :circle, :orange),
            ],
        "lbase"=>4,"cutoff"=>4.001*a/sqrt(3),
        "regulator"=>"exp",
        "scales"=>BN_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"=>0.5
    )

end
