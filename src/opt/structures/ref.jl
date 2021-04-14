function Control_Scales(a0,c,δ,vppσ,vppπ,ϵBB,ϵNN)
    a= a0/sqrt(3) ; ds = δ*ones(4,4)
    as=[0 0 0 0;
        0 a 0 0;
        0 0 c 0;
        0 0 0 a ]
    vs=[0    0 0 0;
        0 vppπ 0 0;
        0 0 vppσ 0;
        0 0 0 vppπ]
    Dict(
        ("N","N")=>(as=as, ds=ds, vs=vs,
        es=[0 0 0 0; 0 ϵNN 0 0 ; 0 0 ϵNN  0; 0 0 0 ϵNN]),
        ("B","B")=>(as=as, ds=ds, vs=vs,
        es=[0 0 0 0; 0 ϵBB 0 0 ; 0 0 ϵBB  0; 0 0 0 ϵBB]),
        ("N","B")=>(as=as, ds=ds, vs=vs,
        es=[0 0 0 0; 0 0 0 0   ; 0 0 0    0; 0 0 0 0]  ),
        ("B","N")=>(as=as, ds=ds, vs=vs,
        es=[0 0 0 0; 0 0 0 0   ; 0 0 0    0; 0 0 0 0]  )
    )
end

export Control1L
function Control1L()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 0.5; ϵNN  = -0.5
    cut  = 1.001*a/sqrt(3)
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0,  a/sqrt(3),  0.0],   (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 1, [0.0, -a/sqrt(3), 0.0],   (0//1,[1]), (1,[2]), :circle, :orange)
                        ],
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end

export Control1LI
function Control1LI()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 0.0; ϵNN  = 0.0
    cut  = 1.001*a/sqrt(3)
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0, a/sqrt(3),  0.0],   (0//1,[1]), (1,[2]), :circle, :black),
                        ("N", 1, [0.0, -a/sqrt(3), 0.0],   (0//1,[1]), (1,[2]), :circle, :black)
                        ],
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end

export Control2LMz
function Control2LMz()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 0.5; ϵNN  = -0.5
    cut  = 1.001*a
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 2*c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :orange),
                        ("B", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle, :orange)
                        ],
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end

export Control2LMzZero
function Control2LMzZero()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.0
    ϵBB  = 0.5; ϵNN  = -0.5
    cut  = 1.001*a
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 2*c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :orange),
                        ("B", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle, :orange)
                        ],
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end

export Control2LI
function Control2LI()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 0.5; ϵNN  = -0.5
    cut  = 1.001*a
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 2*c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :orange),
                        ("B", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle, :orange)
                        ],
        "lbase"=>5, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end

export Control2LMzBernal
function Control2LMzBernal()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 0.5; ϵNN  = -0.5
    cut  = 1.001*c
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 2*c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 1, [0.0, 0.0,        -c/2], (0//1,[1]), (1,[2]), :circle, :orange),
                        ("B", 2, [0.0, 0.0,         c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle, :orange)
                        ],
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end

export Control2LIBernal
function Control2LIBernal()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 0.5; ϵNN  = -0.5
    cut  = 1.001*c
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 2*c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 1, [0.0, 0.0,        -c/2], (0//1,[1]), (1,[2]), :circle, :orange),
                        ("B", 2, [0.0, -a/sqrt(3),         c/2], (0//1,[1]), (1,[2]), :circle,  :blue),
                        ("N", 2, [0.0, 0.0,  c/2], (0//1,[1]), (1,[2]), :circle, :orange)
                        ],
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end

export Control2LSM
function Control2LSM()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 0.0; ϵNN  = 0.0
    cut  = 1.001*c
    Dict(
        "blv"       => [a 0 0 ; a/2 (sqrt(3)/2*a) 0 ; 0 0 2*c ] ,
        "sk"        => Dict("Atom"=>1,"Layer"=>2,"Pos"=>3,"Spin"=>4,"Orbital"=>5,"Glyph"=>6,"Color"=>7),
        "sites"     => [
                        ("B", 1, [0.0, a/sqrt(3),  -c/2], (0//1,[1]), (1,[2]), :circle, :black),
                        ("N", 1, [0.0, -a/sqrt(3), -c/2], (0//1,[1]), (1,[2]), :circle, :black),
                        ("B", 2, [0.0, a/sqrt(3),   c/2], (0//1,[1]), (1,[2]), :circle, :black),
                        ("N", 2, [0.0, -a/sqrt(3),  c/2], (0//1,[1]), (1,[2]), :circle, :black)
                        ],
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
        "filling"   => 0.5
    )
end
