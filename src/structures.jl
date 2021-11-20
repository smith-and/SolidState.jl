##########################################
#### Reference Systems
##########################################

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

export Control1L0
function Control1L0()
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

export Control2LMzWG
function Control2LMzWG()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 4.0; ϵNN  = -4.0
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


export Control4LMz
function Control4LMz()
    a    = 1;   c    = 1; δ    = .1
    vppπ = 1.0; vppσ = 0.1
    ϵBB  = 4.0; ϵNN  = -4.0
    cut  = 1.001*c
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 4*c],
        "sk" => Dict(
            "Atom" => 1,
            "Layer" => 2,
            "Pos" => 3,
            "Spin" => 4,
            "Orbital" => 5,
            "Glyph" => 6,
            "Color" => 7,
        ),
        "sites" => [
            ("B", 1, [0.0, -a / sqrt(3), -3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, a / sqrt(3), -3c / 2],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, -a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, -a / sqrt(3), 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, a / sqrt(3), 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
        "lbase"=>3, "cutoff" => cut,
        "regulator" => "exp",
        "scales"    => Control_Scales(a,c,δ,vppσ,vppπ,ϵBB,ϵNN),
    )
end


##########################################
#### Boron Nitride System
##########################################

function E_Scales(a, c, δBB, δNN, δNB, ϵBB, ϵNN, vppσ, vppπBB, vppπNN, vppπNB)
    vσσ = 0
    Dict(
        ("N", "N") => (
            as = [0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0; 0 0 0 a/sqrt(3)],
            ds = [δNN δNN δNN δNN; δNN δNN δNN δNN; δNN δNN δNN δNN; δNN δNN δNN δNN],
            vs = [
                vσσ 0 0 0
                0 vppπNN 0 0
                0 0 vppσ 0
                0 0 0 vppπNN
            ],
            es = [
                ϵNN 0 0 0
                0 ϵNN 0 0
                0 0 ϵNN 0
                0 0 0 ϵNN
            ],
        ),
        ("B", "B") => (
            as = [0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0; 0 0 0 a/sqrt(3)],
            ds = [δBB δBB δBB δBB; δBB δBB δBB δBB; δBB δBB δBB δBB; δBB δBB δBB δBB],
            vs = [vσσ 0 0 0; 0 vppπBB 0 0; 0 0 vppσ 0; 0 0 0 vppπBB],
            es = [ϵBB 0 0 0; 0 ϵBB 0 0; 0 0 ϵBB 0; 0 0 0 ϵBB],
        ),
        ("N", "B") => (
            as = [0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0; 0 0 0 a/sqrt(3)],
            ds = [δNB δNB δNB δNB; δNB δNB δNB δNB; δNB δNB δNB δNB; δNB δNB δNB δNB],
            vs = [vσσ 0 0 0; 0 vppπNB 0 0; 0 0 vppσ 0; 0 0 0 vppπNB],
            es = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0],
        ),
        ("B", "N") => (
            as = [0 0 0 0; 0 a/sqrt(3) 0 0; 0 0 c 0; 0 0 0 a/sqrt(3)],
            ds = [δNB δNB δNB δNB; δNB δNB δNB δNB; δNB δNB δNB δNB; δNB δNB δNB δNB],
            vs = [vσσ 0 0 0; 0 vppπNB 0 0; 0 0 vppσ 0; 0 0 0 vppπNB],
            es = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0],
        ),
    )
end

function wide_gap_wha()
    a = 0.2512
    c = 0.334
    δBB = 0.022
    δNN = 0.022
    δNB = 0.0475
    vppπBB = 0.7
    vppπNN = 0.15
    vppπNB = 2.3
    vppσ = 0.2
    ϵBB = 3.125
    ϵNN = -3.125
    Dict(
        "lbase" => 10,
        "cutoff" => 20 * a / sqrt(3),
        "regulator" => "exp",
        "scales" => E_Scales(a, c, δBB, δNN, δNB, ϵBB, ϵNN, vppσ, vppπBB, vppπNN, vppπNB),
    )
end

function small_gap_wha()
    a = 0.2512
    c = 0.334
    δBB = 0.022
    δNN = 0.022
    δNB = 0.0475
    vppπBB = 0.7
    vppπNN = 0.15
    vppπNB = 2.3
    vppσ = 0.2
    ϵBB = 0.10
    ϵNN = -0.10
    Dict(
        "lbase" => 10,
        "cutoff" => 20 * a / sqrt(3),
        "regulator" => "exp",
        "scales" => E_Scales(a, c, δBB, δNN, δNB, ϵBB, ϵNN, vppσ, vppπBB, vppπNN, vppπNB),
    )
end

function random_onsite()
    α = 0.1
    a = 0.2512
    c = 0.334
    δBB = 0.022
    δNN = 0.022
    δNB = 0.0475
    vppπBB = 0.7
    vppπNN = 0.15
    vppπNB = 2.3
    vppσ = 0.2
    ϵBB = 4.0*(1+(rand()-1/2)*α)
    ϵNN = -4.0*(1+(rand()-1/2)*α)
    Dict(
        "lbase" => 10,
        "cutoff" => 20 * a / sqrt(3),
        "regulator" => "exp",
        "scales" => E_Scales(a, c, δBB, δNN, δNB, ϵBB, ϵNN, vppσ, vppπBB, vppπNN, vppπNB),
    )
end

function random_hopping()
    α = 0.1
    a = 0.2512
    c = 0.334
    δBB = 0.022
    δNN = 0.022
    δNB = 0.0475
    vppπBB = 0.7*(1+(rand()-1/2)*α)
    vppπNN = 0.15*(1+(rand()-1/2)*α)
    vppπNB = 2.3*(1+(rand()-1/2)*α)
    vppσ = 0.2*(1+(rand()-1/2)*α)
    ϵBB = 4.0
    ϵNN = -4.0
    Dict(
        "lbase" => 10,
        "cutoff" => 20 * a / sqrt(3),
        "regulator" => "exp",
        "scales" => E_Scales(a, c, δBB, δNN, δNB, ϵBB, ϵNN, vppσ, vppπBB, vppπNN, vppπNB),
    )
end


export ASD1
function ASD1()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 c],
        "sk" => Dict(
            "Atom" => 1,
            "Layer" => 2,
            "Pos" => 3,
            "Spin" => 4,
            "Orbital" => 5,
            "Glyph" => 6,
            "Color" => 7,
        ),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange,),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export ASD2
function ASD2()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 2*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2],(0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, -a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export SGASD2
function SGASD2()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 2*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2],(0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, -a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, small_gap_wha())
end


export ASD2RH
function ASD2RH()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 2*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2],(0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, -a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, random_hopping())
end

export ASD2CS
function ASD2CS()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 2*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2],(0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, -a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export ASD2AB
function ASD2AB()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 2*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2],(0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, 2a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, 0, c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end


export SGASD2AB
function SGASD2AB()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 2*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2],(0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, 2a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, 0, c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, small_gap_wha())
end

export ASD2p3
function ASD2p3()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 c],
        "sk" => Dict(
            "Atom" => 1,
            "Layer" => 2,
            "Pos" => 3,
            "Spin" => 4,
            "Orbital" => 5,
            "Glyph" => 6,
            "Color" => 7,
        ),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [1,2,3]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [1,2,3]), :circle, :orange,),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export ASD3
function ASD3()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 3*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3),  -c], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 1, [0.0, -a / sqrt(3),  0],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, a / sqrt(3),   0],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, a / sqrt(3),  c],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, -a / sqrt(3), c],  (0 // 1, [1]), (1, [2]), :circle, :orange),

        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export ASD3AB
function ASD3AB()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 3*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3),  -c], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 1, [0.0, -a / sqrt(3),  0],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, a / sqrt(3),   0],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, 2a / sqrt(3),  c],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, 0, c],  (0 // 1, [1]), (1, [2]), :circle, :orange),

        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export ASD3RO
function ASD3RO()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 3*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3),  -c], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, -a / sqrt(3),  0],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, a / sqrt(3),   0],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 3, [0.0, a / sqrt(3),  c],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 3, [0.0, -a / sqrt(3), c],  (0 // 1, [1]), (1, [2]), :circle, :orange),

        ],
        "filling" => 0.5,
    )
    merge(asd, random_onsite())
end

export ASD3RH
function ASD3RH()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 3*c],
        "sk" => Dict("Atom" => 1, "Layer" => 2,"Pos" => 3,"Spin" => 4,"Orbital" => 5,"Glyph" => 6,"Color" => 7),
        "sites" => [
            ("B", 1, [0.0, a / sqrt(3),  -c], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, -a / sqrt(3),  0],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, a / sqrt(3),   0],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 3, [0.0, a / sqrt(3),  c],  (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 3, [0.0, -a / sqrt(3), c],  (0 // 1, [1]), (1, [2]), :circle, :orange),

        ],
        "filling" => 0.5,
    )
    merge(asd, random_hopping())
end

export ASD4
function ASD4()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 4*c],
        "sk" => Dict(
            "Atom" => 1,
            "Layer" => 2,
            "Pos" => 3,
            "Spin" => 4,
            "Orbital" => 5,
            "Glyph" => 6,
            "Color" => 7,
        ),
        "sites" => [
            ("B", 1, [0.0, -a / sqrt(3), -3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, a / sqrt(3), -3c / 2],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, -a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, -a / sqrt(3), 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, a / sqrt(3), 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export ASD4AB
function ASD4AB()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 4*c],
        "sk" => Dict(
            "Atom" => 1,
            "Layer" => 2,
            "Pos" => 3,
            "Spin" => 4,
            "Orbital" => 5,
            "Glyph" => 6,
            "Color" => 7,
        ),
        "sites" => [
            ("B", 1, [0.0, -a / sqrt(3), -3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, a / sqrt(3), -3c / 2],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, 2.0a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, 0.0, c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, 0.0, 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, 2.0a / sqrt(3), 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end

export ASD4CS
function ASD4CS()
    a = 0.2512
    c = 0.334
    asd = Dict(
        "blv" => [a 0 0; a/2 sqrt(3)/2*a 0; 0 0 4*c],
        "sk" => Dict(
            "Atom" => 1,
            "Layer" => 2,
            "Pos" => 3,
            "Spin" => 4,
            "Orbital" => 5,
            "Glyph" => 6,
            "Color" => 7,
        ),
        "sites" => [
            ("B", 1, [0.0, -a / sqrt(3), -3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, a / sqrt(3), -3c / 2],  (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 1, [0.0, a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 1, [0.0, -a / sqrt(3), -c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, -a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, a / sqrt(3), c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
            ("B", 2, [0.0, a / sqrt(3), 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :blue),
            ("N", 2, [0.0, -a / sqrt(3), 3c / 2], (0 // 1, [1]), (1, [2]), :circle, :orange),
        ],
        "filling" => 0.5,
    )
    merge(asd, wide_gap_wha())
end


##########################################
#### Graphene Systems
##########################################

export G1
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

function G1()
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