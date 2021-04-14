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

#=
function WHA()
    a = 0.2512
    δBB = 0.022
    δNN = 0.022
    δNB = 0.022
    vppπBB = 0.7
    vppπNN = 0.15
    vppπNB = 2.3
    vppσ = 0.2
    ϵBB = 4.0
    ϵNN = -4.0
    c = 0.333
    a = 0.2512
    Dict(
        "lbase" => 10,
        "cutoff" => 4.001 * a / sqrt(3),
        "regulator" => "exp",
        "scales" => E_Scales(a, c, δBB, δNN, δNB, ϵBB, ϵNN, vppσ, vppπBB, vppπNN, vppπNB),
    )
end
=#

function WHA()
    a = 0.2512
    c = 0.333
    δBB = 0.022
    δNN = 0.022
    δNB = 0.0475
    vppπBB = 0.7
    vppπNN = 0.15
    vppπNB = 2.3
    vppσ = 0.2
    ϵBB=4.0
    ϵNN=-4.0
    Dict(
        "lbase" => 10,
        "cutoff" => 4.001 * a / sqrt(3),
        "regulator" => "exp",
        "scales" => E_Scales(a, c, δBB, δNN, δNB, ϵBB, ϵNN, vppσ, vppπBB, vppπNN, vppπNB),
    )
end

export ASD
function ASD()
    a = 0.2512
    c = 0.333
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
    merge(asd, WHA())
end

export ASD2
function ASD2()
    a = 0.2512
    c = 0.333
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
    merge(asd, WHA())
end

export ASD2B
function ASD2B()
    a = 0.2512
    c = 0.333
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
    merge(asd, WHA())
end

export ASD2p3
function ASD2p3()
    a = 0.2512
    c = 0.333
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
    merge(asd, WHA())
end

export ASD4
function ASD4()
    a = 0.2512
    c = 0.333
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
    merge(asd, WHA())
end
