# Code for calling Hop
using Distributed
@everywhere using BSON, OrderedCollections, LinearAlgebra, Plots
@everywhere using Hop, SolidState, SolidStateApps

function getBN()
    lat = [1 1/2 0; 0 √3/2 0; 0 0 2.0];
    site_positions = lat*([1/3 1/3 0; 2/3 2/3 0; 1/3 1/3 1; 2/3 2/3 1]')
    tm = TBModel(lat, site_positions, [[0], [0], [0], [0]], isorthogonal=true)
    # Onsite Terms
    addhopping!(tm, [0, 0, 0], (1, 1), (1, 1), -0.5)
    addhopping!(tm, [0, 0, 0], (2, 1), (2, 1), 0.5)
    addhopping!(tm, [0, 0, 0], (3, 1), (3, 1), -0.5)
    addhopping!(tm, [0, 0, 0], (4, 1), (4, 1), 0.5)

    # In-plane Neatest Neighbor
    addhopping!(tm, [0, 0, 0], (1, 1), (2, 1), 1.0)
    addhopping!(tm, [-1, 0, 0], (1, 1), (2, 1), 1.0)
    addhopping!(tm, [0, -1, 0], (1, 1), (2, 1), 1.0)

    addhopping!(tm, [0, 0, 0], (3, 1), (4, 1), 1.0)
    addhopping!(tm, [-1, 0, 0], (3, 1), (4, 1), 1.0)
    addhopping!(tm, [0, -1, 0], (3, 1), (4, 1), 1.0)

    # Out-of-Plane Neatest Neighbor
    addhopping!(tm, [0, 0, 1], (1, 1), (3, 1), 0.1)
    addhopping!(tm, [0, 0, 1], (2, 1), (4, 1), 0.1)

    return tm
end


# Calculate and Plot the band structures
function make_data(; asd, paths, Npath, Nbz, dtype, indices, priors, base, dom, kargs...)
    hd = SolidState.TightBindingDensity(asd())
    band_dict = SolidStateApps.band_plot(asd(),hd,paths,Npath)

    ds = DataSection(DataMap(dtype,asd(),indices,priors,base),dom);
    ds();
    dsi = SolidState.integrate(ds);

    di = DataIntegral(DataMap(dtype,asd(),indices,priors,base))
    di(Nbz*Nbz)

    Dict(:di=>di,:ds=>ds,:dsi=>dsi,:bands=>band_dict)
end

function make_data(args)
    make_data(; args...)
end

argf(asd, dtype, Nbz) = Dict(
    # Material Inputs
    :asd => asd,
    # Band Inputs
    :paths => ["K'1","Γ","K1","M1","K'1"],
    :Npath => 100,
    :Nbz   => Nbz,
    #Data Inputs
    :dtype => dtype,
    :indices => [(2,2,2)],
    :priors => [(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.1,0.1,1)],
    :base => [(:ω,0.0,4.0,200)],
    # base => (:ω,:bandwidth,100),
    :dom => [(:n1,-1.0,1.0,Nbz),(:n2,-1.0,1.0,Nbz)],
    :f  => x->abs(x),
)

function compute_set(inputs...)
    args = argf.(inputs...)
    datum = make_data.(args)

    (args,datum)
end
