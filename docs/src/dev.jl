# # Reference Check for Second Harmonic Generation
#=
 We need to compare the code code for (Chongs)[] with (Andrew)[]
 in order to confirm that they produce the same result, given the same model.

 In order to demonstrate equivalence of the code base we compare the results
 each calcuates for a variety of reference models.

 The first such model is


=#
# ## Compute the shg signal in Hop.jl
include("app.jl")

# ## Set up Model
model = Dict(
    :tm  => getBN(),
    :ϵ   => 0.1,
    :μ   => 0.0,
    :Nbz => 100,
    :bs  => 100,
)

# ## Define some domains
doms = OrderedDict(
    :bs   => [ 0 0.5 0 ; 0 0 0 ; 0.5 0 0 ],
    :dos  => -4.0:0.04:4.0|>collect,
    :jdos => 0.0:0.04:8.0|>collect,
    :shg  => 0.0:0.03:4.0|>collect,
)

# ## Calculate Model Properties
datum = OrderedDict(
    :bs   => Hop.BandStructure.getbs(  model[:tm], doms[:bs]', model[:bs], connect_end_points=true),
    :dos  => Hop.BandStructure.getdos( model[:tm], doms[:dos], [model[:Nbz],model[:Nbz],1], ϵ = model[:ϵ]),
    :jdos => Hop.BandStructure.getjdos(model[:tm], doms[:jdos], model[:μ], [model[:Nbz],model[:Nbz],1], ϵ = model[:ϵ]),
    :shg  => Hop.Optics.get_shg(model[:tm], 2, 2, 2, doms[:shg],model[:μ], [model[:Nbz],model[:Nbz],1]; ϵ = model[:ϵ]),
)


# # Compute the signal in SolidState

datum_1 = compute_set(Control2LMz,[SHG],100);
base = getindex.(datum_1[2][1][:di].dm.chart.shg.base,1)
data_dsi = datum_1[2][1][:dsi].data[1][:]
data_di  = datum_1[2][1][:di].data[1][:]

# ## Compare the
α1 = abs.(datum[:shg])[1]/abs.(data_dsi)[1];
α2 = abs.(datum[:shg])[1]/abs.(data_di)[1];
plt = plot(
    frame = :box,
)
plot!(plt, base, abs.(data_dsi)*α1,
    label = "SolidState Mesh"
    )
plot!(plt, base, abs.(data_di)*α2,
    label = "SolidState Adaptive"
    )
plot!(plt, doms[:shg],abs.(datum[:shg]),
    label = "Hop.jl"
    )
