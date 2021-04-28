
function hs_spectra_v_twist(asd0,series,ksymbols; force = false, cachedir, outdir, kargs...)
    datafile = "$outdir/hs_spectra_v_twist-$(series[1])-$(series[2])-$(series[3]).bson"
    if !isfile(datafile) || force
        #Check for model information and make it if absent
        comargs = SolidState.twist_series(series...);
        SolidState.make_models(asd0,comargs;cachedir=cachedir)

        #Initialize data structures for plots and spectra
        kenergies = copy(OrderedDict("hi"=>OrderedDict(1.0=>rand(5))))|>empty
        for k ∈ ksymbols
            push!(kenergies,k=>(OrderedDict(1.0=>rand(5))|>empty))
        end

        #Loop over different twisted models
        for mn ∈ comargs
            println("Loading $(mn)");flush(stdout)
            asd = BSON.load("$cachedir/$asd0/asd-$(mn[1])-$(mn[2]).bson")
            asdg = asd|>SolidState.ASDGeometry
            kledger  = Dict(ksymbols|>x->x.=>getindex.(Ref(asdg["bz_hs"]),x))

            hd =  data_import("$cachedir/$asd0/hd-$(mn[1])-$(mn[2]).bson")
            # Looping over different k points
            for k ∈ ksymbols
                key = SolidState.cθ(mn...)*180/π
                println("spectrum calculated at $k");flush(stdout)
                hd(kledger[k])
                push!(kenergies[k],key=>eigvals(Hermitian(hd.h_ops.h)))
            end

        end


        bson(datafile,
            Dict(
                :spectra  =>kenergies,
                :asd      =>asd0,
                :comargs  =>comargs,
                :ksymbols =>ksymbols
                )
            )
    end
    return datafile
end

function com_windows(comargs)
    angles  = Vector{Float64}(undef,length(comargs))
    for (i,mn) ∈ enumerate(comargs)
        angles[i] = SolidState.cθ(mn...)*180/π
    end

    comargs .= comargs[sortperm(angles)];
    angles .= angles[sortperm(angles)];

    leftwindows = angles .+  vcat([0.0],reverse(diff(reverse(angles))))./2
    rightwindows = angles .+ vcat(diff(angles),[0.0])./2

    windows = eachslice(hcat(leftwindows,rightwindows),dims=1)|>collect

    OrderedDict((angles.=>windows)...)
end

function com_windows_loose(comargs;α=0.1)
    angles  = Vector{Float64}(undef,length(comargs))
    for (i,mn) ∈ enumerate(comargs)
        angles[i] = SolidState.cθ(mn...)*180/π
    end

    comargs .= comargs[sortperm(angles)];
    angles .= angles[sortperm(angles)];

    leftwindows = (angles  .+ vcat([0.0],reverse(diff(reverse(angles))))./2)*(1+α)
    rightwindows = angles .+ vcat(diff(angles),[0.0])./2*(1-α)

    windows = eachslice(hcat(leftwindows,rightwindows),dims=1)|>collect

    OrderedDict((angles.=>windows)...)
end

function plot_twist_spectra_levels(datafile; α,  outdir, saving = true, kargs...)
    data      = BSON.load(datafile)
    ksymbols  = data[:ksymbols]
    kenergies = data[:spectra]
    comargs   = data[:comargs]
    windows   = com_windows_loose(comargs; α=α)

    #Initialize data structures for plots and spectra
    kplots = copy(OrderedDict("hi"=>plot()))|>empty
    for k ∈ ["Γ","K1"]
        push!(kplots,
            k=>plot(;
                label="",
                frame = :box,
                title = "$k twist spectrum",
                xguide = "twist angle (deg)",
                yguide = "energy (eV)",
                color=:black,
            )
        )
        for key ∈ kenergies[k]|>keys
            println("Plotting $k $key");flush(stdout)
            plot!(kplots[k],
                fill(windows[key],length(kenergies[k])),
                collect(eachslice(hcat(kenergies[k][key],kenergies[k][key]),dims=1));
                label="",
                color=:black,
                kargs...
            )
            display(kplots[k])
        end
    end
    saving && println("Saving Plots");flush(stdout)
    saving && Plots.pdf.(values(kplots),"$outdir/".*keys(kplots))
    saving && Plots.png.(values(kplots),"$outdir/".*keys(kplots))

    kplots|>values
end

function plot_twist_spectra_points(datafile; α,  outdir, saving = true, kargs...)
    data      = BSON.load(datafile)
    ksymbols  = data[:ksymbols]
    kenergies = data[:spectra]
    comargs   = data[:comargs]
    windows   = com_windows_loose(comargs; α=α)

    #Initialize data structures for plots and spectra
    kplots = copy(OrderedDict("hi"=>plot()))|>empty
    for k ∈ ["Γ","K1"]
        push!(kplots,
            k=>plot(;
                label="",
                frame = :box,
                title = "$k twist spectrum",
                xguide = "twist angle (deg)",
                yguide = "energy (eV)",
                color=:black,
            )
        )
        for key ∈ kenergies[k]|>keys
            scatter!(kplots[k],
            fill(key,length(kenergies[k][key])),kenergies[k][key];
            label="",
            color=:black,
            marker=:circle,
            kargs...
            )
        end
        display(kplots[k])
    end
    saving && println("Saving Plots");flush(stdout)
    saving && Plots.pdf.(values(kplots),"$outdir/".*keys(kplots))
    saving && Plots.png.(values(kplots),"$outdir/".*keys(kplots))

    kplots|>values
end

function spectral_density(energies,damping,kenergies)
    d = length(energies)
    density = Array{Float64,2}(undef,(length(energies),length(keys(kenergies))))
    for (i,en) in enumerate(energies)
        for (j,kens) in enumerate(collect(values(kenergies)))
            z= 0.0
            for ken in kens
                z += imag(1.0 / Complex(ken-en,damping))
            end
            density[i+d*(j-1)] = -z
        end
    end
    return density
end

function plot_twist_spectral_density(datafile; cachedir, damping = 0.02, N=10, outdir, saving = true, kargs...)
    data      = BSON.load(datafile)
    ksymbols  = data[:ksymbols]
    kenergies = data[:spectra]
    comargs   = data[:comargs]
    windows   = com_windows(comargs)

    ucvols = Vector{Float64}(undef,length(comargs))
    for (i,mn)∈ enumerate(comargs)
        ucvols[i] = det(BSON.load("$cachedir/$(data[:asd])/asd-$(mn[1])-$(mn[2]).bson")["blv"][1:2,1:2])
    end

    #Initialize data structures for plots and spectra
    kplots = copy(OrderedDict("hi"=>plot()))|>empty
    for k ∈ ["Γ","K1"]

        emin = min([min(ens...) for ens in values(kenergies[k])]...)
        emax = max([max(ens...) for ens in values(kenergies[k])]...)

        energies = emin:((emax-emin)/N):emax
        angles = keys(kenergies[k])|>collect
        density = spectral_density(energies,damping,kenergies[k])
        spa     = sortperm(angles)
        vols    = ucvols[spa]


        println("Plotting $k");flush(stdout)
        kplots[k] = heatmap(angles[spa],energies,(density./(ucvols'))[:,spa];
            title = "$k twist spectrum",
            xguide = "twist angle (deg)",
            yguide = "energy (eV)",
            kargs...
            )
        display(kplots[k])
    end

    saving && println("Saving Plots");flush(stdout)
    saving && Plots.pdf.(values(kplots),"$outdir/".*keys(kplots))
    saving && Plots.png.(values(kplots),"$outdir/".*keys(kplots))

    kplots|>values
end
