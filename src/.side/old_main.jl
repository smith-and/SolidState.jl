
######################################
#### Calculating High Symmetry Spectra
######################################


function hs_spectra_v_twist(RN, asd0,series,ksymbols, cachedir=ENV["cachedir"], scriptdir=ENV["scriptdir"])
    datafile = "$scriptdir/out/$RN/$asd0-$(series[2])-$(series[3]).bson"

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

    return datafile
end

"""
    bands(asd::Symbol, mn::Tuple{Int,Int}, RN, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...)

Calculate the band structure of a model along the high symmetry points listed in
"""
function twistedbands(RN, asd, mn::Tuple{Int,Int}, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...)
    #Compute Data
    models(asd,[mn])
    ASD       = BSON.load(   "$cachedirr/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd        = data_import( "$cachedirr/$asd/hd-$(mn[1])-$(mn[2]).bson")
    dict = band_data(ASD,hd,pathlist,Npath; title = "$(round(180/π*SolidState.cθ(mn...),digits=3))ᵒ")

    #Export Data
    bson("$(mkpath("$scriptdir/out/$RN/"))/$asd-$(mn[1])-$(mn[2]).bson", deepcopy(dict))

    "$scriptdir/out/$RN/$asd-$(mn[1])-$(mn[2]).bson"
end

"""
    bands(asd::Symbol, mn::Tuple{Int,Int}, RN, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...)

Calculate the band structure of a model along the high symmetry points listed in
"""
function shiftedbands(RN, asd, (n1,n2,N)::Tuple{Int,Int,Int}, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...)
    #Compute Data
    models(asd,[mn])
    ASD       = BSON.load(   "$cachedirr/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd        = data_import( "$cachedirr/$asd/hd-$(mn[1])-$(mn[2]).bson")
    dict = band_data(ASD,hd,pathlist,Npath; title = "$(round(180/π*SolidState.cθ(mn...),digits=3))ᵒ")

    #Export Data
    bson("$(mkpath("$scriptdir/out/$RN/"))/$asd-$(mn[1])-$(mn[2]).bson", deepcopy(dict))

    "$scriptdir/out/$RN/$asd-$(mn[1])-$(mn[2]).bson"
end



function shg_section(;RN,asd,mn,Nbz)
    models(asd,[mn])

    asdmn = BSON.load("$(ENV["scriptdir"])/.cache/$asd/asd-$(mn[1])-$(mn[2]).bson")
    asdg = (asdmn|>SolidState.ASDGeometry)

    hd = data_import("$(ENV["scriptdir"])/.cache/$asd/hd-$(mn[1])-$(mn[2]).bson")

    K = KinematicDensity(hd,[(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)])

    dim = size(K.hd.h_ops.h,1)
    vband = dim/2|>Int
    cband = vband+1

    patch_rng = range(-1.25*asdg["Kmag"],1.25*asdg["Kmag"], length=Nbz)

    data = Dict(
        :patch_rng => patch_rng,
        :asdg => asdg,
        :asd => asd,
        :mn  => mn,
        :T1d  => zeros(Complex{Float64},(Nbz,Nbz)),
        :T2d  => zeros(Complex{Float64},(Nbz,Nbz)),
    );

    for (i,kx) in enumerate(patch_rng)
        for (j,ky) in enumerate(patch_rng)
            ij = i + (j-1)*Nbz
            K([patch_rng[i],patch_rng[j]])
            data[:T1d][ij] = -im/(2.0*K.k_m.dω[cband,vband]^2)*(4*K.k_m.re[2][vband,cband]*(K.k_m.re[2][cband,vband]*K.k_m.Δ[2][cband,vband])+im*(K.hd.h_ops.a[2,2][cband,vband]*K.k_m.re[2][vband,cband]))
            data[:T2d][ij] = K.k_m.re[2][vband,cband]
        end
    end

    datafile = "$(mkpath("$(ENV["scriptdir"])/out/$RN"))/$asd-$(mn[1])-$(mn[2])-$Nbz.bson"
    bson(datafile,data)

    datafile
end

function shg_section(RN,asd,mn,Nbz)
    shg_section(;
        RN=RN,
        asd=asd,
        mn=mn,
        Nbz,
    )
end


######################################
#### Calculating DataIntegrals
######################################
function shg_section_direct(asd,mn,Nbz)
    SolidState.Main.models(asd,[mn])

    asdmn = BSON.load("$(ENV["scriptdir"])/.cache/$asd/asd-$(mn[1])-$(mn[2]).bson")
    asdg = (asdmn|>SolidState.ASDGeometry)

    hd = data_import("$(ENV["scriptdir"])/.cache/$asd/hd-$(mn[1])-$(mn[2]).bson")
    K = KinematicDensity(hd,[(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)])

    dim = size(K.hd.h_ops.h,1)
    vband = dim/2|>Int
    cband = vband+1

    patch_rng = range(-1.25*asdg["Kmag"],1.25*asdg["Kmag"], length=Nbz)

    data = Dict(
        :patch_rng => patch_rng,
        :asdg => asdg,
        :asd => asd,
        :mn  => mn,
        :T1d  => zeros(Complex{Float64},(Nbz,Nbz)),
        :T2d  => zeros(Complex{Float64},(Nbz,Nbz)),
    );

    for (i,kx) in enumerate(patch_rng)
        for (j,ky) in enumerate(patch_rng)
            ij = i + (j-1)*Nbz
            K([patch_rng[i],patch_rng[j]])
            data[:T1d][ij] = -im/(2.0*K.k_m.dω[cband,vband]^2)*(4*K.k_m.re[2][vband,cband]*(K.k_m.re[2][cband,vband]*K.k_m.Δ[2][cband,vband])+im*(K.hd.h_ops.a[2,2][cband,vband]*K.k_m.re[2][vband,cband]))
            data[:T2d][ij] = 2*(K.k_m.re[2][cband,vband]/K.k_m.dω[cband,vband]^2)*(4*K.k_m.re[2][vband,cband]*K.k_m.Δ[2][vband,cband]+im*K.hd.h_ops.a[2,2][cband,vband])
        end
    end

    datafile = "$(mkpath("$(ENV["scriptdir"])/out/$RN"))/$asd-$(mn[1])-$(mn[2])-$Nbz.bson"
    bson(datafile,data)

    datafile
end

function shg_section_direct(RN,asd,mn,Nbz)
    shg_section_direct(;
        RN=RN,
        asd=asd,
        mn=mn,
        Nbz,
    )
end

######################################
#### Calculating Shifted SHG
######################################
shifted_response_tag(asd,integral_info,rngs) = "shifted-$asd-$(hash(integral_info))-$(hash(rngs))"

function shifted_response(RN,asd,integral_info,rngs)
    (dtype,indices,priors,base,Nevals) = integral_info
    bdom    = getindex.(SolidState.range_scope(base),1)
    asd0 = asd()
    asdg = asd0|>SolidState.ASDGeometry

    maxshift = asdg["Lmag"]*sqrt(3)
    lc_z = asdg["xtal"][1][3,3]
    configurations = [ [ 0.0*[n1,n2,n3], [n1,n2,n3]] for n1 in (maxshift.*rngs[1]), n2 in (maxshift.*rngs[2]), n3 in (lc_z.*rngs[3])][:]

    data = Vector{Vector{ComplexF64}}(undef,length(configurations))
    for (i,shifts) in enumerate(configurations)
        shifted_asd = SolidState.LayerShiftASD(deepcopy(asd0),shifts)
        hd = TightBindingDensity(shifted_asd)
        dm = DataMap(()->shifted_asd,hd,dtype,indices,priors,base)
        di = DataIntegral(dm)
        di(Nevals:Nevals)

        data[i] = di.data[1][:]
    end

    bson("$(mkpath("$(ENV["scriptdir"])/out/$RN"))/shifted-$asd-$(hash(integral_info))-$(hash(rngs)).bson",
        Dict(
            :asd => asd,
            :Nevals => Nevals,
            :rngs => rngs,
            :bdom => bdom,
            :configurations => configurations,
            :data => data,
            :ucvol => det(asd0["blv"])
        )
    )
end
######################################
#### Shifted Sections
######################################

function shifted_response_section(RN,asd,dtype,indices,priors,base,Nevals,steps)
    bdom    = getindex.(SolidState.range_scope(base),1)

    asd0 = SolidState.LayerShiftASD(asd(),[zeros(3),[n1/N*asdg["Lmag"]*sqrt(3),n2/N*asdg["Lmag"]*sqrt(3),0.0]])
    asdg = asd0|>SolidState.ASDGeometry

    hd = TightBindingDensity(asd0)
    K = KinematicDensity(hd,[(:T,0.0,0.0,1),(:μ,0.0,0.0,1),(:δ,0.02,0.02,1)])

    dim = size(K.hd.h_ops.h,1)
    vband = dim/2|>Int
    cband = vband+1

    patch_rng = range(-1.25*asdg["Kmag"],1.25*asdg["Kmag"], length=Nbz)

    data = Dict(
        :patch_rng => patch_rng,
        :asdg => asdg,
        :asd => asd,
        :mn  => mn,
        :Td  => zeros(Complex{Float64},(Nbz,Nbz)),
        :Tv  => zeros(Complex{Float64},(Nbz,Nbz)),
    );

    for (i,kx) in enumerate(patch_rng)
        for (j,ky) in enumerate(patch_rng)
            ij = i + (j-1)*Nbz
            K([patch_rng[i],patch_rng[j]])
            for m in 1:dim
                for n in 1:dim
                    nm = n + (m-1)*dim
                    mn = m + (n-1)*dim
                    data[:Td][ij] += Complex(0,-3.0)*K.k_m.df[1][nm]/(K.k_m.dω[nm]^3)*K.k_m.re[2][nm]*(K.k_m.re[2][mn]*K.k_m.Δ[3][mn]+Complex(0.0,K.h_ops.a[2,2][mn]))
                    for l in 1:dim
                        ln = l + (n-1)*dim
                        ml = m + (l-1)*dim
                        data[:Tv][ij] += K.k_m[2][nm]*(K.h_ops.E.values[l]*(1.0/K.k_m.dω[ln]+1.0/K.k_m.dω[ml]+1.0/K.k_m.dω[nm])+0.5)*(K.k_m.re[2][ml]*K.k_m.re[2][ln])
                    end
                end
            end
        end
    end

    datafile = "$(mkpath("$(ENV["scriptdir"])/out/$RN"))/shifted-shg-section-$asd-$(mnN[1])-$(mnN[2])-$(mnN[3])-$Nbz.bson"
    bson(datafile,data)

end


function book_save(plotbook, plotdir)
    for key ∈ keys(plotbook)
        Plots.pdf(getindex(plotbook,key), "$plotdir/$key")
    end
end

function load(description::AbstractDict)
    dict = Dict()
    map(keys(description),values(description)) do name,RN
        push!(dict,name=>load(RN,name))
    end
    dict
end

function series_extract(RN,tags,plotargs::Tuple)
    map(enumerate(tags)) do (i,tag)
        SolidState.Main.extract(RN,tag,plotargs)
    end
end

function series_extract(RN,tags,plotargs::AbstractVector)
    map(enumerate(tags)) do (i,tag)
        SolidState.Main.extract(RN,tag,plotargs[i])
    end
end

function series_extract(RN,plottags::AbstractDict)
    map(collect(keys(plottags)),collect(values(plottags))) do (tag,plotarg)
        SolidState.Main.extract(RN,tag,plotarg...)
    end
end
