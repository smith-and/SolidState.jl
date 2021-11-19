######################################
#### Model Making
######################################
function dictstring(p)
    handle = string(map(keys(p),values(p)) do name,val
        if typeof(val) <: Tuple
            "-$name-"*string((",".*string.(val))...)[2:end]
        else
            "-$name-$val"
        end
    end...)[2:end]

    if length(handle) > 100
        string(hash(handle))
    else
        handle
    end
end

function make_model(p,cachedir=ENV["cachedir"])
    if !isfile("$(mkpath("$cachedir/ASD"))/$(dictstring(p)).bson")
        bson("$(mkpath("$cachedir/ASD"))/$(dictstring(p)).bson",form_model(p))
    end
    BSON.load("$(mkpath("$cachedir/ASD"))/$(dictstring(p)).bson")
end

function make_tb_model(p,cachedir=ENV["cachedir"])
    if isfile("$(mkpath("$cachedir/TBHD"))/$(dictstring(p)).bson")
        println("Already Made!")
    elseif isfile("$cachedir/ASD/$(dictstring(p)).bson")
        data_export("$(mkpath("$cachedir/TBHD"))/$(dictstring(p)).bson",TightBindingDensity(load("ASD",p)))
    else
        make_model(p,cachedir)
        make_tb_model(p,cachedir)
    end
    data_import("$(mkpath("$cachedir/TBHD"))/$(dictstring(p)).bson")
end

function load(tag,p,cachedir=ENV["cachedir"])
    if tag=="ASD"
        BSON.load("$cachedir/$tag/$(dictstring(p)).bson")
    else
        data_import("$cachedir/$tag/$(dictstring(p)).bson")
    end
end

######################################
#### Band Structures
######################################
#Simple Band Structure
function full_band_data(asd::Dict{String,Any}, hd::HamiltonianDensity, pathlist=["K1","Γ","M1","K'1"],Npath = 100; kargs...)
    asdg = SolidState.ASDGeometry(asd);
    tbks, corner_indices, odimeter = SolidState.Main.path_points(asdg["bz_hs"],[pathlist],Npath)

    #Compute Data
    tbEs = zeros(Float64,(size(hd.h_ops.h,1),Npath));
    tbVs = zeros(Complex{Float64},(size(hd.h_ops.h,1),size(hd.h_ops.h,1),Npath));
    for (i,k) ∈ enumerate(tbks)
        hd(k)
        E = eigen(Hermitian(hd.h_ops.h))
        tbEs[:,i] .= E.values
        tbVs[:,:,i] .= E.vectors
    end

    tick_odimeters = [getindex.(Ref(odimeter[i]),corner_indices[i]) for i=1:length(corner_indices)]

    args = (
            frame   = :box,
            grid    = :x,
            label   = "",
            xlims   = (0, 1),
            yrotation=60,
            yguide  = "eV",
            gridalpha      = .5 ,
            gridlinewidth  = 1,
            gridstyle      = :dashdot,
            xtickfonthalign =:center,
            xtickfontsize  = 9,
            yguidefontsize = 10,
            xticks         = (tick_odimeters[1]./max(odimeter[1]...), pathlist),
            minorgrid      = false,
            margin         = 4Plots.mm,
            apsectratio    = 1/4,
            title_pos       = :right,
            titlefontvalign = :bottom,
            kargs...
        );

    return Dict(
            :odimeter=>odimeter./max(odimeter[1]...),
            :Es=>tbEs',
            :Vs=>tbVs,
            :args=>args,
            :ks=>tbks,
            :path=>pathlist,
        )
end

"""
    bands(asd::Symbol, mn::Tuple{Int,Int}, RN, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...)

Calculate the band structure of a model along the high symmetry points listed in
"""
function bands(RN, pmodel::OrderedDict, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...;force = false)
    #Compute Data
    asd       = make_model(pmodel)
    hd        = make_tb_model(pmodel)
    dict      = full_band_data(asd,hd,pathlist,Npath; title = "$(round(180/π*SolidState.cθ(pmodel[:twist]...),digits=3))ᵒ")

    #Export Data
    bson("$(mkpath("$scriptdir/out/$RN/"))/bands$(dictstring(pmodel)).bson",
        Dict(
            :dict   => dict,
        )
    )

    "$scriptdir/out/$RN/bands$(dictstring(pmodel)).bson"
end

######################################
#### Adaptive Integrals
######################################

function integral(RN, pmodel, pchart, Neval, pool=default_worker_pool(), cachedir=ENV["cachedir"], scriptdir=ENV["scriptdir"]; force=false)
    println("Calculating Integral with $Neval points");flush(stdout)
    di = DataIntegral(DataMap(pmodel,pchart))
    if length(pool)==0
        di(Neval)
    else
        di(Neval,pool)
    end

    bson("$(mkpath("$scriptdir/out/$RN"))/$(dictstring(pmodel))-$(dictstring(pchart))-$Neval.bson", Dict(
        :pmodel => pmodel,
        :pchart => pchart,
        :Neval => Neval,
        :scriptdir => scriptdir,
        :cachedir => cachedir,
        :plotdir => mkpath("$scriptdir/out/$RN"),
        :handle => "$(dictstring(pmodel))-$(dictstring(pchart))",
        :data => di.data,
        :err  => di.err,
        :evals => di.evals,
        :base => getindex.(getfield(di.dm.chart,1).base,1),
    ))

    "$(mkpath("$scriptdir/out/$RN"))/$(dictstring(pmodel))-$(dictstring(pchart))-$Neval.bson"
end

######################################
#### Calculating DataSections
######################################

function section(RN, asd, mn, dtype, indices, priors, base, Neval, pool=default_worker_pool(), offset=(0.0,0.0), cachedir=ENV["cachedir"], scriptdir=ENV["scriptdir"])

    println("Calculation Section for $asd $(mn[1])-$(mn[2])");flush(stdout)
    models(asd,[mn])
    args = Dict(
        :asd=>asd,
        :mn => mn,
        :dtype => dtype,
        :N => Neval,
        :indices => indices,
        :priors => priors,
        :base => base,
        :cachedir => cachedir,
        :offset => offset,
        )

    ds = DataSection(; args...)
    ds()

    outdir = mkpath("$scriptdir/out/$RN")
    data_export("$outdir/$asd-$(mn[1])-$(mn[2]).bson", ds);
end
