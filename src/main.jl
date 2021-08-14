
module Main

using Revise
using Distributed, Dates, OrderedCollections, BSON, Plots, PackageCompiler
using LinearAlgebra, SharedArrays, StaticArrays
using CubicSplines, Roots, SpecialFunctions, HCubature
using Base: Threads
using ..SolidState

function test(args...)
    println(args...)
end

###########################################
#### Compilation Methods
###########################################

### Run the precompile script
function test_compile_run()
    include("$(@__DIR__)/precompile_script.jl")
end

"""
    base_compile(cachedir::String=pwd())

Add some basics ([:BSON, :OrderedCollections, :LinearAlgebra, :Plots, :ColorSchemes]) to the current default system image
"""
function base_compile()
    println("compiling a system image")
    stats = @timed PackageCompiler.create_sysimage([:BSON, :OrderedCollections, :LinearAlgebra, :Plots],
        #project                     = "$(ENV["HOME"])/.julia/environments/v1.5/",
        incremental                 = true,
        replace_default = true,
    )
    println("build time \t"*string(stats.time/60)*"m")
    nothing
end


"""
    prestack(cachedir::String=pwd())

Make sysimage for modules [:Plots,:SolidState,:SolidStateApps] and in cachedir
"""
function sysimage(workingdir::String=ENV["cachedir"])
    println("compiling a system image")
    cdir = mkpath("$workingdir/system");
    stats = @timed PackageCompiler.create_sysimage([:BSON, :OrderedCollections, :LinearAlgebra, :Plots, :SolidState],
        sysimage_path               = "$cdir/sysimage.dylib",
        precompile_execution_file   = "$(@__DIR__)/precompile_script.jl",
        #project                     = "$(ENV["HOME"])/.julia/environments/v1.5/",
        incremental                 = true
    )
    println("build time \t"*string(stats.time/60)*"m")
    nothing
end

######################################
#### Make Tightbinding Models
######################################

function first_series(mmin,mmax)
        args = Vector{Tuple{Int64,Int64}}(undef,length(mmin:mmax))
        for (i,n) ∈ enumerate(mmin:mmax)
                args[i] = (1,n)
        end
        args
end

function principal_series(mmin,mmax)
        args = Vector{Tuple{Int64,Int64}}(undef,length(mmin:mmax))
        for (i,m) ∈ enumerate(mmin:mmax)
                args[i] = (m,m+1)
        end
        args
end

dimx((m,n)) = (SolidState.TwistedTriangularGeometry(ASD2()["blv"],(m,n))["blv"]|>det)/(ASD2()["blv"]|>det)
function hull_series(mmin,mmax)
        args = union(vcat(principal_series(mmin,mmax),first_series(mmin,mmax)))
        sizeargs = dimx.(args)

        args[sortperm(sizeargs)]
end

function bulkhead_series(dmin,dmax)
        #need way to relate these to the dcut
        (mm,sm) = (200,200)
        θid = [(m,m+s) for m∈0:mm for s∈0:sm]
        θsp = [SolidState.cθ(m,m+s)*180/π for m∈0:mm for s∈0:sm]
        θLM = [dimx((m,m+s)) for m∈0:mm for s∈0:sm]

        LMmask= 0 .< θLM .< dmax
        unqθ = (union(round.(θsp[LMmask],digits=4)))
        θspargs = [findall(x->round(x,digits=4)==θ,θsp[LMmask]) for θ∈unqθ]
        θLMcopies=getindex.(Ref(θLM[LMmask]),θspargs)
        θidcopies=getindex.(Ref(θid[LMmask]),θspargs)

        θhullid = Vector{typeof(θid[1])}(undef,length(unqθ))
        θhullLM = Vector{typeof(θLM[1])}(undef,length(unqθ))
        for (i,LMs) ∈ enumerate(θLMcopies)
                θhullid[i] = θidcopies[i][argmin(LMs)]
                θhullLM[i] = θLMcopies[i][argmin(LMs)]
        end

        return θhullid[dmin .< θhullLM .< dmax][sortperm(θhullLM[dmin .< θhullLM .< dmax])]
end

### Commensurate Arg Plot
function full_series(dmin,dmax)
        #need way to relate these to the dcut
        (mm,sm) = (200,200)
        θid = [(m,m+s) for m∈-mm:mm for s∈-sm:sm]
        θsp = [SolidState.cθ(m,m+s)*180/π for m∈-mm:mm for s∈-sm:sm]
        θLM = [dimx((m,m+s)) for m∈-mm:mm for s∈-sm:sm]

        LMmask= 0 .< θLM .< dmax
        unqθ = (union(round.(θsp[LMmask],digits=4)))
        θspargs = [findall(x->round(x,digits=4)==θ,θsp[LMmask]) for θ∈unqθ]
        θLMcopies=getindex.(Ref(θLM[LMmask]),θspargs)
        θidcopies=getindex.(Ref(θid[LMmask]),θspargs)

        θhullid = Vector{typeof(θid[1])}(undef,length(unqθ))
        θhullLM = Vector{typeof(θLM[1])}(undef,length(unqθ))
        for (i,LMs) ∈ enumerate(θLMcopies)
                θhullid[i] = θidcopies[i][argmin(LMs)]
                θhullLM[i] = θLMcopies[i][argmin(LMs)]
        end

        return θhullid[dmin .< θhullLM .< dmax][sortperm(θhullLM[dmin .< θhullLM .< dmax])]
end

function twist_series(series,mmin,mmax)
        if series==:first
                return first_series(mmin, mmax)
        elseif series==:principal
                return principal_series(mmin, mmax)
        elseif series==:hull
                return hull_series(mmin,mmax)
        elseif series==:bulkhead
                return bulkhead_series(mmin,mmax)
        elseif series==:full
                return full_series(mmin,mmax)
        elseif series==:mn
                return [(mmin,mmax)]
        end
end

function model_id(x)
        try
                s1,s2 = findall(isequal('-'),x)
                s3,   = findall(isequal('.'),x)
                return (parse(Int,x[(s1+1):(s2-1)]),parse(Int,x[(s2+1):(s3-1)]))
        catch
                return nothing
        end
end

"""
    models(asd, comargs::Vector{Tuple{Int,Int}}, force=false, cachedir=ENV["cachedir"], kargs...)

Makes the `ASD` and `HamiltonianDensity` data objects and saves them to bson files in the directory `cachedir="$(@__FILE__)/cache/asd"`.
"""
function models(asd::Symbol, comargs::Vector{Tuple{Int,Int}}, force=false, cachedir=ENV["cachedir"], kargs...)
    models(eval(quote $asd end), comargs, force, cachedir, kargs...)
end
function models(asd::Function, comargs::Vector{Tuple{Int,Int}}, force=false, cachedir=ENV["cachedir"], kargs...)
        rootdir = isdir("$cachedir/$asd") ? "$cachedir/$asd" : mkpath("$cachedir/$asd");
        hs_asd    = asd()
        made_models = model_id.(readdir(rootdir))
        println("");flush(stdout)
        println("Making $asd Models in $cachedir");flush(stdout)
        for mn ∈ comargs
                if mn ∉ made_models || force
                        com_asd = SolidState.CommensurateASD(hs_asd,mn);
                        hd  = TightBindingDensity(com_asd)
                        bson("$rootdir/asd-$(mn[1])-$(mn[2]).bson",com_asd)
                        data_export("$rootdir/hd-$(mn[1])-$(mn[2]).bson",hd)
                        println("$mn made");flush(stdout)
                else
                        println("$mn already made");flush(stdout)
                end
        end
        bson("$rootdir/comargs.bson",Dict(:cargs=>model_id.(readdir(rootdir))))
        comargs
end

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

######################################
#### Calculating BandStructures
######################################

#Piecewise-Linear Path (1D) Discretization
function path_points(path_corners::Vector{Vector{Float64}}, N0::Int64)::Array{Array{Float64,1},1}
    N=N0-1
    #obtaining the difference vectors between subsequent path points
    relative_vecs::Vector{Vector{Float64}}=[path_corners[i+1]-path_corners[i] for i=1:(length(path_corners)-1)]
    #Obtaining the number of points in each leg for approximately N overall
    p_lengths = LinearAlgebra.norm.(relative_vecs)
    p_length  = sum(p_lengths)
    leg_Ns = max.(floor.(p_lengths ./ p_length .* N ),Ref(1))
    i=0; while N-sum(leg_Ns) > 0
        leg_Ns[i%length(leg_Ns)+1] += 1
    end
    #Calculating the points along each leg of the path
    path_positions = push!(vcat([ [path_corners[leg_idx] + relative_vecs[leg_idx]*pt_idx/leg_Ns[leg_idx] for pt_idx=0:(leg_Ns[leg_idx]-1) ] for leg_idx=1:length(leg_Ns)]...),path_corners[end])

    path_positions
end



function findindex(pathlist,point)
    for i=1:length(pathlist)
        LinearAlgebra.norm(pathlist[i]-point) < 1e-10 ? (return i; break) : nothing
    end
    return 0 #this is bad
end

function findcorners!(corner_indices, pathlist, corner_list)
    j = 1
    for i=1:length(pathlist)
        for corner in unique(corner_list)
            if norm(pathlist[i]-corner) < 1e-10
                corner_indices[j]=i
                j+=1
            end
        end
    end
end

function path_points(path_list::Vector{Vector{Vector{Float64}}}, N::Int64)::(Tuple{Array{Vector{Float64},2},Array{Array{Int64,1},1},Array{Array{Float64,1},1}})
    paths           = Vector{Vector{Float64}}[];
    odometers       = Vector{Float64}[];
    corner_indices  = [zeros(Int64, length(path_list[i])) for i=1:length(path_list)];
    for i ∈ 1:length(path_list)
        positions = path_points(path_list[i], N)
        push!(paths,    positions)
        push!(odometers, vcat([0.0],cumsum(LinearAlgebra.norm.(diff(positions)))))
        #println(indexin.(hash.(path_list[i]), Ref(hash.(positions))))
        # corner_indices[i] .= findindex.(Ref(positions), path_list[i])
        findcorners!(corner_indices[i], positions, path_list[i])
    end

    hcat(paths...), corner_indices, odometers
end

function path_points(bz_pt_ledger::Dict{String,Array{Float64,1}}, path_list, N::Int64)::(Tuple{Array{Vector{Float64},2},Array{Array{Int64,1},1},Array{Array{Float64,1},1}})
    path_points([get.(Ref(bz_pt_ledger), path_corners, Ref([Inf,Inf])) for path_corners ∈ path_list], N)
end

function path_points(asd, path_list, N::Int64)#::(Tuple{Array{Vector{Float64},2},Array{Array{Int64,1},1},Array{Array{Float64,1},1}})
    bz_pt_ledger = SolidState.ASDGeometry(asd)["bz_hs"]
    path_points([get.(Ref(bz_pt_ledger), path_corners, Ref([Inf,Inf])) for path_corners ∈ path_list], N)
end


#Simple Band Structure
function band_data(asd::Dict{String,Any}, hd::HamiltonianDensity, pathlist=["K1","Γ","M1","K'1"],Npath = 100; kargs...)
    asdg = SolidState.ASDGeometry(asd);
    tbks, corner_indices, odimeter = path_points(asdg["bz_hs"],[pathlist],Npath)

    #Compute Data
    tbEs = zeros(Float64,(size(hd.h_ops.h,1),Npath));
    for (i,k) ∈ enumerate(tbks)
        hd(k)
        tbEs[:,i] .= sort(eigvals(Hermitian(hd.h_ops.h)))
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
            :args=>args,
            :ks=>tbks,
            :path=>pathlist,
            )
end


"""
    bands(asd::Symbol, mn::Tuple{Int,Int}, RN, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...)

Calculate the band structure of a model along the high symmetry points listed in
"""
function bands(RN, asd, mn::Tuple{Int,Int}, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...;force = false)
    #Compute Data
    models(asd,[mn],force)
    ASD       = BSON.load(   "$cachedirr/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd        = data_import( "$cachedirr/$asd/hd-$(mn[1])-$(mn[2]).bson")
    dict = band_data(ASD,hd,pathlist,Npath; title = "$(round(180/π*SolidState.cθ(mn...),digits=3))ᵒ")

    #Export Data
    bson("$(mkpath("$scriptdir/out/$RN/"))/$asd-$(mn[1])-$(mn[2]).bson",
        Dict(
                :dict=>dict
        )
    )

    "$scriptdir/out/$RN/$asd-$(mn[1])-$(mn[2]).bson"
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

######################################
#### Calculating 1D path sections with full structure data retained
######################################
function hd_section_1d(asd::Dict{String,Any}, hd::HamiltonianDensity, Npath=100, pathlist=["K1","Γ","M1","K'1"]; kargs...)
    asdg = SolidState.ASDGeometry(asd);
    tbks, corner_indices, odimeter = path_points(asd,[pathlist],Npath)

    #Compute Data
    Es = zeros(Float64,(size(hd.h_ops.h,1),Npath));
    Evs = Vector{typeof(hd.h_ops.E.vectors)}(undef,Npath);
    # Vs = Vector{typeof(hd.h_ops.v)}(undef,Npath)
    # As = Vector{typeof(hd.h_ops.a)}(undef,Npath)
    for (i,k) ∈ enumerate(tbks)
        hd(k)
        E = eigen(Hermitian(hd.h_ops.h))
        Es[:,i] = E.values
        Evs[i] = E.vectors
        # SolidState.cartan_transform(hd.h_ops)
        # Es[:,i] = deepcopy(hd.h_ops.E.values)
        # Evs[i] = deepcopy(hd.h_ops.E.vectors)
        # Vs[i] = deepcopy(hd.h_ops.v)
        # As[i] = deepcopy(hd.h_ops.a)
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

    Dict(
        :odimeter => odimeter./max(odimeter[1]...),
        :Es       => Es',
        :args     => args,
        :ks       => tbks,
        :path     => pathlist,
        :Evs      => Evs,
        # :Vs       => Vs,
        # :As       => As,
    )
end

hds1d_tag(asd,mn,Npath,pathlist) = "bands-$asd-$(mn[1])-$(mn[2])-$Npath-$(hash(pathlist))"
function hd_section_1d(RN, asd, mn, Npath = 300, pathlist=["K1","Γ","M1","K'1"], scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...;force = false)
    #Compute Data
    ASD       = BSON.load(   "$cachedirr/$asd/asd-$(mn[1])-$(mn[2]).bson")
    hd        = data_import( "$cachedirr/$asd/hd-$(mn[1])-$(mn[2]).bson")
    dict      = hd_section_1d(ASD,hd,Npath,pathlist; title = "$(round(180/π*SolidState.cθ(mn...),digits=3))ᵒ")

    #Export Data
    bson("$(mkpath("$scriptdir/out/$RN/"))/$(hds1d_tag(asd,mn,Npath,pathlist)).bson",
        Dict(
            :dict  => dict,
            :tag   => hds1d_tag(asd,mn,Npath,pathlist),
            :angle => SolidState.cθ(mn...)*180/π
        )
    )
    "$(mkpath("$scriptdir/out/$RN/"))/$(hds1d_tag(asd,mn,Npath,pathlist)).bson"
end

######################################
#### Calculating Spectral Sections
######################################

function spectral_section(;RN,asd,mn,Nbz)
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
        :dω  => zeros(Float64,(Nbz,Nbz)),
        :re1 => zeros(Complex{Float64},(Nbz,Nbz)),
        :re2 => zeros(Complex{Float64},(Nbz,Nbz)),
        :Δ1  => zeros(Complex{Float64},(Nbz,Nbz)),
        :Δ2  => zeros(Complex{Float64},(Nbz,Nbz)),
    );

    for (i,kx) in enumerate(patch_rng)
        for (j,ky) in enumerate(patch_rng)
            ij = i + (j-1)*Nbz
            K([patch_rng[i],patch_rng[j]])
            data[:re1][ij] = K.k_m.re[1][vband,cband]
            data[:re2][ij] = K.k_m.re[2][vband,cband]
            data[:dω][ij]  = K.k_m.dω[vband,cband]
            data[:Δ1][ij]  = K.k_m.Δ[1][vband,cband]
            data[:Δ2][ij]  = K.k_m.Δ[2][vband,cband]
        end
    end

    datafile = "$(mkpath("$(ENV["scriptdir"])/out/$RN"))/$asd-$(mn[1])-$(mn[2])-$Nbz.bson"
    bson(datafile,data)

    datafile
end

function spectral_section(RN,asd,mn,Nbz)
    spectral_section(;
        RN=RN,
        asd=asd,
        mn=mn,
        Nbz,
    )
end

######################################
#### Calculating DataIntegrals
######################################


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

function shifted_response(RN,asd,dtype,indices,priors,base,Nevals,rngs)

    bdom    = getindex.(SolidState.range_scope(base),1)
    asd0 = asd()
    asdg = asd0|>SolidState.ASDGeometry

    maxshift = asdg["Lmag"]*sqrt(3)
    lc_z = asdg["xtal"][1][3,3]
    configurations = [ [ -0.5*[n1,n2,n3], 0.5*[n1,n2,n3]] for n1 in (maxshift.*rngs[1]), n2 in (maxshift.*rngs[2]), n3 in (lc_z.*rngs[3])][:]

    data = Array{Vector{ComplexF64}}(undef,length(configurations))
    for (i,shifts) in enumerate(configurations)
        shifted_asd = SolidState.LayerShiftASD(deepcopy(asd0),shifts)
        hd = TightBindingDensity(shifted_asd)
        dm = DataMap(dtype,shifted_asd,hd,indices,priors,base)
        di = DataIntegral(dm)
        di(Nevals:Nevals)

        data[i] = di.data[1][:]
    end

    bson("$(mkpath("$(ENV["scriptdir"])/out/$RN"))/shifted-$asd-$Nevals"*"$(string(("-".*string.(rngs))...)).bson",
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

######################################
#### Calculating DataIntegrals
######################################

function integral(RN, asd, mn, dtype::Type{T} where T <: SolidState.DataChart, indices, priors, base, Neval, pool=default_worker_pool(), cachedir=ENV["cachedir"], scriptdir=ENV["scriptdir"]; force=false)

    println("Calculating Integral for $asd $(mn[1])-$(mn[2]) with $Neval points");flush(stdout)
    models(asd,[mn])
    di0 = DataIntegral(asd, mn, dtype, indices, priors, base)

    di = di0(Neval, pool)

    bson("$(mkpath("$scriptdir/out/$RN"))/$asd-$(mn[1])-$(mn[2])-$dtype-$Neval.bson", Dict(
        :mn => mn,
        :asd => asd,
        :dtype => dtype,
        :indicess => indices,
        :priors =>priors,
        :base => base,
        :Neval => Neval,
        :scriptdir => scriptdir,
        :cachedir => cachedir,
        :plotdir => mkpath("$scriptdir/out/$RN"),
        :handle => "$asd-$(mn[1])-$(mn[2])-$dtype-$Neval",
        :data => di.data,
        :err  => di.err,
        :evals => di.evals,
        :base => getindex.(getfield(di.dm.chart,1).base,1),
        :npool => length(pool)
    ))

    "$scriptdir/out/$RN/$asd-$(mn[1])-$(mn[2])-$dtype-$Neval.bson"
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

######################################
#### Scaling Routines
######################################

function dm_scaling(RN, asd, comargs, datatype, indices, priors, base)

        args = Dict(
                :asd => asd,
                :datatype => datatype,
                :indices => indices,
                :priors  => priors,
                :base => base,
                :comargs => comargs,
                :cachedir => "$(ENV["cachedir"])",
                :datadir => "$(ENV["scriptdir"])/out/$RN"|>mkpath,
                :plotdir => "$(ENV["scriptdir"])/plot/$RN"|>mkpath
        )
        SolidState.Scaling.dm_scaling(;args...)
end

function mem_scaling(RN, asd, comargs::Tuple{Int64,Int64}, datatype, indices, priors, base)

        args = Dict(
                :asd => asd,
                :datatype => datatype,
                :indices => indices,
                :priors  => priors,
                :base => base,
                :comargs => models(asd,twist_series(:bulkhead, comargs...)),
                :cachedir => "$(ENV["cachedir"])",
                :datadir => "$(ENV["scriptdir"])/out/$RN"|>mkpath,
                :plotdir => "$(ENV["scriptdir"])/plot/$RN"|>mkpath
        )
        SolidState.Scaling.mem_test(;args...)
end

function core_scaling(RN, asd, comargs, datatype, indices, priors, base, Neval)
    args = Dict(
            :asd => asd,
            :datatype => datatype,
            :indices => indices,
            :priors  => priors,
            :base => base,
            :comargs => models(asd,twist_series(:bulkhead, comargs...)),
            :cachedir => "$(ENV["cachedir"])",
            :datadir => "$(ENV["scriptdir"])/out/$RN"|>mkpath,
            :plotdir => "$(ENV["scriptdir"])/plot/$RN"|>mkpath,
            :pool => default_worker_pool(),
            :Neval => Neval
    )

    SolidState.Scaling.core_scaling(;args...)
end

function integral_convergence(RN, asd, comargs, datatype, f, sindices, priors, base, evals)
    args = Dict(
            :asd => asd,
            :datatype => datatype,
            :indices => indices,
            :priors  => priors,
            :base => base,
            :comargs => models(asd,twist_series(:bulkhead, comargs...)),
            :cachedir => "$(ENV["cachedir"])",
            :datadir => "$(ENV["scriptdir"])/out/$RN"|>mkpath,
            :plotdir => "$(ENV["scriptdir"])/plot/$RN"|>mkpath,
            :pool => default_worker_pool(),
            :f => f,
            :evals => evals,
    )

    SolidState.Scaling.julia_worker_test(;args...)
end



# function blas_julia_thread_tradeoff(; asd, datatype, indices, priors, base, comargs, cachedir, datadir, plotdir, pool = default_worker_pool(), Neval, blas_thread_max=length(pool), kargs...)
# function blas_thread_map_test(; asd, datatype, indices, priors, base, comargs, cachedirr, scriptdir, blas_thread_max, kargs...)
# function blas_thread_integral_test(asd, datatype, indices, priors, base, args...; comargs, scriptdir, blas_thread_max, Neval, pool, kargs...)

###############################################################
###### Extraction Methods
###############################################################

function book_save(plotbook, plotdir)
    for key ∈ keys(plotbook)
        Plots.pdf(getindex(plotbook,key), "$plotdir/$key")
    end
end

function load(RN,name)
        BSON.load("$(ENV["scriptdir"])/out/$RN/$name.bson")
end

function load(description::AbstractDict)
    dict = Dict()
    map(keys(description),values(description)) do name,RN
        push!(dict,name=>load(RN,name))
    end
    dict
end

function extract(RN,name,plotfunction,args...)
    datafile = "$(ENV["scriptdir"])/out/$RN/$name.bson"
    plotdir = "$(ENV["scriptdir"])/plot/$RN"|>mkpath
    plotfunction(args...; BSON.load(datafile)..., plotdir=plotdir)
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

function b2_extract(RN,name,plotfunction,args...)
    plotdir = "$(ENV["scriptdir"])/plot/b2-$RN"|>mkpath
    source = "$(ENV["b2scriptdir"])/out/$RN/$name.bson"
    target = "$(mkpath("$(ENV["scriptdir"])/out/b2-$RN"))/$name.bson"
    run(`rsync -r --progress asmithc@bridges2.psc.edu:$source $target`)
    plotfunction(args...; BSON.load(target)..., plotdir=plotdir)
end

end
