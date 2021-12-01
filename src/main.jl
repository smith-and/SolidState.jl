
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

function stat_print(stats)
    println("time   (s): $(stats.time)")
    println("bytes (gb): $(stats.bytes/1e9) \n")
    flush(stdout)
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
                println("making $mn");flush(stdout)
                stats = @timed if mn ∉ made_models || force
                        com_asd = SolidState.CommensurateASD(hs_asd,mn);
                        hd  = TightBindingDensity(com_asd)
                        bson("$rootdir/asd-$(mn[1])-$(mn[2]).bson",com_asd)
                        data_export("$rootdir/hd-$(mn[1])-$(mn[2]).bson",hd)
                else
                        println("$mn already made");flush(stdout)
                end
                stat_print(stats)
        end
        nothing
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
    odimeter       = Vector{Float64}[];
    corner_indices  = [zeros(Int64, length(path_list[i])) for i=1:length(path_list)];
    for i ∈ 1:length(path_list)
        positions = path_points(path_list[i], N)
        push!(paths,    positions)
        push!(odimeter, vcat([0.0],cumsum(LinearAlgebra.norm.(diff(positions)))))
        #println(indexin.(hash.(path_list[i]), Ref(hash.(positions))))
        # corner_indices[i] .= findindex.(Ref(positions), path_list[i])
        findcorners!(corner_indices[i], positions, path_list[i])
    end

    hcat(paths...), corner_indices, odimeter
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

#### disordered bands

function get_hsp_ticks(dm::DataMap)
    asdg = SolidState.ASDGeometry(SolidState.CommensurateASD(dm.input.asd(),dm.input.mn));
    tbks, corner_indices, odimeter = SolidState.Main.path_points(asdg["bz_hs"],[dm.input.base[1]],dm.input.base[2])
    tick_odimeters = [getindex.(Ref(odimeter[i]),corner_indices[i]) for i=1:length(corner_indices)]
    odimeter[1]./max(odimeter[1]...),(tick_odimeters[1]./max(odimeter[1]...), dm.input.base[1])
end

function band_average(RN::String, asd::Function, mn::Tuple{Int,Int}, npath::Int, pathlist::Vector{String}, nsample::Int, α, force::Bool=false)
    #Compute Data
    # SolidState.Main.models(asd,[mn],force)
    data = @sync @distributed for i in 1:nsample
        asdmn       = BSON.load(   "$(ENV["cachedir"])/$asd/asd-$(mn[1])-$(mn[2]).bson")
        rasd = SolidState.randomize_hopping!(α,asdmn)
        dm = DataMap(()->rasd,BANDS,(:bandstructure,),[(:μ,0.0,0.0,1)],(pathlist,npath))
        stats = @timed dm()
        stat_print(stats)
        bson("$(mkpath("$(ENV["scriptdir"])/out/$RN/samples-$(mn[1])-$(mn[2])"))/sample-$i.bson",Dict(:data=>dm.chart.data[:,1,1,:]))
        nothing
    end

    dm = DataMap(()->BSON.load(   "$(ENV["cachedir"])/$asd/asd-$(mn[1])-$(mn[2]).bson"),BANDS,(:bandstructure,),[(:μ,0.0,0.0,1)],(pathlist,npath))
    avg = zeros(typeof(dm.chart.data[1]),size(dm.chart.data[:,1,1,:]))
    for i in 1:nsample
        avg .+= BSON.load("$(ENV["scriptdir"])/out/$RN/samples-$(mn[1])-$(mn[2])/sample-$i.bson")[:data]
    end
    avg ./= nsample
    dm.chart.data[:,1,1,:] .= avg

    std = zeros(eltype(dm.chart.data[:,1,1,:]),size(dm.chart.data[:,1,1,:]))
    for i in 1:nsample
        dat = BSON.load("$(ENV["scriptdir"])/out/$RN/samples-$(mn[1])-$(mn[2])/sample-$i.bson")[:data]
        std .+= (dat.-avg).^2
    end
    rm("$(ENV["scriptdir"])/out/$RN/samples-$(mn[1])-$(mn[2])",recursive=true)
    #Package Data
    dict = Dict(
        :chart => dm.chart,
        :tick_info => get_hsp_ticks(dm),
        :ribbon => sqrt.(std./nsample),
        :θ => SolidState.cθ(mn...)*180/π,
    )
    #Export Data
    bson("$(mkpath("$(ENV["scriptdir"])/out/$RN/"))/$asd-$(mn[1])-$(mn[2])-$nsample-$α.bson",dict)
    #Return Extraction Tag
    "$asd-$(mn[1])-$(mn[2])-$nsample-$α"
end

######################################
#### Calculating 1D path sections with full structure data retained
######################################
function hd_section_1d(asd::Dict{String,Any}, hd::HamiltonianDensity, Npath=100, pathlist=["K1","Γ","M1","K'1"])
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
            gridalpha      = .5 ,
            gridlinewidth  = 1,
            gridstyle      = :dashdot,
            xtickfonthalign =:center,
            xticks         = (tick_odimeters[1]./max(odimeter[1]...), pathlist),
            minorgrid      = false,
            apsectratio    = 1/4,
            xguidefontsize = 16,
            title_pos       = :right,
            titlefontvalign = :bottom,
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
    dict      = hd_section_1d(ASD,hd,Npath,pathlist)

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

integral_tag(asd,mn,chart_integral_info) = "$asd-$(mn[1])-$(mn[2])-$(hash(chart_integral_info))"

function integral(RN::String, asd, mn::Tuple, chart_integral_info::Tuple, pool=default_worker_pool(), cachedir=ENV["cachedir"], scriptdir=ENV["scriptdir"]; force=false)

    di0 = DataIntegral(asd, mn, chart_integral_info[1:end-1]...)
    stats = @timed (di0 = di(chart_integral_info[end], pool))
    stat_print(stats)

    bson("$(mkpath("$scriptdir/out/$RN"))/$asd-$(mn[1])-$(mn[2])-$(hash(chart_integral_info)).bson", Dict(
        :chart_integral_info => chart_integral_info,
        :mn => mn,
        :θ => SolidState.cθ(mn...)*180/π,
        :asd => asd,
        :base => getindex.(di0.dm.chart.base,1),
        :data => di0.data,
        :ribbon  => di0.err,
        :npool => length(pool),
    ))

    "$asd-$(mn[1])-$(mn[2])-$(hash(chart_integral_info))"
end

function integral_average(RN::String, asd::Function, mn::Tuple{Int,Int}, chart_integral_info::Tuple, nsample::Int, α::AbstractFloat, pool::AbstractWorkerPool=default_worker_pool(), force::Bool=false)
    data = pmap(1:nsample, batch_size = Int(ceil(nsample/nworkers()))) do _
        asdmn       = BSON.load(   "$(ENV["cachedir"])/$asd/asd-$(mn[1])-$(mn[2]).bson")
        rasd = SolidState.randomize_hopping!(α,asdmn)
        dm = DataMap(()->rasd,chart_integral_info[1:end-1]...)
        di = DataIntegral(dm)
        stats = @timed di(chart_integral_info[end])
        stat_print(stats)
        di.data[1]
    end
    avg = sum(data)./length(data)
    std = sqrt.(sum(map(data) do dat
        (dat.-avg).^2
    end)./length(data))
    dm = DataMap(()->BSON.load(   "$(ENV["cachedir"])/$asd/asd-$(mn[1])-$(mn[2]).bson"),chart_integral_info[1:end-1]...)

    #Package Data
    bson("$(mkpath("$(ENV["scriptdir"])/out/$RN"))/$asd-$(mn[1])-$(mn[2])-$nsample-$α-$(hash(chart_integral_info)).bson", Dict(
        :chart_integral_info => chart_integral_info,
        :mn => mn,
        :θ => SolidState.cθ(mn...)*180/π,
        :asd => asd,
        :base => getindex.(dm.chart.base,1),
        :data => [avg],
        :ribbon => std./2,
        :npool => length(pool),
    ))
    #Return Extraction Tag
    "$asd-$(mn[1])-$(mn[2])-$nsample-$α-$(hash(chart_integral_info))"
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
                :comargs => models(asd,SolidState.(:bulkhead, comargs...)),
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
            :comargs => models(asd,SolidState.(:bulkhead, comargs...)),
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
            :comargs => models(asd,SolidState.(:bulkhead, comargs...)),
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

function load(RN,name)
        BSON.load("$(ENV["scriptdir"])/out/$RN/$name.bson")
end

function extract(RN,name,plotfunction,args...;)
    plotfunction(args...; BSON.load("$(ENV["scriptdir"])/out/$RN/$name.bson")...)
end

function extract(;RN,mainf,extractf,asd,mn,main_info,force=false)
    if isfile("$(ENV["scriptdir"])/out/$RN/$asd-$(mn[1])-$(mn[2])-$(hash(main_info))")&&(!force)
        # SolidState.Main.load(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(integral_info))")
        @assert !force
        SolidState.Main.extract(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(main_info))",extractf)
    else
        mainf(RN,asd,mn,main_info)
        # SolidState.Main.load(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(integral_info))")
        SolidState.Main.extract(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(main_info))",extractf)
    end
end

function extract(extract_info;RN,mainf,extractf,asd,mn,main_info,force=false)
    if isfile("$(ENV["scriptdir"])/out/$RN/$asd-$(mn[1])-$(mn[2])-$(hash(main_info))")&&(!force)

        # SolidState.Main.load(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(integral_info))")
        @assert !force
        SolidState.Main.extract(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(main_info))",extractf,extract_info...)
    else
        mainf(RN,asd,mn,main_info)
        # SolidState.Main.load(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(integral_info))")
        SolidState.Main.extract(RN,"$asd-$(mn[1])-$(mn[2])-$(hash(main_info))",extractf,extract_info...)
    end
end

end
