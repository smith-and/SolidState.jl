
module Main

using Revise
using Distributed, Dates, OrderedCollections, BSON, Plots, PackageCompiler
using LinearAlgebra, SharedArrays, StaticArrays
using Mux, WebIO, Interact, InteractiveUtils
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

function twist_series(series,mmin,mmax)
        if series==:first
                return first_series(mmin, mmax)
        elseif series==:principal
                return principal_series(mmin, mmax)
        elseif series==:hull
                return hull_series(mmin,mmax)
        elseif series==:bulkhead
                return bulkhead_series(mmin,mmax)
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

function path_points(path_list::Vector{Vector{Vector{Float64}}}, N::Int64)::(Tuple{Array{Vector{Float64},2},Array{Array{Int64,1},1},Array{Array{Float64,1},1}})
    paths           = Vector{Vector{Float64}}[];
    odometers       = Vector{Float64}[];
    corner_indices  = [zeros(Int64, length(path_list[i])) for i=1:length(path_list)];
    for i ∈ 1:length(path_list)
        positions = path_points(path_list[i], N)
        push!(paths,    positions)
        push!(odometers, vcat([0.0],cumsum(LinearAlgebra.norm.(diff(positions)))))
        #println(indexin.(hash.(path_list[i]), Ref(hash.(positions))))
        corner_indices[i] .= findindex.(Ref(positions), path_list[i])
    end

    hcat(paths...), corner_indices, odometers
end

function path_points(bz_pt_ledger::Dict{String,Array{Float64,1}}, path_list, N::Int64)::(Tuple{Array{Vector{Float64},2},Array{Array{Int64,1},1},Array{Array{Float64,1},1}})
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
function bands(RN, asd, mn::Tuple{Int,Int}, pathlist=["K1","Γ","M1","K'1"], Npath = 300, scriptdir=ENV["scriptdir"], cachedirr=ENV["cachedir"], args...)
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
#### Calculating DataIntegrals
######################################

function integral(RN, asd, mn, dtype::Type{T} where T <: SolidState.DataChart, indices, priors, base, Neval, pool=default_worker_pool(), cachedir=ENV["cachedir"], scriptdir=ENV["scriptdir"])

    println("Calculating Integral for $asd $(mn[1])-$(mn[2])");flush(stdout)
    models(asd,[mn])
    di = DataIntegral(asd, mn, dtype, indices, priors, base)

    di(Neval)

    bson("$(mkpath("$scriptdir/out/$RN"))/$asd-$(mn[1])-$(mn[2]).bson", Dict(
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
        :handle => "$asd-$(mn[1])-$(mn[2])",
        :data => di
    ))
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
                :comargs => models(asd,twist_series(:bulkhead, comargs...)),
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

function load(RN,name)
        BSON.load("$(ENV["scriptdir"])/out/$RN/$name.bson")
end

function extract(RN,name,plotfunction,args...)
    datafile = "$(ENV["scriptdir"])/out/$RN/$name.bson"
    plotdir = "$(ENV["scriptdir"])/plot/$RN"|>mkpath
    plotfunction(args...; BSON.load(datafile)..., plotdir=plotdir)
end

function b2_extract(RN,name,plotfunction,args...)
    plotdir = "$(ENV["scriptdir"])/plot/$RN"|>mkpath
    source = "$(ENV["b2scriptdir"])/out/$RN/$name.bson"
    target = "$(mkpath("$(ENV["scriptdir"])/out/b2-$RN"))/$name.bson"
    run(`rsync -r --progress asmithc@bridges2.psc.edu:$source $target`)
    plotfunction(args...; BSON.load(target)..., plotdir=plotdir)
end

end
