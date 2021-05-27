module Extraction

using Revise
using Distributed, Dates, OrderedCollections, BSON, Plots
using LinearAlgebra, SharedArrays, StaticArrays
using Mux, WebIO, Interact, InteractiveUtils
using CubicSplines, Roots, SpecialFunctions, HCubature
using Base: Threads
using ..SolidState

############################
### Extraction
############################
function dir_extract(dir, asd, RN; target_dir, mount_dir, kargs...)
    #Extract these directories from bridges
    dir_contents = readdir(mkpath("$target_dir/$asd/$dir/$RN"))
    ext_dirs = dir_contents[isdir.("$target_dir/$asd/$dir/$RN/" .* dir_contents)]

    extract  = false
    if length(ext_dirs)==0
        extract = true
    elseif length(ext_dirs)==1
        if ext_dirs[1]=="Ext"
            extract=true
        end
    elseif length(ext_dirs) >= 2
        extract=false
    end

    if extract
        println("-----Extracting: $target_dir/$asd/$dir/$RN")
        run(`rsync -r --progress asmithc@bridges2.psc.edu:$mount_dir/$asd/$dir/$RN $target_dir/$asd/$dir/`)
    else
        println("-----Already Extracted: $target_dir/$asd/$dir/$RN")
    end
end

function datapull(dict::AbstractDict)
    return dict
end

function datapull(data::DataIntegral)
    (
        data  = data.data,
        chart = data.dm.chart,
        evals = data.evals,
        error = data.err
    )
end

"""
    import_charts(asd, RN; target_dir, kargs...)::AbstractChannel
"""
function import_charts(asd, RN; target_dir, kargs...)::AbstractChannel
    #Accumulate Information from the local directories
    chnl = Channel{Dict{Symbol,Any}}(100)
    for (root, dirs, files) in walkdir("$target_dir/$asd/.out/$RN")
        for file in files
            if occursin(".bson",file) && !occursin("error",file)
                println("Added file: ")
                println(joinpath(root, file)) # path to files
                model = root[(findlast("/",root)[1]+1):end]
                idxs = findall(isequal('-'),model)

                asd = Symbol(model[1:(idxs[1]-1)])
                mn   = (parse(Int, model[(idxs[1]+1):(idxs[2]-1)]),parse(Int, model[(idxs[2]+1):end]))

                ASD = CommensurateASD(eval(Expr(:call,asd)),mn)
                Λ   = ASD["blv"]
                vol = LinearAlgebra.det(Λ)

                println("Loading")
                data_dict = Dict(
                    :filename => "$file",
                    :datafile => ("$root/$file"),
                    :asd   => asd,
                    :mn    => mn,
                    :angle => cθ(mn...)*180/π,
                    :ucvol => vol,
                )

                try
                    di = data_import("$root/$file")
                    push!(data_dict, :data  => datapull(di))
                catch
                    di = BSON.load("$root/$file")
                    push!(data_dict, :data  => datapull(di))
                end

                put!(chnl, data_dict)
            end
        end
    end
    println("-----Sorting by angle")
    #Sort the channel data by twist angle
    sort!(chnl.data,by=(x->x[:angle]))
    println("-----Done Loading")
    chnl
end

"""
    export_collection(chnl::AbstractChannel, asd, RN; target_dir, force = false, kargs...)
"""
function export_collection(chnl::AbstractChannel, asd, RN; target_dir, force = false, kargs...)
    println("-----Starting Collection Export")
    target_dir = mkpath("$target_dir/$asd/.out/$RN");
    if !isfile("$target_dir/$RN.bson") || force
        println("-----Exporting Collection")
        bson("$target_dir/$RN.bson",Dict(:chnl=>chnl))
    else
        println("-----Already Exported Collection")
    end
end

function compose_collection(asd, RN; target_dir, mount_dir, force=false, kargs...)
    #Copying from the mount
    dir_extract(".plot", asd, RN; target_dir=target_dir, mount_dir=mount_dir)
    dir_extract(".out",  asd, RN; target_dir=target_dir, mount_dir=mount_dir)

    #Loading BSON Files & Composing Data
    chnl = import_charts(asd, RN; target_dir=target_dir)

    #Writing Collection to BSON
    export_collection(chnl, asd, RN, target_dir=target_dir, force=force)

    println("-----Done Composing")
end

"""
    import_collection(asd, RN; target_dir, mount_dir, force=false, kargs...)

Import collection data and make new collection if there is no collection.
  - asd: Symbol for the ASD
  - RN: Symbol for the run name
  - mount_dir: path of the directory of the sftp mount
  - target_dir: path of the directory to look for an extraction (and make if none)
"""
function load(asd, RN; datadir=ENV["datadir"], mountdir=ENV["mountdir"], force=false, kargs...)
    collection_path = "$datadir/$asd/.out/$RN/$RN.bson"
    dir_extract(".plot", asd, RN; target_dir=datadir, mount_dir=mountdir)
    dir_extract(".out",  asd, RN; target_dir=datadir, mount_dir=mountdir)

    if isfile(collection_path)
        println("-----Importing Collection:")
        println(collection_path)
        return BSON.load(collection_path)[:chnl]
    else
        println("-----No Collection")
        if typeof(mount_dir)==String
            println("-----Composing Collection")
            compose_collection(asd, RN; target_dir=datadir, mount_dir=mountdir, force=force)
            BSON.load(collection_path)[:chnl]
        end
    end
end


function load(RN,asd,mn)
    target = "$(ENV["scriptdir"])/out/$RN/$asd-$(mn[1])-$(mn[2]).bson"

    try
        data_import(target)
    catch
        BSON.load(target)
    end
end

"""
    dataload(description)

data4 = Dict(
    :asd   => ASD4,
    :shg   => "SHG-1-50k",
    :bands => "BANDS",
)|>dataload
"""
function dataload(description)
    Dict(map(keys(description),values(description)) do name,tag
        if name != :asd
            name=>Extraction.load(description[:asd],tag)
        else
            name=>tag
        end
    end...)
end

end
