function dir_extract(dir, asd, RN; local_dir, mount_dir, kargs...)
    #Extract these directories from bridges
    dir_contents = readdir(mkpath("$local_dir/$asd/$dir/$RN"))
    ext_dirs = dir_contents[isdir.("$local_dir/$asd/$dir/$RN/" .* dir_contents)]

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
        println("-----Extracting: $local_dir/$asd/$dir/$RN")
        run(`rsync -r --progress asmithc@bridges2.psc.edu:$mount_dir/$asd/$dir/$RN $local_dir/$asd/$dir/`)
    else
        println("-----Already Extracted: $local_dir/$asd/$dir/$RN")
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
    import_charts(asd, RN; local_dir, kargs...)::AbstractChannel
"""
function import_charts(asd, RN; local_dir, kargs...)::AbstractChannel
    #Accumulate Information from the local directories
    chnl = Channel{Dict{Symbol,Any}}(100)
    for (root, dirs, files) in walkdir("$local_dir/$asd/.out/$RN")
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
    export_collection(chnl::AbstractChannel, asd, RN; local_dir, force = false, kargs...)
"""
function export_collection(chnl::AbstractChannel, asd, RN; local_dir, force = false, kargs...)
    println("-----Starting Collection Export")
    extraction_dir = mkpath("$local_dir/$asd/.out/$RN");
    if !isfile("$extraction_dir/$RN.bson") || force
        println("-----Exporting Collection")
        bson("$extraction_dir/$RN.bson",Dict(:chnl=>chnl))
    else
        println("-----Already Exported Collection")
    end
end

function compose_collection(asd, RN; local_dir, mount_dir, force=false, kargs...)
    #Copying from the mount
    dir_extract(".plot", asd, RN; local_dir=local_dir, mount_dir=mount_dir)
    dir_extract(".out",  asd, RN; local_dir=local_dir, mount_dir=mount_dir)

    #Loading BSON Files & Composing Data
    chnl = import_charts(asd, RN; local_dir=local_dir)

    #Writing Collection to BSON
    export_collection(chnl, asd, RN, local_dir=local_dir, force=force)

    println("-----Done Composing")
end

"""
    import_collection(asd, RN; local_dir, mount_dir, force=false, kargs...)

Import collection data and make new collection if there is no collection.
  - asd: Symbol for the ASD
  - RN: Symbol for the run name
  - mount_dir: path of the directory of the sftp mount
  - local_dir: path of the directory to look for an extraction (and make if none)
"""
function import_collection(asd, RN; local_dir, mount_dir, force=false, kargs...)
    collection_path = "$local_dir/$asd/.out/$RN/$RN.bson"
    dir_extract(".plot", asd, RN; local_dir=local_dir, mount_dir=mount_dir)
    dir_extract(".out",  asd, RN; local_dir=local_dir, mount_dir=mount_dir)

    if isfile(collection_path)
        println("-----Importing Collection:")
        println(collection_path)
        return BSON.load(collection_path)[:chnl]
    else
        println("-----No Collection")
        if typeof(mount_dir)==String
            println("-----Composing Collection")
            compose_collection(asd, RN; local_dir=local_dir, mount_dir=mount_dir, force=force)
            BSON.load(collection_path)[:chnl]
        end
    end
end

"""
    import_collection(; asd, RN, local_dir, mount_dir, force=false, kargs...)

A keyword argument only interface for the import
"""
function import_collection(; asd, RN, local_dir, mount_dir, force=false, kargs...)
    import_collection(asd, RN; local_dir=local_dir, mount_dir=mount_dir, force=force, kargs...)
end
