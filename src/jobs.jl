module Jobs

using BSON

const JULIA_CALL="julia -O 3"
const MACHINEFILE="--machine-file <(srun hostname -s)"
const JPROJECT="--project=$(ENV["HOME"])/.julia/environments/v1.6/Project.toml"

function julia_call(RN,name,jit,big,job,JID,scriptdir=ENV["scriptdir"])
    if jit
        if big===true
            "      $JULIA_CALL $JPROJECT $MACHINEFILE $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.jl"
        else
            "      $JULIA_CALL $JPROJECT -p $big $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.jl"
        end
    else
        if big===true
            "      $JULIA_CALL $MACHINEFILE $JPROJECT --sysimage=$(ENV["scriptdir"])/.cache/system/sysimage.dylib $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.jl"
        else
            "      $JULIA_CALL $JPROJECT --sysimage=$(ENV["scriptdir"])/.cache/system/sysimage.dylib -p $big $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.jl"
        end
    end
end

function bash_run_file(RN,name,jit,big,job,JID,scriptdir = ENV["scriptdir"])
    open("$(pwd())/$RN/$name/bin/run-$JID/$JID-$name.sh",create=true,write=true) do io
        write(io,
"#!/bin/bash

echo \"************************************************************\" >> $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.o
echo \"run on \$(date) in \$PWD\"  >> $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.o
echo '-host $(ENV["HOME"]) -j $job -p $big'  >> $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.o

{
    time {
            {
                $(julia_call(RN,name,jit,big,job,JID,scriptdir))
            } >> $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.o
    }
} 2>> $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.o
"
        )
    end
end

function julia_run_file(RN,name,job,jobargs,JID,scriptdir = ENV["scriptdir"],cachedir = ENV["cachedir"])
    open("$(pwd())/$RN/$name/bin/run-$JID/$JID-$name.jl",create=true,write=true) do io
        write(io,
"using Distributed
@everywhere begin
    ENV[\"cachedir\"]=\"$cachedir\";
    ENV[\"scriptdir\"]=\"$scriptdir\";
    using SolidState, OrderedCollections, LinearAlgebra, BSON, Measures
end
SolidState.Main.$job$jobargs
"
        )
    end
end

function prep_job(RN,name,jit,big,job,jobargs,scriptdir=ENV["scriptdir"],cachedir=ENV["cachedir"])
    # Counts the number of jobs already executed
    mkpath("$(pwd())/$RN/$name/bin")
    JID = Int(length(readdir("$(pwd())/$RN/$name/bin"))+1)
    # Makes Directory for new job
    mkpath("$(pwd())/$RN/$name/bin/run-$JID")
    # Makes julia file for new job
    julia_run_file(RN,name,job,jobargs,JID,scriptdir,cachedir)
    # Makes shell file for new job
    bash_run_file(RN,name,jit,big,job,JID,scriptdir)
    JID
end

function run(RN,name,jit,big,job,jobargs,slurmargs,scriptdir = ENV["scriptdir"],cachedir=ENV["cachedir"])
    JID = prep_job(RN,name,jit,big,job,jobargs,scriptdir,cachedir)
    # run(
    "bash $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.sh"
end

function queue(RN,name,jit,big,job,jobargs,slurmargs,scriptdir = ENV["scriptdir"],cachedir=ENV["cachedir"])
    JID = prep_job(RN,name,jit,big,job,jobargs,scriptdir,cachedir)
    "sbatch -o $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.o --mail-type=ALL $slurmargs $(pwd())/$RN/$name/bin/run-$JID/$JID-$name.sh"
end

function models(RN,asds,(idxS,idxL),f,jit,ps,slurmargs,pray)
    comargs = BSON.load("$(@__DIR__)/mns.bson")[:mns]

    cmds = map(asds) do asd
        f(RN,"models",false,true,"models",(asd,comargs[idxS:idxL]),slurmargs)
    end
    # open("$(ENV["scriptdir"])/bin/$RN.sh",create=true,write=true) do io
    open("$(pwd())/$RN/models/$asd.sh",create=true,write=true) do io
        write(io,"#!/bin/bash \n")
        map(cmds) do cmd
            write(io,"$cmd \n")
        end
    end

    pray ? Base.run(`bash $(pwd())/$RN/models/$asd.sh`) : `bash $(pwd())/$RN/models/$asd.sh`

end

function spray(RN::String,name::String,asd::Function,(idxS,idxL),job,jobargs,f,jit,ps,slurmargs,pray)

    comargs = BSON.load("$(@__DIR__)/mns.bson")[:mns]
    cmds = map(comargs[idxS:idxL]) do mn
        fulljobargs = ("$RN/$name/$asd",asd,mn,jobargs...)
        f(RN,name,jit,ps,job,fulljobargs,slurmargs)
    end
    # mkpath("$(pwd())$(string(("/".*split(RN,"/")[1:end-1])...))")|>println
    mkpath("$(pwd())/$RN/$name")
    open("$(pwd())/$RN/$name/$asd.sh",create=true,write=true) do io
        write(io,"#!/bin/bash \n")
        map(cmds) do cmd
            write(io,"$cmd \n")
        end
    end

    pray ? Base.run(`bash $(pwd())/$RN/$name/$asd.sh`) : `bash $(pwd())/$RN/$name/$asd.sh`

end

function spray(RN::String,name::String,asds::AbstractVector,rng,job,jobargs,runargs)
    map(asds) do asd
        spray(RN,name,asd,rng,job,jobargs,runargs...)
    end
end

function run_env_setup(p=1)
    if ENV["HOME"]=="/jet/home/asmithc"
        ENV["scriptdir"] = "/ocean/projects/phy190028p/asmithc/scripts"
        ENV["cachedir"]  = "/ocean/projects/phy190028p/asmithc/scripts/.cache"
        (queue,false, true)
    else
        ENV["scriptdir"] = "$(ENV["HOME"])/Dropbox/Graduate/scripts"
        ENV["cachedir"] = "$(ENV["HOME"])/Dropbox/Graduate/scripts/.cache"
        (run,  true, p)
    end
end

function spray(RN::String,asds::AbstractVector,jobs::AbstractDict,pray::Bool;p=1)
    cmds = vcat(map(keys(jobs),values(jobs)) do name,(rng,job,jobargs,slurmargs)
        runargs = (run_env_setup(p)...,slurmargs,false)
        spray(RN,name,asds,rng,job,jobargs,runargs)
    end...)

    open("$(pwd())/$RN/$RN.sh",create=true,write=true) do io
        write(io,"#!/bin/bash \n")
        map(cmds) do cmd
            write(io,"$(string(cmd)[2:end-1]) \n")
        end
    end

    pray ? Base.run(`bash $(pwd())/$RN/$RN.sh`) : `bash $(pwd())/$RN/$RN.sh`

end

########################################################################################
#### Utilities
using Pkg;

function update_package()

    if ENV["HOME"]=="/jet/home/asmithc"
        Base.rm("$(ENV["HOME"])/.julia/dev/DraftMill",force=true,recursive=true)
        Base.rm("$(ENV["HOME"])/.julia/dev/SolidState",force=true,recursive=true)
        Pkg.rm("SolidState")
        Pkg.rm("DraftMill")
        Pkg.develop(url="https://gitlab.com/solidstateapps/SolidState")
        Pkg.develop(url="https://gitlab.com/smith-and/PaperMill.jl")
    else
        dir = pwd()

        cd("$(ENV["HOME"])/Dropbox/Graduate/dev/SolidState")
        Base.run(`git add .`)
        Base.run(`git status`)
        try Base.run(`git commit -m 'update'`); catch end
        Base.run(`git push`)

        cd("$(ENV["HOME"])/Dropbox/Graduate/dev/DraftMill.jl")
        Base.run(`git add .`)
        Base.run(`git status`)
        try Base.run(`git commit -m 'update'`); catch end
        Base.run(`git push`)

        cd(dir)
    end
end

function qcheck()
    Base.run(`squeue -u $(ENV["USER"])`)
end

function glance(job,r)
    Base.run(`cat $(ENV["scriptdir"])/$job/bin/run-$r/$job-$r.o`)
end

function check(job,r)
    Base.run(`cat $(ENV["scriptdir"])/$job/bin/run-$r/$job-$r.jl`)
end

function mount()
    mkpath("$(ENV["HOME"])/Dropbox/Graduate/scripts/b2scripts")
    Base.run(`sshfs asmithc@bridges2.psc.edu:/ocean/projects/phy190028p/asmithc/scripts $(ENV["HOME"])/Dropbox/Graduate/scripts/b2scripts`)
end

function umount()
    Base.run(`umount $(ENV["HOME"])/Dropbox/Graduate/scripts/b2scripts`)
    rm("$(ENV["HOME"])/Dropbox/Graduate/scripts/b2scripts",recursive=true)
end

########################################################################################
#### Retrieving Data

function pull_b2_data(RN)
    Base.run(`rsync -r --progress asmithc@bridges2.psc.edu:/ocean/projects/phy190028p/asmithc/scripts/out/$RN $(ENV["scriptdir"])/out/`)
end

function aggregate_data(dir,RN,force=false)
    if !isfile("$(ENV["scriptdir"])/out/$dir/$RN/$RN.bson")
        files = readdir("$(ENV["scriptdir"])/out/$dir/$RN", join=true)
        data = BSON.load.(files)
        bson("$(ENV["scriptdir"])/out/$dir/$RN.bson",Dict(:data=>data))
    elseif isfile("$(ENV["scriptdir"])/out/$dir/$RN/$RN.bson")&&force
        Base.rm("$(ENV["scriptdir"])/out/$dir/$RN/$RN.bson")
        files = readdir("$(ENV["scriptdir"])/out/$dir/$RN", join=true)
        data = BSON.load.(files)
        bson("$(ENV["scriptdir"])/out/$dir/$RN.bson",Dict(:data=>data))
    else
        println("already aggregated")
    end
end


function aggregate_project(dir,force=false)
    map(readdir("$(ENV["scriptdir"])/out/$dir")) do RN
        aggregate_data(dir,RN,force)
    end
end

function datapull(RN)
    Jobs.pull_b2_data(RN)
    Jobs.aggregate_data(RN, true)
end

end
