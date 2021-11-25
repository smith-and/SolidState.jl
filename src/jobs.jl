module Jobs

using BSON

const JULIA_CALL="julia -O 3"
const MACHINEFILE="--machine-file <(srun hostname -s)"
const JPROJECT="--project=$(ENV["HOME"])/.julia/environments/v1.6/Project.toml"

function julia_call(jit,big,job,JID,scriptdir=ENV["scriptdir"])
    if jit
        if big===true
            "      $JULIA_CALL $JPROJECT $MACHINEFILE $scriptdir/bin/$job/run-$JID/$JID-$job.jl"
        else
            "      $JULIA_CALL $JPROJECT -p $big $scriptdir/bin/$job/run-$JID/$JID-$job.jl"
        end
    else
        if big===true
            "      $JULIA_CALL $MACHINEFILE $JPROJECT --sysimage=$(ENV["scriptdir"])/.cache/system/sysimage.dylib $scriptdir/bin/$job/run-$JID/$JID-$job.jl"
        else
            "      $JULIA_CALL $JPROJECT --sysimage=$(ENV["scriptdir"])/.cache/system/sysimage.dylib -p $big $scriptdir/bin/$job/run-$JID/$JID-$job.jl"
        end
    end
end

function bash_run_file(RN,jit,big,job,JID,scriptdir = ENV["scriptdir"])
    open("$scriptdir/bin/$job/run-$JID/$JID-$job.sh",create=true,write=true) do io
        write(io,
"#!/bin/bash

echo \"************************************************************\" >> $scriptdir/bin/$job/run-$JID/$JID-$job.o
echo \"run on \$(date) in \$PWD\"  >> $scriptdir/bin/$job/run-$JID/$JID-$job.o
echo '-host $(ENV["HOME"]) -j $job -p $big'  >> $scriptdir/bin/$job/run-$JID/$JID-$job.o

{
    time {
            {
                $(julia_call(jit,big,job,JID,scriptdir))
            } >> $scriptdir/bin/$job/run-$JID/$JID-$job.o
    }
} 2>> $scriptdir/bin/$job/run-$JID/$JID-$job.o
"
        )
    end
end

function julia_run_file(job,jobargs,JID,scriptdir = ENV["scriptdir"],cachedir = ENV["cachedir"])
    open("$scriptdir/bin/$job/run-$JID/$JID-$job.jl",create=true,write=true) do io
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

function prep_job(RN,jit,big,job,jobargs,scriptdir=ENV["scriptdir"],cachedir=ENV["cachedir"])
    # Counts the number of jobs already executed
    mkpath("$scriptdir/bin/$job")
    JID = Int(length(readdir("$scriptdir/bin/$job"))+1)
    # Makes Directory for new job
    mkpath("$scriptdir/bin/$job/run-$JID")
    # Makes julia file for new job
    julia_run_file(job,jobargs,JID,scriptdir,cachedir)
    # Makes shell file for new job
    bash_run_file(RN,jit,big,job,JID,scriptdir)
    JID
end

function run(RN,jit,big,job,jobargs,slurmargs,scriptdir = ENV["scriptdir"],cachedir=ENV["cachedir"])
    JID = prep_job(RN,jit,big,job,jobargs,scriptdir,cachedir)
    # run(
    "bash $scriptdir/bin/$job/run-$JID/$JID-$job.sh"
end

function queue(RN,jit,big,job,jobargs,slurmargs,scriptdir = ENV["scriptdir"],cachedir=ENV["cachedir"])
    JID = prep_job(RN,jit,big,job,jobargs,scriptdir,cachedir)
    "sbatch -o $scriptdir/bin/$job/run-$JID/$JID-$job.o --mail-type=ALL $slurmargs $scriptdir/bin/$job/run-$JID/$JID-$job.sh"
end

using SolidState
function spray(RN,asd,idxS,idxL,job,jobargs,f,jit,ps,slurmargs,pray)

    comargs = BSON.load("$(@__DIR__)/mns.bson")[:mns]
    cmds = map(comargs[idxS:idxL]) do mn
        fulljobargs = (RN,asd,mn,jobargs...)
        f(RN,jit,ps,job,fulljobargs,slurmargs)
    end
    mkpath("$(pwd())$(string(("/".*split(RN,"/")[1:end-1])...))")|>println
    open("$(pwd())/$RN.sh",create=true,write=true) do io
        write(io,"#!/bin/bash \n")
        map(cmds) do cmd
            write(io,"$cmd \n")
        end
    end

    pray ? Base.run(`bash $(pwd())/$RN.sh`) : `bash $(pwd())/$RN.sh`

end

function models(f,jit,ps,asds,(idxS,idxL),slurmargs,pray)
    RN = "models"
    comargs = BSON.load("$(@__DIR__)/mns.bson")[:mns]

    cmds = map(asds) do asd
        f(RN,false,true,"models",(asd,comargs[idxS:idxL]),slurmargs)
    end
    # open("$(ENV["scriptdir"])/bin/$RN.sh",create=true,write=true) do io
    open("$(pwd())/$RN.sh",create=true,write=true) do io
        write(io,"#!/bin/bash \n")
        map(cmds) do cmd
            write(io,"$cmd \n")
        end
    end

    pray ? Base.run(pipeline(`bash $(pwd())/$RN.sh`,"$(pwd())/$RN.txt")) : pipeline(`bash $(pwd())/$RN.sh`,"$(pwd())/$RN.txt")

end

function qcheck()
    Base.run(`squeue -u $(ENV["USER"])`)
end

function glance(job,r)
    Base.run(`cat $(ENV["scriptdir"])/bin/$job/run-$r/$job-$r.o`)
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
        bson("$(ENV["scriptdir"])/out/$dirks/$RN.bson",Dict(:data=>data))
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
