module Jobs

using BSON

const JULIA_CALL="julia -O 3"
const MACHINEFILE="--machine-file <(srun hostname -s)"
const JPROJECT="--project=$(ENV["HOME"])/.julia/environments/v1.6/Project.toml"

function julia_call(jit,big,job,JID,scriptdir=ENV["scriptdir"])
    if jit
        if big===true
            "      $JULIA_CALL $JPROJECT $MACHINEFILE $scriptdir/bin/$job/run-$JID/$job-$JID.jl"
        else
            "      $JULIA_CALL $JPROJECT -p $big $scriptdir/bin/$job/run-$JID/$job-$JID.jl"
        end
    else
        if big===true
            "      $JULIA_CALL $MACHINEFILE $JPROJECT --sysimage=$(ENV["scriptdir"])/.cache/system/sysimage.dylib $scriptdir/bin/$job/run-$JID/$job-$JID.jl"
        else
            "      $JULIA_CALL $JPROJECT --sysimage=$(ENV["scriptdir"])/.cache/system/sysimage.dylib -p $big $scriptdir/bin/$job/run-$JID/$job-$JID.jl"
        end
    end
end

function bash_run_file(RN,jit,big,job,JID,scriptdir = ENV["scriptdir"])
    open("$scriptdir/bin/$job/run-$JID/$job-$JID.sh",create=true,write=true) do io
        write(io,
"#!/bin/bash

echo \"************************************************************\" >> $scriptdir/bin/$job/run-$JID/$job-$JID.o
echo \"run on \$(date) in \$PWD\"  >> $scriptdir/bin/$job/run-$JID/$job-$JID.o
echo '-host $(ENV["HOME"]) -j $job -p $big'  >> $scriptdir/bin/$job/run-$JID/$job-$JID.o

{
    time {
            {
                $(julia_call(jit,big,job,JID,scriptdir))
            } >> $scriptdir/bin/$job/run-$JID/$job-$JID.o
    }
} 2>> $scriptdir/bin/$job/run-$JID/$job-$JID.o
"
        )
    end
end

function julia_run_file(job,jobargs,JID,scriptdir = ENV["scriptdir"],cachedir = ENV["cachedir"])
    open("$scriptdir/bin/$job/run-$JID/$job-$JID.jl",create=true,write=true) do io
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
    "bash $scriptdir/bin/$job/run-$JID/$job-$JID.sh"
end

function queue(RN,jit,big,job,jobargs,slurmargs,scriptdir = ENV["scriptdir"],cachedir=ENV["cachedir"])
    JID = prep_job(RN,jit,big,job,jobargs,scriptdir,cachedir)
    "sbatch -o $scriptdir/bin/$job/run-$JID/$job-$JID.o --mail-type=ALL $slurmargs $scriptdir/bin/$job/run-$JID/$job-$JID.sh"
end

using SolidState
function spray(RN,asd,idxS,idxL,job,chartinfo,neval,f,jit,ps,N,npn,pray)

    slurmargs = "-p RM -N $N --ntasks-per-node=$npn -t 48:00:00"

    comargs = BSON.load("$(@__DIR__)/mns.bson")[:mns]
    cmds = map(comargs[idxS:idxL]) do mn
        jobargs = (RN,asd,mn,(chartinfo...,neval))
        f(RN,jit,ps,job,jobargs,slurmargs)
    end
    mkpath("$(pwd())/$(split(RN,"/")[1])")
    # open("$(ENV["scriptdir"])/bin/$RN.sh",create=true,write=true) do io
    open("$(pwd())/$RN.sh",create=true,write=true) do io
        write(io,"#!/bin/bash \n")
        map(cmds) do cmd
            write(io,"$cmd \n")
        end
    end

    pray ? Base.run(`bash $(pwd())/$RN.sh`) : `bash $(pwd())/$RN.sh`
end

function models(f,jit,ps,asds,(idxS,idxL),ncpus,pray)
    RN = "models"
    comargs = BSON.load("$(@__DIR__)/mns.bson")[:mns]

    cmds = map(asds) do asd
        f(RN,false,true,"models",(asd,comargs[idxS:idxL]),"-p RM-shared -n $ncpus -t 10:00:00")
    end
    # open("$(ENV["scriptdir"])/bin/$RN.sh",create=true,write=true) do io
    open("$(pwd())/$models.sh",create=true,write=true) do io
        write(io,"#!/bin/bash \n")
        map(cmds) do cmd
            write(io,"$cmd \n")
        end
    end

    pray ? Base.run(`bash $(pwd())/$RN.sh`) : `bash $(pwd())/$RN.sh`

end

########################################################################################
#### Retrieving Data

function pull_b2_data(RN)
    Base.run(`rsync -r --progress asmithc@bridges2.psc.edu:/ocean/projects/phy190028p/asmithc/scripts/out/$RN $(ENV["scriptdir"])/out/`)
end

function aggregate_data(RN,force=false)
    if !isfile("$(ENV["scriptdir"])/out/$RN/$RN.bson")
        files = readdir("$(ENV["scriptdir"])/out/$RN", join=true)
        data = BSON.load.(files)
        bson("$(ENV["scriptdir"])/out/$RN/$RN.bson",Dict(:data=>data))
    elseif isfile("$(ENV["scriptdir"])/out/$RN/$RN.bson")&&force
        Base.rm("$(ENV["scriptdir"])/out/$RN/$RN.bson")
        files = readdir("$(ENV["scriptdir"])/out/$RN", join=true)
        data = BSON.load.(files)
        bson("$(ENV["scriptdir"])/out/$RN/$RN.bson",Dict(:data=>data))
    else
        println("already aggregated")
    end
end

function datapull(RN)
    Jobs.pull_b2_data(RN)
    Jobs.aggregate_data(RN, true)
end

end
