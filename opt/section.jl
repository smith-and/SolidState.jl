function prestack(batchdir::String=pwd())
    println("compiling a system image")
    mkpath(batchdir*"/system")
    stats = @timed PackageCompiler.create_sysimage([:GSK,:TensorCharts,:Lattice,:Model,:Analysis,:Program,:FlakeStacks],
        sysimage_path               = batchdir*"/system/out/sysimage.dylib",
        precompile_execution_file   = batchdir*"/system/bin/startup.jl",
        #project                     = "$(ENV["HOME"])/.julia/environments/v1.5/Project.toml",
        incremental                 = true
    )
    println("build time \t"*string(stats.time/60)*"m")
    nothing
end


#Expands a scope (named tuple of collectable objects e.g. ranges, arrays) into the product span of its ranges (Array of NamedTuples)
function scope_span(scope)
    raw_span = collect(Iterators.product(collect.(collect(values(scope)))...))
    tupletype= NamedTuple{Tuple(collect(keys(scope))),eltype(raw_span)}
    tupletype.(raw_span)
end

function programargs(;
        asds=[:BNAB], mrange=1:1, srange=0:0, shifts=[:none], layers=[1], comargs=:none, mmax=2,
        paths = [
            ["K'1", "Γ", "K1", "M1"],
            ["K'2", "Γ", "K2", "M2"],
            ["K'3", "Γ", "K3", "M3"]
            ],
        mesh = [["K1","K2"],[0.0,0.0,0.0]],
        Σinfo::Tuple{Symbol,Float64,Float64,Int64} = (:δ,0.02,0.02,1),
        priors = [(:T, 0.0,0.0,1),(:ν,0.5,0.5,1)],
        dos  = (10,), jdos = (10,), ldos = (10,), shdos= (10,),
        shg  = (10, [(1,1,1),(2,2,2),(2,1,2),(1,2,1)]),
        cpge = (10, [(1,1,1),(2,2,2),(2,1,2),(1,2,1)]),
        optσ = (10, [(1,1),(2,2),(2,1),(1,2)]),
        quadrange=10:10:10, pathrange=50:10:50, meshrange=10:10:10,
        job="warmup",
        exportpath = (pwd()*"/"*job*"/out"),
        commpath = (pwd()*"/jobs/"*job),
        batchdir   = pwd(),
        timestamp  = Dates.format(now(),"/T-yyyy-mm-dd-HH-MM-SS"),
        batch=true, reporting=false, serving=false, presenting=false,
        saving=true, plotting=true, sections=sections,
        kargs...
    )

    if comargs!=:none
        model_span = vcat([scope_span(  (shifts=shifts, layers=layers, asd=asds, m=comargs[i][1], s=comargs[i][2]))[:] for i∈1:length(comargs)]...)
    else
        model_span = scope_span((shifts=shifts, layers=layers, asd=asds, m=mrange, s=srange))[:]
    end

    return (
        m= model_span,
        a=(psN=pathrange, msN=meshrange, qN=quadrange),
        ex=(
            spaths=paths, smeshes=mesh, priors=priors,
            #The Chart Section Information
            spectral=(
                    ("dos" , Complex{Float64},   [(:ω,-10.0,10.0,dos[1]::Int),  Σinfo], [(1,)], 1),
                ),
            moredos =(
                    ("jdos" , Complex{Float64},  [ (:ω,-10.0,10.0,jdos[1]::Int), Σinfo], [(1,)], 1),
                    ("ldos" , Complex{Float64},  [ (:ω,-10.0,10.0,ldos[1]::Int), Σinfo], [(1,)], 1),
                    ("shdos" ,Complex{Float64},  [(:ω,-10.0,10.0,shdos[1]::Int), Σinfo], [(1,)], 2)
            ),
            optical=(
                    ("shg"  , Complex{Float64},   [ (:ω,0.0,10.0,shg[1]::Int), Σinfo], shg[2] , 2),
                    #("cpge" , Complex{Float64},   [(:ω,0.0,10.0,cpge[1]::Int), Σinfo], cpge[2], 1),
                    ("optcond", Complex{Float64}, [(:ω,0.0,10.0,optσ[1]::Int), Σinfo], optσ[2], 1)
            ),
            expath=exportpath, commpath = commpath,
            batch=batch, reporting=reporting, serving=serving,
            presenting=presenting, plotting=plotting, saving=saving, timestamp=timestamp,
            batchdir=batchdir, job=job, sections=sections
        )
    )
end

const analysis_map = Dict{Symbol,DataType}(
    :dos     => Analysis.DosData{Int64,Float64},
    :mdos    => Analysis.MoreDosData{Int64,Float64},
    :optical => Analysis.OpticalData{Int64,Float64}
)

function jobargs(   name="warmup"; batchdir::String=pwd(), timestamp = Dates.format(now(),"/T-yyyy-mm-dd-HH-MM-SS"),
                    plotting=true, saving=true,
                    asds=[:BNAB], ms=1:1, ss=1:1,
                    layers=[1], shifts=[:none],
                    series=:none, mmax=2, drng= (0,10),
                    sc=2, dsc=3,
                    reginfo=(0.02,0.02,1),
                    qsc=2,  qstep=10, qlength=1,
                    psc=10, pstep=10, plength=1,
                    msc=10, mstep=10, mlength=1,
                    paths= [["Q'1","Q1","Γ","Q2","Q'1"],
                            ["K'2", "M2", "Γ","M'3","K1"],
                            ["K'2","M2","K2","Γ","K'3","M'3","K1"],
                            ],
                    mesh      = [["K1","K2"],[0.0,0.0,0.0]],
                    priors    = [(:T, 0.0,0.0,1),(:ν,0.5,0.5,1)],
                    Σinfo     = (:δ, reginfo...),
                    dos       = (dsc,),
                    jdos      = (sc,),
                    ldos      = (sc,),
                    shdos     = (sc,),
                    shg       = (sc, [(2,2,2)]),
                    cpge      = (2,    [(3,2,1)]),
                    optσ      = (sc, [(1,1),(1,2)]),
                    sections  = Symbol[],
                    args...)
    section_types = getindex.(Ref(analysis_map),sections)

    pargs = (
            quadrange=(qsc):qstep:(qsc+qstep*(qlength-1)),
            pathrange=(psc):pstep:(psc+pstep*(plength-1)),
            meshrange=(msc):mstep:(msc+mstep*(mlength-1)),
            asds = asds, mrange=ms, srange=ss, layers=layers, shifts=shifts, drng=drng,
            paths     = paths,
            mesh      = mesh,
            priors    = priors,
            Σinfo     = Σinfo,
            dos       = dos,
            jdos      = jdos,
            ldos      = ldos,
            shdos     = shdos,
            shg       = shg,
            cpge      = cpge,
            optσ      = optσ,
            sections  = section_types,
            timestamp = timestamp,
            batchdir   = batchdir,
            exportpath  = batchdir*"/$name/out",
            job         = name,
            batch=true, reporting=false, serving=false, presenting=false,
            saving=saving, plotting=plotting,
        )

        if series==:first
            margs = (comargs = [(1:1, 2:mmax)],)
        elseif series==:principal
            margs = (comargs = [(1:mmax, 1:1)],)
        elseif series==:full
            margs = (comargs = [(0:0, 0:0),(mmax:-1:1, 1:1),(1:1, 2:mmax),(0:0, 1:1)],)
        elseif series==:hull
            margs = (comargs = Lattice.hull_comargs(drng...),)
        else
            margs = (arg=1,)
        end

        return programargs(;merge(margs,pargs)...)
end

function build(args,idx=0)
    #Build the Program
    program = Program.build(args,idx)
    Program.save(program)
    return program
end

function sample!(program, sampling::Symbol,idx::Int64)
    #Sampling
    Analysis.do_sampling(sampling::Symbol,program,idx)
    Program.save(program)
    nothing
end

function integrate!(program, datatype::(Type{T} where T <: Analysis.QuadratureData), idx::Int64)
    #Integration of Model Specific Quantities (projected densities of states, topological indices)
    Analysis.do_quad(datatype, program,idx)
    Program.save(program)
    nothing
end

function thermo!(program,idx)
    #Equilibrium Thermodynamic Analysis,
    Analysis.get_μs(program, idx)
    Program.save(program)
    nothing
end

#Flat Out Program Executation
function batch(args)
    #Build the Program
    program = Program.build_program(args)

    for (i,ex)∈ enumerate(program["executables"])
        println("\nBuilding ASD $(ex["utilities"]["id"])"); flush(stdout)
        Program.build_asd(program,i)
        println("Building Model $(ex["utilities"]["id"])"); flush(stdout)
        Program.build_model(program,i)
        #Path Sampling
        Analysis.do_sampling(:path,program,i)
        flush(stdout)
        #Mesh Sampling
        #Analysis.do_sampling(:mesh,program)
        #Integration of Model Specific Quantities (projected densities of states, topological indices)
        Analysis.do_quad(Analysis.DosData{Int64,Float64}, program, i)
        flush(stdout)
        #Equilibrium Thermodynamic Analysis,
        Analysis.get_μs(program,i)
        #Do Non-Equilibrium Optical Response
        for section ∈ args.ex.sections
            Analysis.do_quad(section,  program, i)
            flush(stdout)
        end
        save(program)
        parse(Bool, open(readlines, "$(program["utilities"]["exportpath"])/state.txt")[1]) ? nothing : break
    end
    println("\n")
    flush(stdout)
    #Model Overview Analysis
    #Analysis.post_process_models(program)
    #Run Program Overview Analysis
    #Analysis.post_process_program(program)
    #Run Metric Analysis routines and do the exports
    #Analysis.basic_metric_analysis(program)
    #Export the program structure to a BSON file
    #Return the computed program
    program
end

#Wrapper API Handle for args tuple

function do_job(args)
    println("starting job")
    println("name\t calc\t plot\t export\t ")
    flush(stdout)
    #Execute the calculation
    mkpath("$(args.ex.batchdir)/$(args.ex.job)/out")
    stats  = @timed batch(args)
    #Timing Information
    println("time: "*string(round(stats.time;digits=3))*" s, or "*string(round(stats.time/60;digits=3))*" m, or "*string(round(stats.time/3600;digits=3))*" h")
    println("or "*string(round(stats.time*nprocs()/3600;digits=3))*" core hours")
    flush(stdout)
    #Return Evaluated Program
    return stats.value
end

#=
#Wrapper API Handle for args tuple
function do_job(arg0=(arg="default",))
    args=merge((batchdir=pwd(),job="job",sections=Symbol[],), arg0)
    println("starting job")
    println("name\t calc\t plot\t export\t ")
    #Execute the calculation
    mkpath("$(args.batchdir)/$(args.job)/out")
    stats  = @timed Program.batch(; args...)
    #Timing Information
    println("time: "*string(round(stats.time;digits=3))*" s, or "*string(round(stats.time/60;digits=3))*" m, or "*string(round(stats.time/3600;digits=3))*" h")
    println("or "*string(round(stats.time*nprocs()/3600;digits=3))*" core hours")
    #Return Evaluated Program
    return stats.value
end
=#
#=
const a=.2512
const c=.323
shiftv(n1,n2) = n1*[a, 0, 0 ]+n2*[ a/2, sqrt(3)/2*a, 0 ]
shiftset(N) = (N==0 ? [[[0.0,0.0,0.0],[0.0,0.0,0.0]]] : [[[0.0,0.0,0.0], shiftv(i/(sqrt(3)*N),i/(sqrt(3)*N))] for i=0:N][:])
=#

#=

```bash
#!/bin/bash

echo \"************************************************************\" >> $BATCH_DIR/$1/bin/run-$JID.o
echo \"run on \$(date) in \$PWD\"  >> $BATCH_DIR/$1/bin/run-$JID.o
echo '-host $MACHINE -j $1 -p $2 -jobargs $3 -qargs "${@:4}"'  >> $BATCH_DIR/$1/bin/run-$JID.o

export GKSwstype=\"100\"
export GRDIR=\"\"

println(\"Load_PATH: \"*string(LOAD_PATH))
println(\"Compiling Modules \")
flush(stdout)
@everywhere using GSK, TensorCharts, Lattice, Model, Analysis, Program
using FlakeStacks
flush(stdout)
do_job(jobargs(\"warmup\";batchdir=\"$BATCH_DIR\"))
flush(stdout)
do_job(jobargs(\"$1\"; batchdir=\"$BATCH_DIR\", $3))

```


Goal here is to establish a job stepper and
then use this running loop to continuously step through jobs

a part of this that is relevant is that I think I should strip out program as it is and start using the
=#

function hull_stepper(args0)
    open("bool-run.txt","w") do io
        write(io, "true")
    end



    flag = parse(Bool, open(readlines, "bool-run.txt")[1])
    i=1
    args = args0
    while flag

        print("hi")
        flag = parse(Bool, open(readlines, "bool-run.txt")[1])
        timedwait(()->false,1)
        args = next_args(args)
        p = jobargs()|>do_job

        i+=1
    end

end
