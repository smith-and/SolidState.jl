function ASDMap(; asd, m, s, shifts, layers)
    ls_asd  = Lattice.LayerShiftASD(LayerASD(eval(asd)(),layers),shifts)
    if (s==0) & (m==0)
        return eval(asd)()
    else
        return Lattice.CommensurateASD(ls_asd,(m,m+s))
    end
end

## Executables

#Creates a Model Executable
function build_model_executable(idx, model_input, analysis_scope, model_const, program)
   Dict{String,Any}(
      "scope" => Dict{String,Any}(
          "path"                                 => [(nbz=nbz,) for nbz=analysis_scope.psN],#building quad input scope
          "mesh"                                 => [(nbz=nbz,) for nbz=analysis_scope.msN],#building quad input scope
          "Analysis.DosData{Int64,Float64}"      => [(nbz=nbz,) for nbz=analysis_scope.qN], #building mesh sampling input scope
          "Analysis.MoreDosData{Int64,Float64}"  => [(nbz=nbz,) for nbz=analysis_scope.qN], #building mesh sampling input scope
          "Analysis.OpticalData{Int64,Float64}"  => [(nbz=nbz,) for nbz=analysis_scope.qN], #building mesh sampling input scope
      ),
      "input" => Dict{String,Any}(
          "asd"                                  => model_input,
          "model"                                => (Model.TightBindingInfo),
          "analysis"                             => analysis_scope.qN,
          "priors"                               => model_const[:priors],
          "paths"                                => model_const[:spaths],
          "mesh"                                 => model_const[:smeshes],
          "Analysis.DosData{Int64,Float64}"      => model_const[:spectral],
          "Analysis.MoreDosData{Int64,Float64}"  => model_const[:moredos],
          "Analysis.OpticalData{Int64,Float64}"  => model_const[:optical]
      ),
      "utilities" => Dict{String,Any}(
          "id"      => "$(model_input.asd)-m-$(model_input.m)-s-$(model_input.s)-l-$(model_input.layers)",
          "idx"     => idx,
          "plot?"   => model_const[:plotting],
          "save?"   => model_const[:saving],
          "workers" => program["utilities"]["workers"],
          "program" => program
      )
   )
end

#Program Executable
function build_program_executable(args)
    m_span  = args.m
    ex_info = args.ex
    Dict{String,Any}(
        "executables" => fill(Dict{String,Any}(), size(m_span)), #execuatable cache
        "state"             => Dict( #state variables
            "base_built?"   => fill(false, size(m_span)),
            "asd_built?"    => fill(false, size(m_span)),
            "model_built?"  => fill(false, size(m_span))
        ),
        "metrics" => Dict{String,Any}(
            "build" => Dict( #metric caches
                "inputs" => m_span, #enumerated inputs
                "metric" => Array{ NamedTuple{(:time, :bytes, :gctime, :gcstats),Tuple{Float64,Int64,Float64,Base.GC_Diff}},length(size(m_span))}(undef,size(m_span)),
                "plot"   => Dict{Symbol,String}()
            ),
            "model" => Dict( #metric caches
                "inputs" => m_span, #enumerated inputs
                "metric" => Array{ NamedTuple{(:time, :bytes, :gctime, :gcstats),Tuple{Float64,Int64,Float64,Base.GC_Diff}},length(size(m_span))}(undef,size(m_span)),
                "plot"   => Dict{Symbol,String}()
            ),
        ),
        "utilities" => Dict(
            "is_init?"    => true,
            "inputs"      => m_span, #enumerated inputs
            "workers"     => CachingPool(workers()),
            "job"         => ex_info.job,
            "timestamp"   => ex_info.timestamp,
            "exportpath"  => homedir(),
            "jobdir"      => ex_info[:expath],
            "exportpath"  => ex_info[:expath]*ex_info[:timestamp],
        )
    )

end

## Executable Management

function timestampdir(dir)
    mkpath(dir*Dates.format(now(),"/T-yyyy-mm-dd-HH-MM-SS"))
    return dir*Dates.format(now(),"/T-yyyy-mm-dd-HH-MM-SS")
end

function key_merge(p1,p2,key )
    Dict(map((x,y,z)->z=>(typeof(x) <: Dict ? merge(x,y) : vcat(x,y) ), values(p1[key]), values(p2[key]), keys(p1[key]))...)
end

export save, load, merge_program
include("utilities.jl")

#Functions on the Program executable
#Build a model executable for each element of the Model Span
function build_program_base(p_args, program::Dict{String,Any}, idx::Int64)
    program["state"]["base_built?"][idx] || begin
    program["executables"][idx] = build_model_executable(idx,p_args.m[idx], deepcopy(p_args.a), deepcopy(p_args.ex), program)
        program["state"]["base_built?"][idx] = true
    end
end

#ASD forming at an executable level
function build_asd(ex::Dict{String,Any})
   stats = @timed begin push!(ex,"asd"=> ASDMap(; ex["input"]["asd"]...)); nothing end
   #Assign Basis Dimension
   push!(ex,"Hdim"=>(ex["asd"]|>ASDBasics|>x->x["sqo"]|>length))
   for i=1:length(ex["input"]["priors"])
       if ex["input"]["priors"][i][1]==:ν
           ex["input"]["priors"][i]=(:ν,ex["asd"]["filling"],ex["asd"]["filling"],1)
       end
   end
   #Assign Proper indicies for the ldos TensorChart (i.e. {(i) : i ∈ [i,L] ⊂ ℕ} the layers)
   ldos_indicies = [(i,) for i=1:max(getindex.(ex["asd"]["sites"],Ref(ex["asd"]["sk"]["Layer"]))...)]
   for i=1:length(ex["input"]["Analysis.DosData{Int64,Float64}"])
       if (ex["input"]["Analysis.DosData{Int64,Float64}"][i][1]=="ldos")
           empty!(ex["input"]["Analysis.DosData{Int64,Float64}"][i][4])
           for j=1:length(ldos_indicies)
               push!(ex["input"]["Analysis.DosData{Int64,Float64}"][i][4],ldos_indicies[j] )
           end
       end
   end

   ex["export/path"] = "$(ex["utilities"]["program"]["utilities"]["exportpath"])/$(ex["utilities"]["id"])"
   mkpath(ex["export/path"])


   Base.structdiff(stats,(value=1,))
end

function build_asd(program::Dict{String,Any},idx::Int64)
    if !program["state"]["asd_built?"][idx]
        #Build the Atomic Site Description
        program["metrics"]["build"]["metric"][idx] = build_asd(program["executables"][idx])
        program["state"]["asd_built?"][idx] = true
    end
end

#Model Making at an executable level
function build_model(ex::Dict{String,Any})
   stats = @timed push!(ex,"model"=>ex["input"]["model"](ex["asd"]));
   Base.structdiff(stats,(value=1,))
end


#build or update base executables asds
function build_model(program::Dict{String,Any}, idx::Int64)::Nothing
    if !program["state"]["model_built?"][idx]
        #Build the Model
        program["metrics"]["model"]["metric"][idx] = build_model(program["executables"][idx])
        program["state"]["model_built?"][idx] = true
    end
    nothing
end

function build_program(args)::Dict{String, Any}
    #Builds the program executable
    program = build_program_executable(args);
    #Builds the model executables and their Analysis executables
    for idx ∈ eachindex(program["utilities"]["inputs"])
        build_program_base(args, program, idx)
    end

    mkpath(program["utilities"]["exportpath"])
    open("$(program["utilities"]["exportpath"])/state.txt","w") do io
        write(io, "true")
    end

    return program
end

function build(args, idx::Int64=0)::Dict{String, Any}
    #Builds the program executable
    program = build_program(args);

    #Builds the model executables and their Analysis executables
    ex_range = (idx==0 ? eachindex(program["utilities"]["inputs"]) : idx:idx)
    for idx ∈ ex_range
        build_asd(program, idx);
        build_model(program, idx);
    end
    return program
end


function save(program::Dict{String,Any})
    bson(program["utilities"]["exportpath"]*"/program.bson",program = program)
end

function load(path::String,)
    BSON.load(path*"/program.bson")[:program]
end

function load(job::String,timestamp::String, batchdir::String=pwd())
    BSON.load(batchdir*"/$job/out/"*timestamp*"/program.bson")[:program]
end

function load(job::String, id::Int64, batchdir::String=pwd())
    timestamp = readdir(batchdir*"/$(job)/out")[end-id]
    return BSON.load(batchdir*"/$job/out/"*timestamp*"/program.bson")[:program]
end

function merge_program(p1,p2)
    timestamp = Dates.format(now(),"T-yyyy-mm-dd-HH-MM-SS")
    Dict{String,Any}(
        "is_init?"    => true,
        "workers"     => CachingPool(),
        "job"         => p1["job"],
        "jobdir"      => p1["jobdir"],
        "commpath"    => p1["commpath"],
        "timestamp"   => timestamp,
        "exportpath"  => p1["jobdir"]*"/"*timestamp,
        "executables" => vcat(p1["executables"],p2["executables"]),
        "inputs"      => vcat(p1["inputs"],p2["inputs"]),
        "state"       => key_merge(p1,p2,"state"),
        "model:metrics" => key_merge(p1,p2,"model:metrics"),
        "build:metrics" => key_merge(p1,p2,"build:metrics"),
        "plot"          => merge(p1["plot"],p2["plot"])
    )
end
