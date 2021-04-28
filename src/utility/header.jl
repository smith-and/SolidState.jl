using SolidState: cθ, DOS, LP, SHG, make_models, DataChart

#Print Args
"""
    dict_print(args)

print the different keys and calues of a NamedTuple...bad name
"""
function dict_print(args)
    for  (i,key) ∈ enumerate(keys(args))
        println("$key: $(args[i])")
    end
    flush(stdout)
end

#Program Args
"""
    program_args(jobname::Symbol; headdir = pwd(), cachedir="\$headdir/.cache", pool::Union{AbstractWorkerPool,Symbol}=default_worker_pool(), key="R", kargs...)

Returns a NamedTuple:
```julia
kargs = (
    jobname          = "\$headdir/\$jobname",
    cachedirr       = mkpath(cachedir),
    scriptdir       = mkpath("\$headdir/\$jobname/out/\$key\$keyN"),
    plotdir         = mkpath("\$headdir/\$jobname/plot/\$key\$keyN"),
    pool            = poolz, #:none,
)
```
"""
function program_args(jobname::Symbol; headdir = pwd(), cachedir="$headdir/.cache", pool::Union{AbstractWorkerPool,Symbol}=default_worker_pool(), key="R", kargs...)
     # asd header dir

    if typeof(pool)==Symbol
        poolz=pool
    elseif length(pool)==0
        addprocs(1)
        poolz = default_worker_pool()
    else
        poolz = pool
    end

    keyN = length(filter(x->occursin(key,x),readdir(mkpath("$headdir/$jobname/out"))))+1

    kargs = (
        jobname          = "$headdir/$jobname",
        cachedirr       = mkpath(cachedir),
        scriptdir       = mkpath("$headdir/$jobname/out/$key$keyN"),
        plotdir         = mkpath("$headdir/$jobname/plot/$key$keyN"),
        pool            = poolz, #:none,
    )
end

#Model Args
"""
    model_args(asd::Symbol,(series,minn,maxn)::Tuple{Symbol,Int,Int},pargs)

    Calls make models for the series information and appends its output and the asd symbol to pargs, intended as some program args
"""
function model_args(asd::Symbol,(series,minn,maxn)::Tuple{Symbol,Int,Int},pargs)
    margs = (
        asd             = asd,
        comargs         = make_models(asd,series,minn,maxn,cachedir=pargs.cachedirr),
    )
    args = merge(pargs, margs)
    dict_print(args)

    args
end

export logstep
"""
    function logstep(p,n)::AbstractRange
"""
function logstep(p,n)::AbstractRange
    p^(n-1):p^(n-1):p^n
end

export log2step
"""
    function log2step(p::Int,n::Int)::AbstractRange
"""
function log2step(p::Int,n::Int)::AbstractRange
    p^(n-1):9*p^(n-1):p^n
end

export log3step
function log3step(p::Int,n::Int)::AbstractRange
    p^(n-1):3*p^(n-1):p^n
end
