using SolidState: cθ, TwistedTriangularGeometry, BNAB

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

dimx((m,n)) = (TwistedTriangularGeometry(BNAB()["blv"],(m,n))["blv"]|>det)/(BNAB()["blv"]|>det)
function hull_series(mmin,mmax)
        args = union(vcat(principal_series(mmin,mmax),first_series(mmin,mmax)))
        sizeargs = dimx.(args)

        args[sortperm(sizeargs)]
end

function bulkhead_series(dmin,dmax)
        #need way to relate these to the dcut
        (mm,sm) = (200,200)
        θid = [(m,m+s) for m∈0:mm for s∈0:sm]
        θsp = [cθ(m,m+s)*180/π for m∈0:mm for s∈0:sm]
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


function make_models(asds::Symbol,comargs::Vector{Tuple{Int,Int}}; cachedir="$(@__FILE__)/cache")
        rootdir = isdir("$cachedir/$asds") ? "$cachedir/$asds" : mkpath("$cachedir/$asds");
        asd0    = eval(Expr(:call,asds))
        made_models = model_id.(readdir(rootdir))
        println("");flush(stdout)
        println("Making $asds Models in $cachedir");flush(stdout)
        for mn ∈ comargs
                if mn ∉ made_models
                        asd = SolidState.CommensurateASD(asd0,mn);
                        hd  = TightBindingDensity(asd)
                        bson("$rootdir/asd-$(mn[1])-$(mn[2]).bson",asd)
                        data_export("$rootdir/hd-$(mn[1])-$(mn[2]).bson",hd)
                        println("$mn made");flush(stdout)
                else
                        println("$mn already made");flush(stdout)
                end
        end
        bson("$rootdir/comargs.bson",Dict(:cargs=>comargs))
        comargs
end


function make_models(asds::Symbol,series::Symbol,mmin::Int,mmax::Int; cachedir="$(@__FILE__)/cache")
        comargs = twist_series(series,mmin,mmax)
        make_models(asds,comargs;cachedir=cachedir)
end
