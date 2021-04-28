### Run the precompile script
function test_compile_run()
    include("$(@__DIR__)/precompile_script.jl")
end

"""
    base_compile(cachedir::String=pwd())

Add some basics ([:BSON, :OrderedCollections, :LinearAlgebra, :Plots, :ColorSchemes]) to the current default system image
"""
function base_compile(cachedir::String=pwd())
    println("compiling a system image")
    cdir = mkpath("$cachedir/.cache/system");
    stats = @timed PackageCompiler.create_sysimage([:BSON, :OrderedCollections, :LinearAlgebra, :Plots],
        sysimage_path               = "$cdir/sysimage0.dylib",
        precompile_execution_file   = "$(@__DIR__)/precompile_script.jl",
        #project                     = "$(ENV["HOME"])/.julia/environments/v1.5/",
        incremental                 = true
    )
    println("build time \t"*string(stats.time/60)*"m")
    nothing
end


"""
    prestack(cachedir::String=pwd())

Make sysimage for modules [:Plots,:SolidState,:SolidStateApps] and in cachedir
"""
function prestack(cachedir::String=pwd())
    println("compiling a system image")
    cdir = mkpath("$cachedir/.cache/system");
    stats = @timed PackageCompiler.create_sysimage([:BSON, :OrderedCollections, :LinearAlgebra, :Plots, :SolidState],
        sysimage_path               = "$cdir/sysimage.dylib",
        precompile_execution_file   = "$(@__DIR__)/precompile_script.jl",
        #project                     = "$(ENV["HOME"])/.julia/environments/v1.5/",
        incremental                 = true
    )
    println("build time \t"*string(stats.time/60)*"m")
    nothing
end
