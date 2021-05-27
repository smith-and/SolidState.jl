function serial_test(A::AbstractArray,B::AbstractArray)
    z = 0.0
    for i ∈ eachindex(A)
        z += A[i]/B[i]
    end
    z
end

function dist_test(A::AbstractArray,B::AbstractArray)
    @distributed (+) for i ∈ eachindex(A)
        A[i]/B[i]
    end
end


function thread_test(A::AbstractArray,B::AbstractArray)
    z = Threads.Atomic{Float64}(0.0)
    Threads.@threads for i ∈ eachindex(A)
        Threads.atomic_add!(z, A[i]/B[i])
    end
    z
end

function test_loops(dim::Int64)
    A = SharedArray(rand(dim,dim))
    B = SharedArray(rand(dim,dim))

    @assert (serial_test(A,B) - dist_test(A,B)) < 1e-5
    @assert (serial_test(A,B) - thread_test(A,B).value)  < 1e-5

    @btime(serial_test($A,$B))
    @btime(dist_test($A,$B))
    @btime(thread_test($A,$B))

end

function test_loops(dims::AbstractRange)

    rootdir0 = mkpath("$(@__DIR__)")
    rootdir  = mkpath("$(rootdir0)/R$(length(readdir(rootdir0))+1)")

    stimes = Vector{Float64}(undef,length(dims))
    dtimes = Vector{Float64}(undef,length(dims))
    ttimes = Vector{Float64}(undef,length(dims))

    sstd = Vector{Float64}(undef,length(dims))
    dstd = Vector{Float64}(undef,length(dims))
    tstd = Vector{Float64}(undef,length(dims))

    for (i,dim) ∈ enumerate(dims)
        A = SharedArray(rand(dim,dim))
        B = SharedArray(rand(dim,dim))

        @assert (serial_test(A,B) - dist_test(A,B)) < 1e-5
        @assert (serial_test(A,B) - thread_test(A,B).value)  < 1e-5

        s_stat = @benchmark(serial_test($A,$B))
        d_stat = @benchmark(dist_test($A,$B))
        t_stat = @benchmark(thread_test($A,$B))

        stimes[i] = mean(s_stat.times)
        dtimes[i] = mean(d_stat.times)
        ttimes[i] = mean(t_stat.times)

        sstd[i] = sqrt(sum((s_stat.times.-stimes[i]).^2)/length(s_stat.times))/2
        dstd[i] = sqrt(sum((d_stat.times.-dtimes[i]).^2)/length(d_stat.times))/2
        tstd[i] = sqrt(sum((t_stat.times.-ttimes[i]).^2)/length(t_stat.times))/2

        println("---------Testing $dim dimension----------");flush(stdout)
        println("Serial: $(stimes[i])");flush(stdout)
        println("Dist: $(dtimes[i])");flush(stdout)
        println("Thread: $(ttimes[i])");flush(stdout)
    end

    plt = plot(title="loop scaling", margins = 10Plots.mm, legend=:topleft)
    plot!(plt, dims,stimes, ribbon=sstd, label="serial")
    plot!(plt, dims,dtimes, ribbon=dstd, label="dist")
    plot!(plt, dims,ttimes, ribbon=tstd, label="thread")

    Plots.pdf(plt, "$(rootdir)/scaling-sdt-$(dims[1])-$(dims[end])")

    plt2 = plot(title="loop scaling", margins = 10Plots.mm, legend=:topleft)
    plot!(plt2, dims,stimes, ribbon=sstd, label="serial")
    plot!(plt2, dims,dtimes, ribbon=dstd, label="dist")
    Plots.pdf(plt2, "$(rootdir)/scaling-sd-$(dims[1])-$(dims[end])")

end

function test_2loops(dims::AbstractRange)

    rootdir0 = mkpath("$(@__DIR__)")
    rootdir  = mkpath("$(rootdir0)/R$(length(readdir(rootdir0))+1)")

    stimes = Vector{Float64}(undef,length(dims))
    dtimes = Vector{Float64}(undef,length(dims))

    sstd = Vector{Float64}(undef,length(dims))
    dstd = Vector{Float64}(undef,length(dims))


    for (i,dim) ∈ enumerate(dims)
        A = SharedArray(rand(dim,dim))
        B = SharedArray(rand(dim,dim))

        @assert (serial_test(A,B) - dist_test(A,B)) < 1e-5
        @assert (serial_test(A,B) - thread_test(A,B).value)  < 1e-5

        s_stat = @benchmark(serial_test($A,$B))
        d_stat = @benchmark(dist_test($A,$B))

        stimes[i] = mean(s_stat.times)
        dtimes[i] = mean(d_stat.times)

        sstd[i] = sqrt(sum((s_stat.times.-stimes[i]).^2)/length(s_stat.times))/2
        dstd[i] = sqrt(sum((d_stat.times.-dtimes[i]).^2)/length(d_stat.times))/2

        println("---------Testing $dim dimension----------");flush(stdout)
        println("Serial: $(stimes[i])");flush(stdout)
        println("Dist: $(dtimes[i])");flush(stdout)
    end

    plt = plot(title="loop scaling", margins = 10Plots.mm, lengend=:topleft)
    plot!(plt, dims,stimes, ribbon =sstd,  label="serial")
    plot!(plt, dims,dtimes, ribbon =dstd,  label="dist")

    Plots.pdf(plt, "$(rootdir)/loop2_scaling-$(dims[1])-$(dims[end])")

    nothing
end
