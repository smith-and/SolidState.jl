
function snipe_input(executables)

    w = Dict{Symbol,Any}(
        :update => button("update"),
        :command   => dropdown(["run","walk","slurp"],value="slurp"),
        :comp_type => dropdown(["jit","sys"]),
        :scripts   => dropdown(executables),
        :workers => dropdown(string.([:big,collect(1:100)...])),
        :script_call => textbox(value = "\'main(some code goes here...)\'"),
        :slurm_args  => textbox(value = ""),
        :run         => button("run")
    )

    on(w[:workers]) do i
        w[:update][]+=1
    end

    on(w[:comp_type]) do i
        w[:update][]+=1
    end

    on(w[:script_call]) do i
        w[:update][]+=1
    end

    on(w[:slurm_args]) do i
        w[:update][]+=1
    end

    on(w[:command]) do i
        w[:update][]+=1
    end

    on(w[:scripts]) do dir
        w[:update][]+=1
    end

    on(w[:run]) do i
        run(`sh -c $(w[])`)
    end

    output = map(w[:update]) do i
        string((getindex.(getindex.(Ref(w),[:command,:comp_type,:scripts,:workers,:script_call,:slurm_args])).*" ")...)
    end

    w = Widget(w,output=output)

    @layout! w vbox(hbox(:run,:command,:comp_type,:scripts,:workers), :script_call, :slurm_args, observe(_))
end


function snipe_machine()

    w = OrderedDict{Symbol,Any}(
        :machines => dropdown(["local","dirac","bridges"])
    )

    output = map(w[:machines]) do dir
         scripts = filter(readdir("/home/smitha/Dropbox/Graduate/packages/.scripts/$(dir)")) do file
            file[1]=='.' && file!=".cache" && file!=".side"
        end|> x -> isempty(x) ? [""] : x

        snipe_input(scripts)
    end

    w = Widget(w,output=output)

    @layout! w vbox(:machines,observe(_))
end
