module Report

    using Plots
    using ..SolidState


    function pull(name)
        target = "$(ENV["HOME"])/Dropbox/Graduate/reports/$(name)"
        link = "$(ENV["HOME"])/Dropbox/Graduate/completion/reports/$(name)"

        isfile(link) ? Base.symlink(target,link) : println("$name already pulled")
    end

    function share(name)
        target = "$(ENV["HOME"])/Dropbox/Graduate/completion/reports/$(name)"
        dest = "$(ENV["HOME"])/Dropbox/Graduate/reports/$(name)"

        Base.mv(target, dest)
        Base.symlink(dest,target)
    end

    function new(name)
        template = "$(ENV["HOME"])/Dropbox/Graduate/reports/Template"
        dest = "$(ENV["HOME"])/Dropbox/Graduate/completion/reports/$(name)"
        if isdir(dest)
            @info("report $name already exists")
        else
            @info("creating report at $dest")
            Base.cp(template,dest, force=true)
        end
    end

    function rm(name)
        Base.rm("$(ENV["HOME"])/Dropbox/Graduate/completion/reports/$name")
        Base.rm("$(ENV["HOME"])/Dropbox/Apps/Overleaf/$name", recursive=true)
    end

    function pdf(plt, project, name)
        target = "$(ENV["HOME"])/Dropbox/Graduate/completion/reports/$project/images"|>mkpath
        @info("saved $target/$name")
        Plots.pdf(plt, "$target/$name")
    end

    function book_save(plotbook, plotdir)
        for key ∈ keys(plotbook)
            Plots.pdf(getindex(plotbook,key), "$plotdir/$key")
        end
    end

    function compile(project)
        target = "$(ENV["HOME"])/Dropbox/Graduate/completion/reports/$project"


    end

    # Testing
    # Report.new("TEST")
    # Report.share("TEST")
    # Report.pull("TEST")
    # Report.rm("TEST")

end
