# Add Save Button Component
function add_plot_save!(d,args)
    push!(d[:buttons],:savebutton)

    push!(d,:savebutton  => button("Save"))
    push!(d,:save_path   => textbox(; value = "plot"))

    p_path = pwd()

    try
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])")
    catch
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])")
    end

    on(d[:savebutton]) do b
        Plots.pdf(d[:w0][],"$(p_path)/"*d[:save_path].output[]*"-$b")
        Plots.png(d[:w0][],"$(p_path)/"*d[:save_path].output[]*"-$b")
        d[:update][]+=1
    end

    nothing
end

# Add Figure Tex Component
function add_plot_figure!(d,args)
    push!(d,:fig_caption => Widgets.textarea("figure caption";rows=5,cols=30))

    p_path = pwd()
    try
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])/$(args[:asd])-$(args[:mn][1])-$(args[:mn][2])")
    catch
        p_path = mkpath("$(args[:outdir])/$(args[:asd])/plot/$(args[:dtype])")
    end

    push!(d,:fig_tex => Widgets.textarea("caption"))
    on(d[:fig_caption]) do b
        d[:fig_tex].output[] = "\\begin{figure}
            \\includegraphics[width=1\\textwidth]{$(p_path)/$(d[:save_path].output[])-$(d[:savebutton][])}
            \\caption{\\label{test} $(d[:fig_caption].output[])}
        \\end{figure}"
        d[:update][]+=1
    end

    nothing
end

# Add Beamer Frame Tex Component
function add_fig_slide!(d)
    push!(d,:slide_text => Widgets.textarea("slide text";rows=5,cols=30))
    push!(d,:slide_tex  => Widgets.textarea("tex output", rows=10,cols=30))
    on(d[:slide_text]) do b
        d[:slide_tex].output[] = "
\\begin{frame}
    \\begin{columns}
        \\begin{column}{0.5\\textwidth}
            $(d[:slide_text].output[])
        \\end{column}
        \\begin{column}{0.5\\textwidth}
            $(d[:fig_tex].output[])
        \\end{column}
    \\end{columns}
\\end{frame}
        "
        d[:update][]+=1
    end

    nothing
end

#Add Callback and return a Widget

#Slide Creation tools
function slide_wrapper(w0; args...)
    # Creates a set of components which need to be
    components = OrderedDict{Symbol,Any}(
            :buttons => Symbol[],
            :update  => button("update"),
            :w0      => w0
            )

    add_plot_save!(components, args)
    add_plot_figure!(components, args)
    add_fig_slide!(components)

    output = map(components[:update]) do x
        components[:w0]
    end

    w = Widget(components, output=output)


    # Set layout of widget
    @layout! w hbox(
            observe(_),
            vbox(
                hbox(:savebutton,:save_path),
                :fig_caption,
                :slide_text
            ),
            :slide_tex
        )


    w
end
