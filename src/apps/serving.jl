### Widget Composition

"""
    component_map(inputs)

Pass an OrderedDict containing dictionaries which are as widget descriptions and can be called as

    inputs[i][:widget](intputs[i])

The widgets are then displayed vertically in the order they are listed, in one vbox
"""
function component_map(inputs)
    components = Dict{Symbol,Any}()
    for key in keys(inputs)
        push!(components, key => inputs[key][:widgetf](inputs[key]))
    end

    components
end

### Mux Utilities.
### Currently really only at a sinlge page

"""
    appserve(app::Function)

Starts a firefox tab with a single page webio server displaying the output of app
"""
function appserve(app::Function)
    port = rand(8000:9000)
    webio_serve(page("/", app), port)
    run(`firefox localhost:$port`)
end



"""
    widget_serve(widget)

Serve widget in firefox tab on single page
"""
function widget_serve(widget)
    SolidStateApps.appserve() do req
        widget
    end
end
