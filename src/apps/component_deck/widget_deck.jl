
function component_deck()
    comps = Dict{Symbol,Any}(
            :update           => button("update"),
            :add_tab          => button("add tab"),
            :clear_tab        => button("clear tab"),
            :add_row          => button("add row"),
            :clear_component  => button("clear"),
            :export           => button("export"),
            :widgetf          => dropdown([abs_angle_selector,state_projector_widget,extraction_widget]),
            :components       => OrderedDict{Int64,Any}()
    )

    on(comps[:add_tab]) do i
        push!(comps[:components], i => vbox())
        comps[:tabs][] = length(comps[:components]|>keys)
        comps[:update][]+=1
    end

    # on(comps[:inputs]) do i
    #     comps[:update][]+=1
    # end

    on(comps[:clear_component]) do i
        empty!(comps[:components])
        comps[:add_tab][]=1
        comps[:update][]+=1
    end

    comps[:tabs]=tabulator(comps[:components])

    on(comps[:clear_tab]) do i
        tabN = comps[:tabs][]
        comps[:components][tabN] = vbox()
        comps[:update][]+=1
    end

    on(comps[:add_row]) do i
        tabN = comps[:tabs][]
        if comps[:inputs][][:widgetf] == extraction_widget
            comps[:components][tabN] = vbox(comps[:components][tabN],comps[:inputs][][:widgetf](comps[:inputs][]))
        else
            comps[:components][tabN] = vbox(comps[:components][tabN],comps[:widgetf][]())
        end
        comps[:update][]+=1
    end

    on(comps[:export]) do i

    end

    output = map(comps[:update]) do i
        tabN = comps[:tabs][]
        comps[:tabs]=tabulator(comps[:components])
        comps[:tabs][] = tabN

        comps[:tabs]
    end

    w = Widget(comps, output=output)

    @layout! w vbox(hbox(:add_tab,:clear_tab,:add_row,:clear_component,:widgetf),observe(_))

    w
end
