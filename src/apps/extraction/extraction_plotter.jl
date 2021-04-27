

function ext_args(asd,RN,f,brng;)
   #Extraction Parameters

   #Collection Importing and Post-Process Path
   args = (
      outdir = mkpath("/home/smitha/Dropbox/Graduate/packages/.extract/bridges2/$(asd)/.plot/$(RN)/Ext"),
      local_dir  = "/home/smitha/Dropbox/Graduate/packages/.extract/bridges2",
      mount_dir  = "/ocean/projects/phy190028p/asmithc/scripts",
      asd = asd,
      RN  = RN,
      f   = f,
      brng= brng,
      dtype = "Ext",
   )

end

function chart_main_widget(asd,RN,f,brng; plotargs...)
   #Extraction Parameters
   kargs = (
      local_dir  = "/home/smitha/Dropbox/Graduate/completion/bridges2",
      mount_dir  = "/ocean/projects/phy190028p/asmithc/scripts",
      asd = asd,
      RN  = RN,
      f   = f,
      brng= brng
   )

   #Collection Importing and Post-Process Path
   args = (
      outdir = mkpath("$(kargs.local_dir)/$(kargs.asd)/.plot/$(kargs.RN)/Ext"),
      chnl = import_collection(; kargs...),
      dtype = "Ext",
      asd = asd,
      plotargs ...
   )

   #Run Plot Methods on Data
    pdict = Dict{Symbol,Any}(
        :heatmapf=>heatmapf(      kargs.f; args...),
       :cutheatmapf=>cutheatmapf(      kargs.f; xlims=(1.0,3.1), args...),
       :waterfall=>waterfall(     kargs.f; args...),
       :isowaterfall=>isowaterfall(  kargs.f; args...),
       :overlay=>overlay(       kargs.f; args...),
       :scaled_overlay=>scaled_overlay(kargs.f; args...),
       :angle_dep=>angle_dep(     kargs.f, kargs.brng;  args...),
    )

    comps = Dict(
        :selector => dropdown(collect(keys(pdict))),
        :args => args
        )

    output = map(comps[:selector]) do plt
        pdict[plt]
    end

    w = Widget(comps,output=output)

    @layout! w vbox(:selector, observe(_))

    w
end

function extraction_widget(kargs...)
    comps = Dict(
        :extract => button("extract"),
        :asd     => dropdown([:ASD2,:ASD4,:ASD2B]),
        :dtype   => textarea("Run"; value = "SHG-50k", cols =10, rows=1),
        :f       => dropdown([abs,real,angle,imag]),
        :brng    => dropdown([(0.0,1.0),(1.0,3.1),(0.0,4.0)])
    )

    output = map(comps[:extract]) do i
        args = ext_args(Symbol(comps[:asd][]),comps[:dtype][],comps[:f][],comps[:brng][])
        SolidStateApps.slide_wrapper(
            chart_main_widget(Symbol(comps[:asd][]),comps[:dtype][],comps[:f][],comps[:brng][])
            ; args...)
    end

    w = Widget(comps, output=output)

    @layout! w vbox(hbox(:extract,:asd,:f,:brng,:dtype), observe(_))

end
