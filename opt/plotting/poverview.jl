
ndiff1(data,step) = data./step|>diff
ndiff2(data,step) = ((data|>diff)[2:end].+(data|>reverse|>diff|>reverse)[1:end-1])./(step^2)

function crit_analysis(dom,data)
    #Spline Domain
    sΔω   = 0.005
    sdom  = dom[1]:sΔω:dom[end]
    #Spline Function
    println(sdom)
    tc_spline = CubicSpline(dom,data)
    sdata   = tc_spline[sdom]
    tc_f(x) = tc_spline[x]

    #First Derivative
    tc_d1_spline = CubicSpline(sdom[1:end-1], ndiff1(sdata,sΔω))
    tc_d1_f(x) = tc_d1_spline[x]
    #Critical Points
    tc_crit_pts=find_zeros(tc_d1_f,sdom[3],sdom[end-1])
    extremal_vals=tc_f.(tc_crit_pts)
    tc_crits=hcat(tc_crit_pts,extremal_vals)
    #Second Derivative
    tc_d2_spline = CubicSpline(sdom[2:end-1], ndiff2(sdata,sΔω))
    tc_d2_f(x) = tc_d2_spline[x]
    #Inflection Points
    tc_inflection_pts=find_zeros(tc_d2_f,sdom[2],sdom[end-1])
    inflection_pts=tc_f.(tc_inflection_pts)
    inflection_slopes=tc_d1_f.(tc_inflection_pts)
    tc_infl=hcat(tc_inflection_pts,inflection_pts)
    #Results
    Dict{Symbol,Any}(
        :sdom        =>sdom,
        :crits       =>tc_crits,
        :inflections =>tc_infl,
        :infl_slope  =>inflection_slopes,
        :f           =>tc_f,
        :d1          =>tc_d1_f,
        :d2          =>tc_d2_f
    )
end

function waterfall_plot(tcs,labels,idx,linecutidx; title="",xaxis="")
    #make the waterfall plot
    waterfall=plot(framestyle=:box,title=title,xaxis=xaxis,yticks=:none)
    dashmax=Int(floor(length(tcs[1].base[:,1])*.3))
    wf_spacing = max([max(abs.(tcs[i].data[idx...])...) for i=1:length(tcs)]...)
    wf_height  = wf_spacing*((getfield.(tcs,1)|>length))|>x->[0,x]
    for (i,tc) ∈ enumerate(tcs)
        [(labels[i]/100) for _=1:Int(floor(length(tcs[i].base[:,1])*.3))]
        plot!(waterfall,getindex.(tcs[i].base[:,1],1),[(i-1)*wf_spacing for _=1:length(tcs[1].base[:,1])],label="",linestyle=:dot,color=:black,opacity=0.5)
        plot!(waterfall,getindex.(tcs[i].base[:,1],1),abs.(tcs[i].data[idx...]).+(wf_spacing*(i-1)),label="",color=i)
    end
    plot!(waterfall,tcs[1].base[:,1][1][1]|>x->[x,x],wf_height,linestyle=:dashdot,opactiy=0.2,label="",color=:black)
    for idx0 ∈ linecutidx
        plot!(waterfall,tcs[1].base[:,1][idx0][1]|>x->[x,x],wf_height,linestyle=:dashdot,opactiy=0.2,label="",color=:black)
    end
    annotate!([(tcs[1].base[1,1][1],wf_spacing*(i-1),text("$(round(θ,digits=3))ᵒ   ",6,:bottom,:left)) for (i,θ) ∈ enumerate(labels)])
    display(waterfall)
    return waterfall
end

function overlay_plot(tcs,labels,idx,(mindom,maxdom); title="",xaxis="",yaxis="")
    #make the waterfall plot
    overlay=plot(framestyle=:box,xaxis=xaxis, yaxis=yaxis,title=title)
    dashmax=max(1,Int(floor(length(tcs[1].base[:,1])*.3)))
    for (i,tc) ∈ enumerate(tcs)
        plot!(overlay,getindex.(tc.base[max(1,Int(floor(dashmax/3))):max(2,dashmax),1],1),abs.(tc.data[1,1,1,max(1,Int(floor(dashmax/3))):max(2,dashmax),1]),label="",color=i)
        #plot!(overlay,getindex.(tcs[i].base[mindom:maxdom,1],1),abs.(tcs[i].data[idx...])[mindom:maxdom],label="")
    end
    return overlay
end

function angle_dep_plot(tcs,labels,idx,linecutidx; title="",xaxis="",yaxis="")
    data0 = tcs.|>x->(abs(x.data[1,1,1,1,1]))
    angle_dep = plot(labels, data0, frame=:box,title=title, xaxis=xaxis, yaxis=yaxis, linestyle=:auto, color=:black,label="ω=0",opacity=0.5,legend=:outertopright,legendtitle="ω₀(eV)")
    scatter!(labels, data0,label="", color=collect(1:length(labels)))
    for i∈linecutidx
        #round((getfield.(tcs,name)[1].base[i,1][1]),digits=3)
        linecutdata = tcs.|>x->(abs(x.data[1,1,1,i,1]))
        scatter!(angle_dep, labels, linecutdata, label="", color=collect(1:length(labels)))
        plot!(angle_dep, labels, linecutdata,opacity=0.5, linestyle=:auto,color=:black,label="ω=$(round(tcs[1].base[i,1][1],digits=2))")
    end
    return angle_dep
end

function critical_pt_plot(tcs,labels)
    crit_plot=plot()
    for i∈1:length(tcs)
        dom=getindex.(tcs[i].base[:,1],1)
        data=abs.(tcs[i].data[1,1,1,:,1])
        crits = crit_analysis(dom,data)
        plot!(crits[:sdom],crits[:f].(crits[:sdom]),label="",opacity=0.3,color=:gray)
        scatter!(crit_plot,crits[:crits][:,1],crits[:crits][:,2],opacity=1,label="$(round(labels[i],digits=3))ᵒ")
        #scatter!(crits[:inflections][:,1],crits[:inflections][:,2],opacity=0.2,label="inflection")
    end
    return crit_plot
end

function sub_peak_linecuts(tcs,key)
    #domain and data info
    dom=getindex.(getfield.(tcs,key)[1].base[:,1],1)
    data=abs.(getfield.(tcs,key)[1].data[1,1,1,:,1])
    crits = crit_analysis(dom,data)
    idx0 = dom[dom .< crits[:crits][1,1]]|>length
    linecutidx= (Int(-floor((idx0/num))) == 0 ? (1:1:1) : idx0:Int(-floor((idx0/num))):Int(floor(idx0/2)))

    return linecutidx
end

function single_val_linecuts(tcs,key,val)
    #domain and data info
    dom=getindex.(getfield.(tcs,key)[1].base[:,1],1)
    idxplus = length(dom[dom .< val])
    idxminus = length(dom)-length(dom[dom .> val])

    return max(1,idxminus):max(2,(idxminus+1))
end

function overview_angle_dependence(program,section,name, idx=(1,1,1,:,1),num=6)
    tcs = Vector{OpticalData{Int64,Float64}}(undef,length(program["executables"]))
    angles = Vector{Float64}(undef,length(program["executables"]))
    for (i,ex) ∈ enumerate(program["executables"])
        angles[i] = (ex["model_input"].m!=0 && ex["model_input"].s!=0) ? Lattice.cθ(ex["model_input"].m,ex["model_input"].m+ex["model_input"].s)*180/π : 0.0
        tcs[i] = TensorCharts.load(section, ex["ex-"*string(section)]["result"][1])
    end

    #linecutidx = sub_peak_linecuts(tcs,name)
    linecutidx = single_val_linecuts(tcs,name,1.17)

    waterfall = waterfall_plot(getfield.(tcs,name),angles,idx,linecutidx, title="$name abs",xaxis="Driving Frequency ħω (eV)")
    overlay = overlay_plot(getfield.(tcs,name),angles,idx,(linecutidx[1],linecutidx[end]);xaxis="ħω (eV)", yaxis="|χ| V/pm")
    angle_dep = angle_dep_plot(getfield.(tcs,name),angles,idx,linecutidx, xaxis="θᵒ", yaxis="|χ| V/pm")
    #crit_plot = critical_pt_plot(getfield.(tcs,name),angles,idx,linecutidx)

    plot0 = plot(waterfall,angle_dep,overlay,layout=@layout [ a{0.5w} [ b{0.5h} ; c{0.5h}] ])
    display(plot0)

    plotdict = Dict(
        :waterfall => waterfall,
        :overlay   => overlay,
        :linecuts  => angle_dep,
        :panel     => plot0
        )

    mkpath(program["exportpath"]*"/$(name)_angle_dep")
    for plt ∈ keys(plotdict)
        Plots.pdf(plotdict[plt],program["exportpath"]*"/$(name)_angle_dep/"*string(plt))
    end

    return plotdict
end

function post_process_program(program)
    if program["job"]=="sequences"
        overview_angle_dependence(program,OpticalData{Int64,Float64},:shg)
        #overview_angle_dependence(program,OpticalData{Int64,Float64},:shg)
    end
end
