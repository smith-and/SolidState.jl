
## Finding just the chemical potential

function fermi_dirac(ϵ,μ,T)
    if T<1e-10
        return (ϵ-μ) > 0 ? 0.0 : 1.0
    else
        return 1/(1+exp((ϵ-μ)/(8.617*1e-5*T)))
    end
end

function occupational(ωs, dos)
    function occupation(μ,T)
        dos .* fermi_dirac.(ωs, Ref(μ), Ref(T))
    end
end

function chemical_potential(dos)
    #Get the original energy domain
    ωs = getindex.(dos.base[:,1],1)
    r_max = max(ωs...)
    r_min = min(ωs...)
    #make a finer mesh
    base = LinRange(r_min, r_max,10000)
    Δω = base[2] - base[1]

    #get the dos data & spline it
    dosdata = getindex.(-1 .* imag.(dos.data[1,1,1,:,1]),1)
    dos_spline = CubicSpline(ωs,dosdata)

    #Sample the spline and get finite temperature function
    sdos = dos_spline[base]|> x-> x ./ (sum(x)*Δω)

    #Do the numerical Legendre Transform and return the function
    occ = occupational(base,sdos)
    occ_lvlset(T,ν) = μ-> (sum(occ(μ, T))*Δω-ν)
    μ(T,ν) = find_zero(occ_lvlset(T,ν),(r_min,r_max))

    return μ
end

#executable interface
function get_ex_μs(ex)
    #Get Priors
    series_tag = string("-s",length(ex["ex-Analysis.DosData{Int64,Float64}"]["result"]))
    ce_priors = ex["input"]["priors"]
    priors = TensorCharts.prior_base_domain(ce_priors)
    new_priors = fill((0.0,0.0), size(priors))
    #Load dos data and
    try
        dos_chart = TensorCharts.load(DosData{Int64,Float64}, ex["ex-Analysis.DosData{Int64,Float64}"]["result"][end]).dos
        μ_f = chemical_potential(dos_chart)

        for i=eachindex(priors)
            new_priors[i] = (priors[i][1], μ_f(priors[i]...))
        end
    catch
        for i=eachindex(priors)
            new_priors[i] = (priors[i][1], 0.0)
        end
    end
    push!(ex["input"], "ce_priors"=>ce_priors)
    push!(ex["input"], "priors"=>new_priors)

    nothing
end

function get_μs(program,i)
    get_ex_μs(program["executables"][i])
end


function get_μs(program)
    for ex ∈ 1:length(program["executables"])
        get_ex_μs(ex,i)
    end
end
