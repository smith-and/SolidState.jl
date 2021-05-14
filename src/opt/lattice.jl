#export sqo,  R2D, R3D, vrd, LatticePoints, LatticePointsSym, BoundaryRegionQ, BoundaryEdgeQ
#export BoundaryRegionQ, BoundaryEdgeQ, TwistedTriangularGeometry
#export ASDGeometry, ASDBasics, ASDGSK
#export CommensurateASD, ShiftASD, LayerASD


function polar_angles(r3d)
    atan(sqrt(r3d[1]^2+r3d[2]^2),r3d[3]),atan(r3d[2],r3d[1])
end

#Some simple rotations
"""
    function R2D(α::Real,θ::Real)::Matrix{Float64}
"""
function R2D(α::Real,θ::Real)::Matrix{Float64}
    α*[round(cos(θ),digits=10) round(-sin(θ),digits=10) ; round(sin(θ),digits=10) round(cos(θ),digits=10)]
end

function R3D(α::Real,θ::Real)::Matrix{Float64}
    α*[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
end

#Functions to round floating point vectors
function zeroflip(x::Real)::Real
    (x==-0.0 ? 0.0 : x)
end

function vrd(d::Int64)
    function vr(vec::Vector{T} where T <: Real)::Vector{Float64}
        vec.=zeroflip.(round.(vec;digits=d))
    end
end

#Functions to generate sets of lattice points
"""
    LatticePoints(xtal::Tuple{Array{Float64,2},Vector{Vector{T}}} where T <: Real, bnds::Array{Int64,2})::Vector{Vector{Vector{Float64}}}
"""
function LatticePoints(xtal::Tuple{Array{Float64,2},Vector{Vector{T}}} where T <: Real, bnds::Array{Int64,2})::Vector{Vector{Vector{Float64}}}
    vr=vrd(12)
    BravaisPoints=Vector{Float64}[]
    for n1=bnds[1,1]:bnds[1,2], n2=bnds[2,1]:bnds[2,2], n3=bnds[3,1]:bnds[3,2]
        push!(BravaisPoints, (n1*xtal[1][1,:]+n2*xtal[1][2,:]+n3*xtal[1][3,:]))
    end
    unqBravaisPoints=unique(vr,BravaisPoints);
    lattice_points = Vector{Vector{Vector{Float64}}}(undef, length(xtal[2]))
    slv_points = Vector{Vector{Float64}}(undef,length(unqBravaisPoints))
    for i ∈ eachindex(xtal[2])
        for j ∈ eachindex(unqBravaisPoints)
            slv_points[j] = xtal[2][i]+unqBravaisPoints[j]
        end
        lattice_points[i] = copy(slv_points)
    end
    lattice_points
end

function LatticePointsSym(xtal::Tuple{Array{Float64,2},Vector{Vector{T}}} where T <: Real, bnds::Array{Int64,2}, n::Int64)::Vector{Vector{Vector{Float64}}}
    vr=vrd(12)
    BravaisPoints=Vector{Float64}[]
    for k=0:(n-1), n1=bnds[1,1]:bnds[1,2], n2=bnds[2,1]:bnds[2,2], n3=bnds[3,1]:bnds[3,2]
        push!(BravaisPoints, R3D(1,k*2π/n)*(n1*xtal[1][1,:]+n2*xtal[1][2,:]+n3*xtal[1][3,:]))
    end
    unqBravaisPoints=unique(vr,BravaisPoints);
    lattice_points = Vector{Vector{Vector{Float64}}}(undef, length(xtal[2]))
    slv_points = Vector{Vector{Float64}}(undef,length(unqBravaisPoints))
    for i ∈ eachindex(xtal[2])
        for j ∈ eachindex(unqBravaisPoints)
            slv_points[j] = xtal[2][i]+unqBravaisPoints[j]
        end
        lattice_points[i] = copy(slv_points)
    end
    lattice_points
end;

#Functions to obtain certain stars of points
function c6_star(v)
    C6 = R3D(1,π/3)
    (v,C6*v,C6*C6*v,C6*C6*C6*v,C6*C6*C6*C6*v,C6*C6*C6*C6*C6*v)
end

function c6_star_ledger(v, name)
    C6 = R3D(1,π/3)
    vecs  = (v,C6*v,C6*C6*v,C6*C6*C6*v,C6*C6*C6*C6*v,C6*C6*C6*C6*C6*v)
    names = (x->(name*x)).(("1","'1","2","'2","3","'3"))
    Dict{String,Vector{Float64}}(names[i]=>vecs[i] for i=1:length(vecs))
end

#Functions to test bulk or boundary inclusion
function BoundaryRegionQ(boundary)
    let proj=[1 0 0 ; 0 1 0 ; 0 0 0];
        function bulktest(u::AbstractVector)::Bool
            ifelse(reduce(*,
            [ifelse(dot(cross([0,0,1],boundary[i+1]-boundary[i]),
                (proj*u-boundary[i])-(dot(proj*u-boundary[i],(boundary[i+1]-boundary[i]))/(norm(boundary[i+1]-boundary[i]))^2*(boundary[i+1]-boundary[i]))
            )>0,1,0) for i=1:(length(boundary)-1)])>0,true&(!reduce(|,[norm(proj*u-v)<10^-8 for v=boundary]) ),false)
        end
        return bulktest
    end
    return bulktest
end;

function BoundaryEdgeQ(boundary)
    proj=[1 0 0 ; 0 1 0 ; 0 0 0];
    let
        function boundarytest(u::AbstractVector)::Bool
        ifelse(reduce(+,[(ifelse(
                    (norm((proj*u-boundary[i])-(dot(proj*u-boundary[i],(boundary[i+1]-boundary[i]))/(norm(boundary[i+1]-boundary[i]))^2*(boundary[i+1]-boundary[i])))<10^-8),
                1,0)) for i=1:(length(boundary)-1)])>0,
            (true&&(norm(proj*u-last(boundary))>10^-8))&&(norm(proj*u)<(1.1*norm(last(boundary)))),false
        );
        end
        return boundarytest
    end
    return boundarytest
end;

#METHODS TO CREATE NEW ASDs FROM EXISTING ONES
#shift the first length(shift) elements of the sites list by the vector in shifts at the same index
#shifts can be shorter but not longer
function UCShiftASD(ASD,shifts::Vector{Vector{Float64}})
    for i ∈ 1:length(shifts)
        ASD["sites"][i][3] -= shifts[i]
    end
    ASD
end

#shift each layer by shift[l]
function LayerShiftASD(ASD,shifts::Symbol)
    shifts==:none ? (return ASD) : nothing;
end

function LayerShiftASD(ASD,shifts::Vector{Vector{Float64}})
    for i ∈ 1:length(ASD["sites"])
        ASD["sites"][i][3] .+= shifts[ASD["sites"][i][2]]
    end
    ASD
end

#Layer two ASD on top of each other
function LayerASD(ASD,n)
    n==1 && return ASD;
    new_sites = [copy(ASD["sites"]) for _=1:n]
    for i=1:n
        for j ∈ 1:length(new_sites[i])
            l_s= [new_sites[i][j]...]
            new_sites[i][j]=(l_s[1],i,l_s[3].+[0.0,0.0,i*ASD["blv"][3,3]],l_s[4:end]...)
        end
    end
    ASD["blv"][3,3]*=n;
    ASD["sites"]=vcat(new_sites...)
    ASD
end

#METHODS TO EXTRACT INFORMATION FROM ASDs
#Second Quantized Operator (sqo) type
struct sqo{Int_Type <: Integer, Float_type <: Real }
    atom::String
    layer::Int_Type
    pos::Array{Float_type,1}
    spin::Array{Rational{Int_Type},1}
    orbital::Array{Int_Type,1}
    glpyh::Symbol
end

##ASD Manipulations
function ASDGeometry(ASD::Dict{String, Any})::Dict{String, Any}
    #Real Space
        #Lattice Generators
        Λr = transpose(ASD["blv"])
        pos=getindex.(ASD["sites"], Ref(ASD["sk"]["Pos"]));
        xtal=(Λr,pos)
        #Wigner Sietz Unit Cell Corners
        C12 = R3D(1/sqrt(3.0),-π/6)
        uc_corners        = c6_star_ledger(C12*Λr[:,1],"K")
        uc_side_midpoints = c6_star_ledger( 0.5*Λr[:,1],"M")
        uc_high_symmetry = merge(merge(uc_corners, uc_side_midpoints), Dict("Γ"=> [0.0,0.0,0.0]))
        #WS unit cell corner scale
        Lmag=sqrt(dot(xtal[1][1,:],xtal[1][1,:])/3);

        l_idx  = getindex.(ASD["sites"], Ref(ASD["sk"]["Layer"]));
        glyphs = getindex.(ASD["sites"], Ref(ASD["sk"]["Glyph"]));
        colors = getindex.(ASD["sites"], Ref(ASD["sk"]["Color"]));
        atom = getindex.(ASD["sites"], Ref(ASD["sk"]["Atom"]));
        layer = getindex.(ASD["sites"], Ref(ASD["sk"]["Layer"]));
        labels = string.(atom,layer)
    #Momentum Space
        #Lattice Generators
        Λk                  = 2*π*inv(transpose(Λr))
        dxtal               = (Λk,([0.0, 0.0, 0.0]));
        #Brillouin Zone Unit Cell Corners
        C12                 = R3D(1/sqrt(3.0),-π/6)
        bz_corners          = c6_star_ledger(C12*Λk[:,1],"K")
        bz_Q_pts            = c6_star_ledger(C12*Λk[:,1]./2.0,"Q")
        bz_side_midpoints   = c6_star_ledger( 0.5.*Λk[:,1],"M")
        bz_high_symmetry    = merge(merge(merge(bz_corners, bz_side_midpoints), bz_Q_pts), Dict("Γ"=> [0.0,0.0,0.0]))
        #BZ Corner Scale
        Kmag=sqrt(dot(dxtal[1][1,:],dxtal[1][1,:])/3);
    #Output
    Dict{String,Any}(
        "xtal"  =>xtal,
        "Lmag"  =>Lmag,
        "uc_c"  =>c6_star(C12*Λr[:,1]),
        "uc_hs" =>uc_high_symmetry,
        "dxtal" =>dxtal,
        "Kmag"  =>Kmag,
        "bz_c"  =>c6_star(C12*Λk[:,1]),
        "bz_hs" =>bz_high_symmetry,
        "layer" =>l_idx,
        "glyph" =>glyphs,
        "colors"=>colors,
        "labels"=>labels
        )
end

#Local Hopping parameterization from ASD
function ASDGSK(ASD::Dict{String,Any})::Dict{NTuple{N,String} where N,Any}
    if ASD["regulator"]=="cov"
        return Dict(x=>gsk(covregfunctions(ASD["scales"][x])) for x=keys(ASD["scales"]))
    elseif  ASD["regulator"]=="exp"
        return Dict(x=>gsk(expregfunctions(ASD["scales"][x])) for x=keys(ASD["scales"]))
    elseif ASD["regulator"] == "const"
        return Dict(x=>gsk(constregfunctions(ASD["scales"][x])) for x=keys(ASD["scales"]))
    end
end

#Provides the basic extractions from the atomic site description
function ASDBasics(ASD::Dict{String, Any})::Dict{String, Any}
    sname    = getindex.(ASD["sites"], Ref(ASD["sk"]["Atom"]))
    lidx     = getindex.(ASD["sites"], Ref(ASD["sk"]["Layer"]))
    pos      = getindex.(ASD["sites"], Ref(ASD["sk"]["Pos"]))
    spins    = getindex.(ASD["sites"], Ref(ASD["sk"]["Spin"]))
    orbitals = getindex.(ASD["sites"], Ref(ASD["sk"]["Orbital"]))
    glyphs   = getindex.(ASD["sites"], Ref(ASD["sk"]["Glyph"]))
    #Need to
    xtal=(ASD["blv"],pos);
    spinlabels=[[[sd[1],proj] for proj=range(-sd[1],sd[1];step=1)][sd[2]] for sd=spins ];
    orbitlabels=[[[sd[1],proj] for proj=range(-sd[1],sd[1];step=1)][sd[2]] for sd=orbitals ];
    stateInfo=[[[sname[i],lidx[i],pos[i],spin,orbit,glyphs[i]] for spin=spinlabels[i], orbit=orbitlabels[i]] for i=1:length(spinlabels)]
    sqoBasis::Array{sqo{Int64,Float64},1}=cat([cat([sqo(sqod...) for sqod=siteSQD]...,dims=2) for siteSQD=stateInfo]...,dims=2)[1,:];
    #need to convert the sites info to the single particle basis as list of sqo opjects
    Dict{String, Any}(
    "sqo"=>sqoBasis,
    "sld"=>stateInfo,
    "Layers"=>lidx,
    "xtal"=>xtal,
    "snames"=>sname,
    "regs"=>"reg",
    "cutfunction"=>(x::Array{Float64,1}-> LinearAlgebra.norm(x)<=ASD["cutoff"]),
    "gsk"=>ASDGSK(ASD),
    "lbase"=>LatticePointsSym((xtal[1],[[0, 0, 0 ]]),ASD["lbase"]*[0 1 0 ;-1 0 0 ; 0 0 0 ],3)[1]
    )
end;

#Twist Commensurate Structures
function cθ(m,n)
    acos((m^2+n^2+4m*n)/(2*(m^2+n^2+m*n)))
end

function TwistedTriangularGeometry(blv::Array{Float64,2},mn::Tuple{Int64,Int64})
    m=mn[1]; n=mn[2];
    angle=cθ(m,n);
    g1=R3D(1,-angle/2);
    g2=R3D(1,angle/2);
    n1= [-n -(m+n) 0 ; m+n m 0 ; 0 0 1];
    n2= [-m -(m+n) 0 ; m+n n 0 ; 0 0 1];
    ac  = similar(blv)
    bc  = similar(blv)
    ac .= transpose(g2*transpose(blv)*n2);
    bc .= 2π*inv(transpose(ac));
    Dict("gs"=>(g1,g2), "θ"=>angle,"blv"=>ac,"dblv"=>bc,"mn"=>(m,n),"n1"=>n1,"n2"=>n2)
end

function CommensurateUnitCell(ASD,(m,n))
    ASDB=ASDBasics(ASD);
    bnds=max(m,n)*abs(m-n)*[0 1 ; -1 0; 0 0];
    ttg=TwistedTriangularGeometry(ASD["blv"],(m,n));
    #form boundary test functions
    vboundary=[R3D(1/sqrt(3),π/6+k*π/3)*ttg["blv"][1,:] for k=0:6];
    vedge=[R3D(1/sqrt(3),π/6+k*π/3)*ttg["blv"][1,:] for k=-2:0];
    vBndQ=BoundaryRegionQ(vboundary);
    vEdgeQ= BoundaryEdgeQ(vedge);
    voronoiQ = x-> vBndQ(x)||vEdgeQ(x);
    #form lattice cover
    bravaisCover=LatticePointsSym(ASDB["xtal"],bnds,3);
    LayerRotations=map(x->ttg["gs"][x],ASDB["Layers"]);
    twistedCover=map((x,y)->map(u->y*u,x),bravaisCover,LayerRotations);
    vuc=map(x->filter(voronoiQ,x),twistedCover)
end

#Generates a commensurate ASD from Bravais one
function CommensurateASD(ASD,(m,n))
    basedim=length(ASDBasics(ASD)["sqo"])
    if m==0
        newsiteinfo = deepcopy(ASD["sites"])
        for (j,site) ∈ enumerate(ASD["sites"])
            newsite = [deepcopy(site)...]
            newsite[ASD["sk"]["Pos"]] = R3D(1,(site[ASD["sk"]["Layer"]]-1)*π/3)*site[ASD["sk"]["Pos"]]
            newsiteinfo[j] = (newsite...,)
        end

        return Dict("cutoff"=>ASD["cutoff"], "lbase"=>ASD["lbase"], "blv"=>ASD["blv"],
        "sk"=>ASD["sk"], "sites"=>newsiteinfo,
        "regulator"=>ASD["regulator"],
        "scales"=>ASD["scales"], "bd"=>basedim,
        "filling"=>ASD["filling"]
        )
    elseif abs(m-n)==0
        return ASD
    else
        ttg=TwistedTriangularGeometry(ASD["blv"],(m,n));
        #form new site information
        vuc=CommensurateUnitCell(ASD,(m,n));
        sitefunctions=[newposition->([ifelse(i==ASD["sk"]["Pos"],newposition,ASD["sites"][j][i]) for i=1:length(ASD["sites"][j])]...,) for j=1:length(ASD["sites"])];
        newsiteinfo= cat([map(x->sitefunctions[i](x),vuc[i]) for i=1:length(vuc)]...,dims=1);
        #Assemble the ASD info for the new sites
        return Dict("cutoff"=>ASD["cutoff"], "lbase"=>((m<3&(-1 <= m-n <=1)) ? 4 : 2), "blv"=>ttg["blv"],
        "sk"=>ASD["sk"], "sites"=>newsiteinfo,
        "regulator"=>ASD["regulator"],
        "scales"=>ASD["scales"], "bd"=>basedim,
        "filling"=>ASD["filling"]
        )
    end
end

θindex(mmax,smax) = [(m,s) for m∈0:mmax for s∈0:smax]
θcomargs(mmax,smax) = [(m:m,s:s) for m∈0:mmax for s∈0:smax]
θseries(mmax,smax) = [s for m∈0:mmax for s∈0:smax]
θspectra(mmax,smax) = [Lattice.cθ(m,m+s)*180/π for m∈0:mmax for s∈0:smax]
LMspectra(mmax,smax) = [(Lattice.TwistedTriangularGeometry(BNAB()["blv"],(m,m+s))["blv"]|>det)/(BNAB()["blv"]|>det) for m∈0:mmax for s∈0:smax]

function hull_comargs(dmin,dmax)
        #need way to relate these to the dcut
        (mm,sm) = (200,200)
        θid = θindex(mm,sm)
        θca = θcomargs(mm,sm)
        θsp = θspectra(mm,sm)
        θLM = LMspectra(mm,sm)
        LMmask= 0 .< θLM .< dmax
        unqθ = (union(round.(θsp[LMmask],digits=4)))
        θspargs = [findall(x->round(x,digits=4)==θ,θsp[LMmask]) for θ∈unqθ]
        θspcopies=getindex.(Ref(θsp[LMmask]),θspargs)
        θLMcopies=getindex.(Ref(θLM[LMmask]),θspargs)
        θcacopies=getindex.(Ref(θca[LMmask]),θspargs)

        θhullca = Vector{typeof(θca[1])}(undef,length(unqθ))
        θhullsp = Vector{typeof(θsp[1])}(undef,length(unqθ))
        θhullLM = Vector{typeof(θLM[1])}(undef,length(unqθ))
        for (i,LMs) ∈ enumerate(θLMcopies)
                θhullca[i] = θcacopies[i][argmin(LMs)]
                θhullLM[i] = θLMcopies[i][argmin(LMs)]
        end

        return θhullca[dmin .< θhullLM .< dmax]
end
