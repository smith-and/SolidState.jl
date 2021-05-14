#export gsk, sYlm

#Generalize Slater Koster Hopping Functions
#radial functions
function rbinomial(n::Int,k::Int)
    n<0 ? (-1)^k*gamma(k-n)/(gamma(k+1)*gamma(-n)) : gamma(n+1)/(gamma(k+1)*gamma(n-k+1))
end

function amid(l::Int,m::Int)
    (l,m)==(0,0) ? 1 : sum(2*s+1 for s=0:max((l-1),0))+m+l+1
end

function constreg(a,d,v,e)
    function(r)
        norm(r)>1e-10 ? v*a : e
    end
end

function constregfunctions(energyscales)
    refcouplingfs = map((a,d,v,e)->constreg(a,d,v,e),energyscales.as,energyscales.ds,energyscales.vs,energyscales.es);
    function(l1, l2, m1, m2,r)
        refcouplingfs[amid(l1,m1),amid(l2,m2)](r)
    end
end


function expreg(a,d,v,e)
    function(r)
        norm(r)>1e-10 ? v*exp(-(norm(r)-a)/d) : e
    end
end

function expregfunctions(energyscales)
    refcouplingfs = map((a,d,v,e)->expreg(a,d,v,e),energyscales.as,energyscales.ds,energyscales.vs,energyscales.es);
    function(l1, l2, m1, m2,r)
        refcouplingfs[amid(l1,m1),amid(l2,m2)](r)
    end
end

function covreg(a,d,v,e)
    function(r)
        norm(r)>1e-10 ? v/(1+(norm(r)-a)/d)^3 : e
    end
end

function covregfunctions(energyscales)
    refcouplingfs = map((a,d,v,e)->covreg(a,d,v,e),energyscales.as,energyscales.ds,energyscales.vs,energyscales.es);
    function(l1, l2, m1, m2,r)
        refcouplingfs[amid(l1,m1),amid(l2,m2)](r)
    end
end;

#spin weighted spherical harmonics
function infTest(x)
    abs(x)==Inf ? 0 : x
end

function swhNorm(l,m,s)
    sqrt((2l+1)/(4π)*gamma(l+m+1)*gamma(l-m+1)/gamma(l+s+1)/gamma(l-s+1))
end

function swhPolar(l,m,s,θ)
    (abs(m)>l)||(abs(s)>l) ? 0 :
        begin
            x=0
            for r=0:1:(l-s)
                x+=binomial(l-s,r)*binomial(l+s,r+s-m)*(-1)^r*infTest((sin(θ/2))^(2l-(2r+s-m)))*(cos(θ/2))^(2r+s-m)
            end
            x
        end
end

function sYlm(l,m,s)
    function(θ,ϕ)
        ((abs(m)>l)||(abs(s)>l) ? 0 : (exp(complex(0,m*ϕ))*swhNorm(l,m,s)*(-1)^(l-s+m)*swhPolar(l,m,s,θ)))
    end
end
#Generalize Slater Koster Functions
function gskAngular(m1,m2,θ,ϕ,l,lp,m,mp)
    4π/sqrt((2*l+1)*(2*lp+1))*(-1)^(m+m1)* exp(complex(0,(m1-m2)*ϕ))*sYlm(l,-m,m1)(θ,ϕ)*sYlm(lp,mp,-m2)(θ,ϕ)
end

function gsk(radials)
    function(l,m,lp,mp)
        function(r3d)
            x=0; for  m1=-l:1:l, m2=(-lp):1:lp
                x+=radials(l,lp,m1,m2,norm(r3d))*gskAngular(m1,m2,atan(sqrt(r3d[1]^2+r3d[2]^2),r3d[3]),atan(r3d[2],r3d[1]),l,lp,m,mp)
            end
            x
        end
    end
end

#=

function cJ(J,m1,m2,l,lp,m,mp)
    (abs(mp-m)>J)|| (abs(m2-m1)>J) ? 0 : (-1)^(m+m1)*sqrt(4*π*(2*J+1))*wigner3j(l,lp,J,-m,mp,mp-m)*wigner3j(l,lp,J,m1,-m2,m2-m1)
end

function gskAngular0(m1,m2,θ,ϕ,l,lp,m,mp)
    x=0.0
    for  J=abs(l-lp):1:(l+lp)
        x+=cJ(J,m1,m2,l,lp,m,mp)*sYlm(J,mp-m,m1-m2)(θ,ϕ)
    end
    x*exp(complex(0,(m1-m2)*ϕ))
end


=#
