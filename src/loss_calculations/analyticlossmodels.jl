

function graphenetwoplasmonemission(ω::Real, μ::Real; mesh = 100, histogram_width =10, δ=1, conesize=1)
    F2 = 0 
    ucellarea = 5.238760872572826 #Graphene unit cell area in angstrom^2 
    kF = μ/6
    q= 2π*2.5*137/(4*pi*ħ*c*μ)*ω^2 #ev^2/ev^2/angstrom -> inverse angstrom units
    println(q)
    #prefactor = (μ^3/kF^2)*(π^2/ucellarea)*(ω/(6*q))^6
    prefactor = (μ^3/kF^2)*(π^2)*(ω/(6*q))^6
    for kiter in 1:mesh
        k = μ/conesize*kiter/mesh
        for θiter in 1:mesh
            θ = 2π*θiter/mesh
            for bandi in [-1, 1]
                ϵi = 6*bandi*k
                fi = heaviside(μ-ϵi)
                fi ≈ 0 && continue
                for bandf in [1]
                    ϵf = 6*bandf*sqrt((k*cos(θ)+2*q)^2+k^2*sin(θ)^2)
                    ff = heaviside(μ-ϵf)
                    (1-ff) ≈ 0 && continue 

                    ϵm = 6*sqrt((k*cos(θ)+q)^2+k^2*sin(θ)^2)

                    ϕi = θ
                    ϕm = atan((k*sin(θ))/(k*cos(θ)+q))
                    ϕf = atan((k*sin(θ))/(k*cos(θ)+2q))

                    overlapim = 1+bandi*cis(ϕi-ϕm)
                    overlapim *= 1/2
                    overlapmf = 1+bandf*cis(ϕm-ϕf)
                    overlapmf *= 1/2

                    m1 = (overlapim*overlapmf)/(ϵm-ϵi-ω-1im*δ)
        
                    ϵm = -6*sqrt((k*cos(θ)+q)^2+k^2*sin(θ)^2)

                    overlapim = 1-bandi*cis(ϕi-ϕm)
                    overlapim *= 1/2
                    overlapmf = 1-bandf*cis(ϕm-ϕf)
                    overlapmf *= 1/2
                    
                    m2 = overlapim*overlapmf/(ϵm-ϵi-ω-1im*δ)

                    if abs(ϵf-ϵi-2*ω)*histogram_width<1
                        #F2 +=  fi*(1-ff)*histogram_width/mesh^2*(abs(m1+m2))^2
                        F2 +=  (1/(4*π))*(μ/conesize)^2*fi*(1-ff)*k*histogram_width/mesh^2*(abs(m1+m2))^2
                    end
                end
            end
        end
    end
    F2 *= prefactor
    return F2
end