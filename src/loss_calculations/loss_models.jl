"""
$(TYPEDSIGNATURES)

Calculate the two plasmon absorption rate of graphene within the Dirac cone approximation. 

This returns the nondimensional F2(ω) as defined in the paper Multi-plasmon absorption in graphene by Jablan and Chang
"""
function graphenetwoplasmonemission(ω::Real, μ::Real; mesh::Integer = 100, histogram_width::Integer = 10, δ=1, conesize=1, verbose::Bool = true)
    F2 = 0 
    kF = μ/6
    q = 2π*2.5*137/(4*pi*ħ*c*μ)*ω^2 #ev^2/ev^2/angstrom -> inverse angstrom units
    # Effective dielectric constant of 2.5
    verbose && println(q)
    prefactor = (μ^3/kF^2)*(π^2)*(ω/(6*q))^6
    for (kiter, θiter) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        k = μ/conesize*kiter/mesh
        θ = 2π*θiter/mesh
        for bandi in [-1, 1] #Initial band may be in either lower or upper bands
            ϵi = 6*bandi*k
            fi = heaviside(μ-ϵi)
            fi ≈ 0 && continue
            for bandf in [1] #All electronic final states are in upper band
                ϵf = 6*bandf*sqrt((k*cos(θ)+2*q)^2+k^2*sin(θ)^2)
                ff = heaviside(μ-ϵf)
                (1-ff) ≈ 0 && continue 
                ϵm = 6*sqrt((k*cos(θ)+q)^2+k^2*sin(θ)^2)
                fm = 0
                ϕi = θ
                ϕm = atan((k*sin(θ))/(k*cos(θ)+q)) #We take q to lie on x axis WLOG 
                ϕf = atan((k*sin(θ))/(k*cos(θ)+2q))
                overlapim = 1/2*(1+bandi*cis(ϕi-ϕm))
                overlapmf = 1/2*(1+bandf*cis(ϕm-ϕf))
                m1 = (1-fm)*(overlapim*overlapmf)/(ϵm-ϵi-ω-1im*δ)
                ϵm = -6*sqrt((k*cos(θ)+q)^2+k^2*sin(θ)^2)
                fm = 0 
                overlapim = 1/2*(1-bandi*cis(ϕi-ϕm))
                overlapmf = 1/2*(1-bandf*cis(ϕm-ϕf))
                m2 = (1-fm)*overlapim*overlapmf/(ϵm-ϵi-ω-1im*δ)
                abs(ϵf-ϵi-2*ω)*histogram_width<0.5 || continue #Energy conserving delta function
                F2 +=  (1/(4*π^2))*(μ/conesize)*(2π)*fi*(1-ff)*k*histogram_width/mesh^2*(abs(m1+m2))^2
                #Explanation: the factors involving spin and valley degeneracy are included in the prefactor, so we have a 1/(2pi)^2 
                #factor from the conversion of the sum to a 2d brillouin zone integral. 
                #Note the interference between m1 and m2- which may add constructively or destructively in their contribution to the 
                #overall loss. 
            end
        end
    end
    F2 *= prefactor
    return F2
end
