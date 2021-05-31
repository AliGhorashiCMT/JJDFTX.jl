"""
$(TYPEDSIGNATURES)

Calculate the two plasmon absorption rate of graphene within the Dirac cone approximation. 

This returns the nondimensional F2(ω) as defined in the paper Multi-plasmon absorption in graphene by Jablan and Chang
"""
function graphenetwoplasmonemission(ω::Real, μ::Real; mesh = 100, histogram_width =10, δ=1, conesize=1, verbose::Bool = true)
    F2 = 0 
    kF = μ/6
    q = 2π*2.5*137/(4*pi*ħ*c*μ)*ω^2 #ev^2/ev^2/angstrom -> inverse angstrom units
    #Note that we are taking an effective dielectric constant of 2.5- corresponding to a silica layer below the graphene 
    #layer and vacuum on top of the graphene layer. 
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
                #fm = heaviside(μ-ϵm) #This should be the contribution, but we remain consistent with the paper cited
                fm = 0
                ϕi = θ
                ϕm = atan((k*sin(θ))/(k*cos(θ)+q)) #We take q to lie on x axis WLOG 
                ϕf = atan((k*sin(θ))/(k*cos(θ)+2q))
                overlapim = 1/2*(1+bandi*cis(ϕi-ϕm))
                overlapmf = 1/2*(1+bandf*cis(ϕm-ϕf))
                m1 = (1-fm)*(overlapim*overlapmf)/(ϵm-ϵi-ω-1im*δ)
                ϵm = -6*sqrt((k*cos(θ)+q)^2+k^2*sin(θ)^2)
                #fm = heaviside(μ-ϵm) #This should be the contribution, but we remain consistent with the paper cited
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

"""
$(TYPEDSIGNATURES)
"""
function graphenetwoplasmonemission(ω::Real, μ::Real, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}; diracpoint = 0, mesh = 100, offset::Vector{<:Real}=[2/3, -1/3, 0], subsample::Real=3, histogram_width =10, δ=1)
    F2 = 0 
    a= 1.42*sqrt(3)
    lat = [[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 20]]
    ucellarea = 5.238760872572826 #Graphene unit cell area in angstrom^2 
    kF = (μ-diracpoint)/6
    q= 2π*2.5*137/(4*pi*ħ*c*(μ-diracpoint))*ω^2 #ev^2/ev^2/angstrom -> inverse angstrom units
    qun = normalize_kvector(lat, (q, 0, 0))
    println(q)
    println(qun)
    prefactor = ((μ-diracpoint)^3/kF^2)*(π^2/ucellarea)*(ω/(6*q))^6
    for (xiter, yiter) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kx, ky = xiter/(subsample*mesh), yiter/(subsample*mesh) 
        ϵis = wannier_bands(HWannier, cellmap, [kx, ky, 0]+offset, 8)
        vis = wannier_vectors(HWannier, cellmap, [kx, ky, 0] +offset)

        ϵms = wannier_bands(HWannier, cellmap, [kx, ky, 0]+offset+ qun, 8)
        vms = wannier_vectors(HWannier, cellmap, [kx, ky, 0] +offset+ qun)

        ϵfs = wannier_bands(HWannier, cellmap, [kx, ky, 0]+offset+2*qun, 8)
        vfs = wannier_vectors(HWannier, cellmap, [kx, ky, 0]+offset+2*qun)

        for (indexi, ϵi) in enumerate(ϵis)
            vi = vis[:, indexi]
            fi = heaviside(μ-ϵi)
            fi ≈ 0 && continue
            for (indexf, ϵf) in enumerate(ϵfs)
                vf = vfs[:, indexf]
                ff = heaviside(μ-ϵf)
                (1-ff) ≈ 0 && continue 
                ifpair = 0 
                for (indexm, ϵm) in enumerate(ϵms)
                    fm = heaviside(μ-ϵm)
                    (1-fm) ≈ 0 && continue
                    vm = vms[:, indexm]
                    overlapim = np.dot(conj(vm), vi)
                    overlapmf = np.dot(conj(vf), vm)
                    m1 = (overlapim*overlapmf)/(ϵm-ϵi-ω-1im*δ)
                    ifpair += m1 
                end
                if abs(ϵf-ϵi-2*ω)*histogram_width<1
                    F2 +=  fi*(1-ff)*1/(subsample)^2*histogram_width/mesh^2*(abs(ifpair))^2
                end
            end
        end
    end
    F2 *= prefactor
    return F2
end