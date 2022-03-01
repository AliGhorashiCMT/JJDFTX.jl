
#The methods here are designed after "Plasmonics in Argentene by Shankar"

#Below, we try to calculate 1/τ(ω) as defined in cited paper above. 


"""
$(TYPEDSIGNATURES)

Returns the Fermi surface averaged velocity squared. Note, the current assumption is that the velocity is isotropic. 

This corresponds to the coefficient of 1/ω (ω in eV units) in the Drude formula
"""
function drude_conductivity(lattice::Vector{<:Vector{Float64}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, 
    PWannier::Array{Float64, 4}, nbands::Integer, μ::Real; mesh =10, histogram_width=10)
    area = unit_cell_area(lattice);
    σ = 0
    gs = 2
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        energies = wannier_bands(HWannier, cellmap, [xmesh/mesh, ymesh/mesh, 0],  nbands);
        np.any((energies .- μ)*histogram_width .< 0.5) || continue
        vnks =  imag(momentum_matrix_elements(HWannier, cellmap, PWannier, [xmesh/mesh, ymesh/mesh, 0])[1, :, :])
        for (n, energy) in enumerate(energies)
            abs(energy-μ)*histogram_width < 0.5 || continue  
            #println(n)
            vnk = vnks[n, n]/mₑ;
            σ += gs*(4*ħ^2)*histogram_width/(mesh^2)*(1/area)*abs(vnk)^2
        end
    end
    return σ
end


"""
$(TYPEDSIGNATURES)

To be used by function tauinverse. Returns a list of points in the brillouin zone with energy dispersions that intersect the Fermi surface. 
Also returns the subsampling fraction as the second return value. 
"""
function returnfermikpoint2d(HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, nbands::Integer, μ::Real, histogram_width::Real=10; mesh::Integer=1000)
    fermikpoints = Vector{Vector{Real}}()
    Nkfermi = 0 
    for _ in 1:mesh
        k = rand(2)
        energies = wannier_bands(HWannier, cellmap, [k..., 0], nbands)
        np.any(abs.((μ .- energies)*histogram_width) .< 0.5) || continue
        push!(fermikpoints, [k..., 0])
        Nkfermi += 1
    end
    return fermikpoints, Nkfermi/mesh
end

function btomega(ω::Real, fracroom::Real)
#fracroom is fraction of room finite_temperature_chemical_potential
    return ω/(1-exp(-ω*40/fracroom)) #Note that 1/40 eV is room temperature so β = 40 at room temperature
end

function tauinverse(HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, forcematrix::Array{<:Real, 3}, cellmapph::Array{<:Real, 2},
    heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, ωs::Vector{<:Real}, μ::Real, nbands::Integer; histogram_width::Integer=10, supplysampling::Union{Nothing, Tuple{<:Vector{<:Vector{<:Real}}, <:Real}}=nothing, 
    supplydos::Union{Nothing, Real}=nothing, mesh::Integer=10, mesh2::Integer=10000, fracroom::Real=0.3)
    tauinv = zeros(length(ωs))
    gμ = 0 
    dosmesh = 300
    if isnothing(supplydos)
        for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(dosmesh, dosmesh)))
            k = [xmesh/dosmesh, ymesh/dosmesh, 0]
            eks = wannier_bands(HWannier, cellmap, k, nbands)
            for ek in eks
                histogram_width*(abs(ek-μ)) < 0.5 || continue 
                gμ += histogram_width/dosmesh^2
            end
        end
    else
        gμ = supplydos
    end
    print("DOS at Fermi Energy is: ", gμ)
    kpoints, subsampling = isnothing(supplysampling) ? returnfermikpoint2d(HWannier, cellmap, nbands, μ, histogram_width, mesh=mesh2) : supplysampling
    nrelevantks = length(kpoints)
    println("Sampling from ", nrelevantks, "kpoints")

    for _ in 1:mesh
        k = kpoints[rand(1:nrelevantks)]
        eks = wannier_bands(HWannier, cellmap, k, nbands)
        vks= imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, k))
        for _ in 1:mesh
            kprime = kpoints[rand(1:nrelevantks)]
            q = kprime - k # Phonon wavevector
            phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
            ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
            ekprimes = wannier_bands(HWannier, cellmap, kprime, nbands)
            vkprimes= imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, kprime))
            for (n, ek) in enumerate(eks)
                vk = vks[1:2, n, n]
                vknorm = sqrt(sum(vk.^2))
                for (m, ekprime) in enumerate(ekprimes)
                    vkprime = vkprimes[1:2, m, m]
                    vkprimenorm = sqrt(sum(vk.^2))
                    for (α, phononomega) in enumerate(phononomegas)
                        (abs(ek-μ)*histogram_width < 0.5 && abs(ekprime-μ)*histogram_width < 0.5) || continue
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        gsquared = abs(ephmatrixelements[α, n, m])^2
                        tauinv += (btomega.(ωs .- phononomega, fracroom) ./ btomega.(ωs, fracroom)) * gsquared*velocityterm*histogram_width*histogram_width*(subsampling)^2/(mesh^2)
                    end
                end
            end
        end
    end
    return (2π/(ħ*gμ))*tauinv
end

"""
Interband conductivity as defined in the paper Plasmonics in argentene by Shankar. Instead of using the bloch eigenvectors, this
calculates the local conductivity by interpolating the momentum matrix elements. 

"""
function interbandsigma(lattice::Vector{<:Vector{Float64}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, 
    PWannier::Array{Float64, 4}, nbands::Integer, μ::Real; mesh::Integer =10, histogram_width::Integer=10, energy_range::Real=20)
    area = unit_cell_area(lattice)
    prefactor = 8*π*ħ^2/area
    mass = 0.5*1e6/(3*1e18)^2 #Mass of electron in eV/(angstrom^2/s^2)
    conds = zeros(round(histogram_width*energy_range) +1 )
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        energies = wannier_bands(HWannier, cellmap, [xmesh/mesh, ymesh/mesh, 0], nbands);
        vnks =  momentum_matrix_elements(HWannier, cellmap, PWannier, [xmesh/mesh, ymesh/mesh, 0])[1:2, :, :] #Only consider x and y components of momentum
        for (n, ϵ1) in enumerate(energies)
            f1 = ϵ1 > μ ? 1 : 0
            iszero(f1) && continue
            for (m, ϵ2) in enumerate(energies)
                f2 = ϵ2 < μ ? 1 : 0
                ω = ϵ1 - ϵ2
                iszero(f2) && continue
                ω < 0 && continue
                idx = round(Int, histogram_width*ω) + 1
                (idx > length(conds)) && continue
                vk = vnks[:, n, m]/mass
                velocityterm = sum((abs.(vk)).^2)/2 #Assume isotropy. 
                conds[idx] += prefactor*histogram_width/(mesh^2)*(1/ω)*velocityterm
            end
        end
    end
    return conds
end