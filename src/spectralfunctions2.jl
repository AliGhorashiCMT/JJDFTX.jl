"""
$(TYPEDSIGNATURES)
This function returns kpoints at the Fermi level and also returns the subsampling.
"""
function returnfermikpoint(HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, nbands::Integer, μ::Real, histogram_width::Real=10; mesh::Integer=1000)
    fermikpoints = Vector{Vector{Real}}()
    Nkfermi = 0 
    for _ in 1:mesh
        k=rand(3)
        energies = wannier_bands(HWannier, cellmap, k, nbands)
        atFermi = false
        for energy in energies
            abs(μ-energy)*histogram_width<1 ? atFermi=true : nothing
        end
        if atFermi==true
            push!(fermikpoints, k)
            Nkfermi += 1
        end
    end
    return fermikpoints, Nkfermi/mesh
end

"""
$(TYPEDSIGNATURES)
"""
function returnfermikpoint_lorentzian(HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, nbands::Integer, μ::Real; esmearing::Real=1, mesh::Integer=1000)
    fermikpoints = Vector{Vector{Real}}()
    Nkfermi = 0 
    for _ in 1:mesh
        k=rand(3)
        energies = wannier_bands(HWannier, cellmap, k, nbands)
        atFermi = false
        for energy in energies
            abs(1/π*imag(1/(energy-μ+esmearing*1im))*esmearing) > 0.001 ? atFermi=true : nothing
        end
        if atFermi==true
            push!(fermikpoints, k)
            Nkfermi += 1
        end
    end
    return fermikpoints, Nkfermi/mesh
end

"""
$(TYPEDSIGNATURES)
"""
function returnfermikpoint_gaussian(HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, nbands::Integer, μ::Real; esmearing::Real=1, mesh=1000)
    fermikpoints = Vector{Vector{Real}}()
    Nkfermi = 0 
    for _ in 1:mesh
        k=rand(3)
        energies = wannier_bands(HWannier, cellmap, k, nbands)
        atFermi = false
        for energy in energies
            1/(esmearing*sqrt(2π))*exp(-0.5*((energy-μ)/esmearing)^2)*esmearing > 0.1 ? atFermi=true : nothing
        end
        if atFermi==true
            push!(fermikpoints, k)
            Nkfermi += 1
        end
    end
    return fermikpoints, Nkfermi/mesh
end

#TODO: Test this function
"""
$(TYPEDSIGNATURES)
Much faster version of eliashberg2 and eliashberg. The purpose of this function is to only look at relevant k points near the Fernu energy. This is substantially better in higher dimensions. For the eliashberg function in 3d, this translates to 6d monte carlo integrations
"""
function eliashberg3(lattice::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, forcematrix::Array{Float64, 3}, cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer, μ::Real; mesh::Integer=10, histogram_width::Real=10, histogram_width2::Real=3, energyrange::Real=1)
    #Find the relevant k points near the Fermi energy 
    relevantks, subsamplingfraction = returnfermikpoint(HWannier, cellmap, nbands, μ, histogram_width2, mesh=50^3) ##Sample 1 million points
    nrelevantks = length(relevantks)

    omegas = zeros(Int(energyrange*histogram_width))
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    println("Number of phonons: ", nphononmodes)
    println("Number of electron bands: ", nbands)
    println("Number of relevant kpoints (make sure this is relatively large...I don't know use your judgement): ", nrelevantks)

    gs = 1 ## In our DOS function we don't take spin into account 
    gμ = dosatmu(HWannier, cellmap, lattice, nbands, μ) ##Density of states at fermi level (in units of 1/eV*1/angstrom^3)
    for _ in 1:mesh ##Sample over mesh number of initial kvectors
        k = relevantks[rand(1:nrelevantks)] # Monte Carlo sampling. Choose an index of a k point at the Fermi level 
        eks = wannier_bands(HWannier, cellmap, k, nbands) 
        vks = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, k)) ##Find momentum matrix elements for k 
        for _ in 1:mesh #Sample over mesh number of initial kvectors
            kprime = relevantks[rand(1:nrelevantks)] # Monte Carlo sampling
            vkprimes = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, kprime))
            q = kprime - k ## Phonon Wavevector
            ekprimes = wannier_bands(HWannier, cellmap, kprime, nbands)
            phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
            ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
            for b in 1:nbands
                ek = eks[b] 
                vk = vks[:, b, b] #Only look at diagonal velocity components  
                vknorm = sqrt(sum(vk.^2))
                for bprime in 1:nbands
                    ekprime = ekprimes[bprime]
                    vkprime = vkprimes[:, bprime, bprime]
                    vkprimenorm = sqrt(sum(vkprime.^2))
                    for α in 1:nphononmodes
                        phononomega = phononomegas[α]
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        abs(ek-μ)*histogram_width2<1 && abs(ekprime-μ)*histogram_width2<1 ? omegas[round(Int, phononomega*histogram_width)+1]  += (gs/gμ)^2*abs(ephmatrixelements[α, b, bprime])^2*histogram_width2*histogram_width2*velocityterm*1/mesh^2*histogram_width : nothing # Use Lorentzian representation of delta function 
                    end
                end
            end
        end
    end
    return omegas*subsamplingfraction*subsamplingfraction #Because we only looked at Fermi kvectors and not arbitrary kvectors   
end

"""
$(TYPEDSIGNATURES)
Same as eliashberg3 but with different representation of delta function
"""
function eliashberg4(lattice::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, forcematrix::Array{Float64, 3}, cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer, μ::Real; mesh::Integer=10, esmearing::Real=.005, histogram_width::Real=1000, energyrange::Real=1)
    #Find the relevant k points near the Fermi energy 
    relevantks, subsamplingfraction = returnfermikpoint_gaussian(HWannier, cellmap, nbands, μ, esmearing=esmearing, mesh=60^3)
    nrelevantks = length(relevantks)

    omegas = zeros(Int(energyrange*histogram_width))
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    println("Number of phonons: ", nphononmodes)
    println("Number of electron bands: ", nbands)
    println("Number of relevant kpoints (make sure this is relatively large...I don't know use your judgement): ", nrelevantks)

    gs = 1 ## In our DOS function we don't take spin into account 
    #gμ = dosatmu(HWannier, cellmap, lattice, nbands, μ) ##Density of states at fermi level (in units of 1/eV*1/angstrom^3)
    gμ = dosatmugaussian(HWannier, cellmap, lattice, nbands, μ, esmearing=esmearing, mesh=30) # Density of states has to have a similar meshing as fermi level sampling
    for _ in 1:mesh ##Sample over mesh number of initial kvectors
        k = relevantks[rand(1:nrelevantks)] # Monte Carlo sampling. Choose an index of a k point at the Fermi level 
        eks = wannier_bands(HWannier, cellmap, k, nbands) 
        vks = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, k))
        for _ in 1:mesh #Sample over mesh number of initial kvectors
            kprime = relevantks[rand(1:nrelevantks)] # Monte Carlo sampling
            vkprimes = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, kprime))
            q = kprime - k ## Phonon Wavevector
            ekprimes = wannier_bands(HWannier, cellmap, kprime, nbands)
            phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
            ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
            for b in 1:nbands
                ek = eks[b] 
                vk = vks[:, b, b] #Only look at diagonal velocity components  
                vknorm = sqrt(sum(vk.^2))
                for bprime in 1:nbands
                    ekprime = ekprimes[bprime]
                    vkprime = vkprimes[:, bprime, bprime]
                    vkprimenorm = sqrt(sum(vkprime.^2))
                    for α in 1:nphononmodes
                        phononomega = phononomegas[α]
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        omegas[round(Int, phononomega*histogram_width)+1]  += (gs/gμ)^2*abs(ephmatrixelements[α, b, bprime])^2*1/(2*pi*esmearing^2)*exp(-0.5*((ek-μ)/esmearing)^2-0.5*((ekprime-μ)/esmearing)^2)*velocityterm*1/mesh^2*histogram_width # Use Lorentzian representation of delta function 
                    end
                end
            end
        end
    end
    return omegas*subsamplingfraction*subsamplingfraction #Because we only looked at Fermi kvectors and not arbitrary kvectors    #*subsampling2(HWannier, cellmap, nbands, μ, histogram_width2)^2
end

"""
$(TYPEDSIGNATURES)
Same as eliashberg3 but with different representation of delta function
"""
function eliashberg5(lattice::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, forcematrix::Array{Float64, 3}, cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer, μ::Real; mesh::Integer=10, esmearing::Real=.005, histogram_width::Real=1000, energyrange::Real=1)
    println("Calculating Eliashberg spectral function with Lorentzian representation. If this is not what you want, consider eliashberg3 and 4")
    #Find the relevant k points near the Fermi energy 
    relevantks, subsamplingfraction = returnfermikpoint_lorentzian(HWannier, cellmap, nbands, μ, esmearing=esmearing, mesh=20^3)
    nrelevantks = length(relevantks)
    omegas = zeros(Int(energyrange*histogram_width))
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    println("Number of phonons: ", nphononmodes)
    println("Number of electron bands: ", nbands)
    println("Number of relevant kpoints (make sure this is relatively large...I don't know use your judgement): ", nrelevantks)

    gs = 1 ## In our DOS function we don't take spin into account 
    #gμ = dosatmu(HWannier, cellmap, lattice, nbands, μ) ##Density of states at fermi level (in units of 1/eV*1/angstrom^3)
    gμ = dosatmulorentzian(HWannier, cellmap, lattice, nbands, μ, esmearing=esmearing, mesh=30) # In order to get correct results, the density of states must be calculated the same way as the fermi level sampling
    for _ in 1:mesh ##Sample over mesh number of initial kvectors
        k = relevantks[rand(1:nrelevantks)] # Monte Carlo sampling. Choose an index of a k point at the Fermi level 
        eks = wannier_bands(HWannier, cellmap, k, nbands) 
        vks = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, k))
        for _ in 1:mesh #Sample over mesh number of initial kvectors
            kprime = relevantks[rand(1:nrelevantks)] # Monte Carlo sampling
            vkprimes = imag.(momentum_matrix_elements(HWannier, cellmap, PWannier, kprime))
            q = kprime - k ## Phonon Wavevector
            ekprimes = wannier_bands(HWannier, cellmap, kprime, nbands)
            phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
            ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
            for b in 1:nbands
                ek = eks[b] 
                vk = vks[:, b, b] #Only look at diagonal velocity components  
                vknorm = sqrt(sum(vk.^2))
                for bprime in 1:nbands
                    ekprime = ekprimes[bprime]
                    vkprime = vkprimes[:, bprime, bprime]
                    vkprimenorm = sqrt(sum(vkprime.^2))
                    for α in 1:nphononmodes
                        phononomega = phononomegas[α]
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        omegas[round(Int, phononomega*histogram_width)+1]  += (gs/gμ)^2*abs(ephmatrixelements[α, b, bprime])^2*(1/π)^2*imag(1/(ek-μ+esmearing*1im))*imag(1/(ekprime-μ+esmearing*1im))*velocityterm*1/mesh^2*histogram_width # Use Lorentzian representation of delta function 
                    end
                end
            end
        end
    end
    return omegas*subsamplingfraction*subsamplingfraction #Because we only looked at Fermi kvectors and not arbitrary kvectors    #*subsampling2(HWannier, cellmap, nbands, μ, histogram_width2)^2
end



