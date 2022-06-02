"""
$(TYPEDSIGNATURES)
This function returns kpoints at the Fermi level and also returns the subsampling.
This is designed to be used to pass relevant kpoints to the eliashberg spectral function methods. When evaluating formula that involve electron-phonon matrix elements, 
an extremely dense sampling of the Brillouin zone is required. In order to circumvent this issue, quantities that primarily rely on the Fermi surface can be preprocessed
in order to not needlessly sample over regions of the Brillouin zone that give a null contribution. 
"""
function returnfermikpoint(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, μ::Real, ::Val{D}=Val(2); 
    weight::Union{Val{:gaussian}, Val{:lorentzian}, Val{:histogram}} = Val(:histogram),
    histogram_width::Real=10, mesh::Integer=10, num_blocks::Integer=10, esmearing::Real=1) where D

    fermikpoints = Array{Real, 2}[]
    Nkfermi = 0 
    for _ in 1:num_blocks
        kpoints = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        Es, _ = wannier_bands(Hwannier, cell_map, kpoints)
        Fermi_Surface = 
            if weight == Val(:histogram)
                np.einsum("ij -> i", abs.((μ.- Es)*histogram_width) .< 0.5)
            elseif weight == Val(:gaussian)
                np.einsum("ij -> i", exp.(-0.5*((Es .- μ)/esmearing).^2) .> sqrt(2*π)*0.001)
            elseif weight == Val(:lorentzian)
                np.einsum("ij -> i", abs.(-1/π*imag.(esmearing ./ (Es .- μ .+ esmearing*1im))) .> 0.001)
            end
        push!(fermikpoints, kpoints[:, Fermi_Surface])
        Nkfermi += sum(Fermi_Surface)
    end
    return hcat(fermikpoints...), Nkfermi/mesh^D/num_blocks
end

"""
$(TYPEDSIGNATURES)
Much faster version of eliashberg2 and eliashberg. The purpose of this function is to only look at relevant k points near the Fermi energy. This is substantially better in higher dimensions. For the eliashberg function in 3d, this translates to 6d monte carlo integrations
"""
function eliashberg(lattice::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, 
    forcematrix::Array{Float64, 3}, cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer, μ::Real; 
    mesh::Integer=10, histogram_width::Real=10, histogram_width2::Real=3, energyrange::Real=1, verbose::Bool=true, subsamp::Union{Tuple{<:Vector{<:Vector{<:Real}}, Real}, Nothing}=nothing, dosmu::Union{Real, Nothing}=nothing)

    #Find the relevant k points near the Fermi energy 
    #We offer the user the chance to give the fermikpoints initially- however for simple tests it might be easier to just let the function go through
    #the whole process of calculating the spectral function
    relevantks, subsamplingfraction = isnothing(subsamp) ? returnfermikpoint(HWannier, cellmap, nbands, μ, histogram_width2, mesh=10^3) : subsamp##Sample 1 million points
    nrelevantks = length(relevantks)

    omegas = zeros(Int(energyrange*histogram_width))
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    verbose && println("Number of phonons: ", nphononmodes)
    verbose && println("Number of electron bands: ", nbands)
    verbose && println("Number of relevant kpoints (make sure this is relatively large...I don't know use your judgement): ", nrelevantks)
    verbose && println("Subsampling fraction is: ", subsamplingfraction )
    gs = 1 ## In our DOS function we don't take spin into account 
    gμ = isnothing(dosmu) ? dosatmu(HWannier, cellmap, lattice, nbands, μ, mesh=40) : dosmu##Density of states at fermi level (in units of 1/eV*1/angstrom^3)
    verbose && println("Density of States is: ", gμ)
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
            for (band1, ek) in enumerate(eks)
                vk = vks[:, band1, band1] #Only look at diagonal velocity components  
                vknorm = sqrt(sum(vk.^2))
                for (band2, ekprime) in enumerate(ekprimes)
                    vkprime = vkprimes[:, band2, band2]
                    vkprimenorm = sqrt(sum(vkprime.^2))
                    for (α, phononomega) in enumerate(phononomegas)
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        (abs(ek-μ)*histogram_width2 < 0.5 && abs(ekprime-μ)*histogram_width2 < 0.5) || continue
                        omegas[round(Int, phononomega*histogram_width)+1]  += (gs/gμ)^2*abs(ephmatrixelements[α, band1, band2])^2*histogram_width2*histogram_width2*velocityterm*1/mesh^2*histogram_width  # Use Lorentzian representation of delta function 
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
function eliashberg_gaussian(lattice::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, 
    forcematrix::Array{Float64, 3}, cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer, μ::Real; mesh::Integer=10, 
    esmearing::Real=.005, histogram_width::Real=1000, energyrange::Real=1)

    #Find the relevant k points near the Fermi energy 
    relevantks, subsamplingfraction = returnfermikpoint_gaussian(HWannier, cellmap, nbands, μ, esmearing=esmearing, mesh=60^3)
    nrelevantks = length(relevantks)

    omegas = zeros(Int(energyrange*histogram_width))
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    println("Number of phonons: ", nphononmodes)
    println("Number of electron bands: ", nbands)
    println("Number of relevant kpoints (make sure this is relatively large...I don't know use your judgement): ", nrelevantks)

    gs = 1 ## In our DOS function we don't take spin into account 
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
            for (band1, ek) in enumerate(eks)
                vk = vks[:, band1, band1] #Only look at diagonal velocity components  
                vknorm = sqrt(sum(vk.^2))
                for (band2, ekprime) in enumerate(ekprimes)
                    vkprime = vkprimes[:, band2, band2]
                    vkprimenorm = sqrt(sum(vkprime.^2))
                    for (α, phononomega) in enumerate(phononomegas)
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        omegas[round(Int, phononomega*histogram_width)+1]  += (gs/gμ)^2*abs(ephmatrixelements[α, band1, band2])^2*1/(2*pi*esmearing^2)*exp(-0.5*((ek-μ)/esmearing)^2-0.5*((ekprime-μ)/esmearing)^2)*velocityterm*1/mesh^2*histogram_width # Use Lorentzian representation of delta function 
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
function eliashberg_lorentzian(lattice::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, forcematrix::Array{Float64, 3}, cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer, μ::Real; mesh::Integer=10, esmearing::Real=.005, histogram_width::Real=1000, energyrange::Real=1)
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
            for (band1, ek) in enumerate(eks)
                vk = vks[:, band1, band1] #Only look at diagonal velocity components  
                vknorm = sqrt(sum(vk.^2))
                for (band2, ekprime) in enumerate(ekprimes)
                    vkprime = vkprimes[:, band2, band2]
                    vkprimenorm = sqrt(sum(vkprime.^2))
                    for (α, phononomega) in enumerate(phononomegas)
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        omegas[round(Int, phononomega*histogram_width)+1]  += (gs/gμ)^2*abs(ephmatrixelements[α, b, bprime])^2*(1/π)^2*imag(1/(ek-μ+esmearing*1im))*imag(1/(ekprime-μ+esmearing*1im))*velocityterm*1/mesh^2*histogram_width # Use Lorentzian representation of delta function 
                    end
                end
            end
        end
    end
    return omegas*subsamplingfraction*subsamplingfraction #Because we only looked at Fermi kvectors and not arbitrary kvectors    #*subsampling2(HWannier, cellmap, nbands, μ, histogram_width2)^2
end



