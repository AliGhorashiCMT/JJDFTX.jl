"""Returns the eliashberg spectral function. This function is modeled after http://jdftx.org/EphMatrixElements.html
"""
function eliashberg(lattice::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, PWannier::Array{Float64, 4}, forcematrix::Array{Float64, 3}, cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer,  μ::Real; mesh::Integer=10, histogram_width::Real=10, energyrange::Real=1)
#=We have to sum over the brillouin zone twice, and over two electronic band indices and one phonon band index
We have three delta functions. One which enforces the frequency to be equal to the phonon energy, 
=#

#= Units check: 
=#
    omegas = zeros(Int(energyrange*histogram_width))
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    println("Number of phonons: ", nphononmodes)
    println("Number of electron bands: ", nbands)
    gs = 1 ## In our DOS function we don't take spin into account 
    gμ = dosatmu(HWannier, cellmap, lattice, nbands, μ) ##Density of states at fermi level (in units of 1/eV*1/angstrom^3)
    for kiter in 1:mesh
        k = rand(3) # Monte Carlo sampling
        eks = wannier_bands(HWannier, cellmap, k, nbands)
        vks = abs.(momentum_matrix_elements(HWannier, cellmap, PWannier, k))
        for kprimiter in 1:mesh
            kprime = rand(3) # Monte Carlo sampling
            vkprimes = abs.(momentum_matrix_elements(HWannier, cellmap, PWannier, kprime))
            q = kprime - k ## Phonon Wavevector
            ekprimes = wannier_bands(HWannier, cellmap, kprime, nbands)
            phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
            ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
            for b in 1:nbands
                ek = eks[b]
                vk = vks[:, b, b]
                vknorm = sqrt(sum(vk.^2))
                for bprime in 1:nbands
                    ekprime = ekprimes[bprime]
                    vkprime = vkprimes[:, bprime, bprime]
                    vkprimenorm = sqrt(sum(vkprime.^2))
                    for α in 1:nphononmodes
                        phononomega = phononomegas[α]
                        velocityterm = (1-dot(vk, vkprime)/(vknorm*vkprimenorm))
                        omegas[round(Int, phononomega*histogram_width)+1]  += (gs/gμ)^2*abs(ephmatrixelements[α, b, bprime])^2*(1/π)^2*imag(1/((ek-μ)+1im))*imag(1/((ekprime-μ)+1im))*velocityterm*1/mesh^2*histogram_width # Use Lorentzian representation of delta function 
                    end
                end
            end
        end
    end
    return omegas
end

"For use by the Eliashberg spectral function method above"
function dosatmu(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice::Vector{<:Vector{<:Real}}, nbands::Integer, μ::Real; mesh::Integer = 10, histogram_width::Real=3)
    volume = unit_cell_volume(lattice)
    dos = 0 
    for x_mesh in 1:mesh^3
        ϵs = wannier_bands(Hwannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            if abs(μ-ϵ)*histogram_width < 1
                dos = dos + histogram_width*(1/mesh)^3*(1/volume)
            end
        end
    end
    return dos
end

function vFsquaredatmu(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, nbands::Integer, μ::Real; mesh::Integer=10, histogram_width::Real=3)
    vFsquared = 0 
    numintersections = 0 
    for x_mesh in 1:mesh^3
        arandk = rand(3)
        ϵs = wannier_bands(Hwannier, cell_map, arandk, nbands)
        ps = momentum_matrix_elements(Hwannier, cell_map, Pwannier, arandk)
        for (index, ϵ) in enumerate(ϵs)
            if abs(μ-ϵ)*histogram_width < 1
                numintersections +=1
                vFsquared = vFsquared + sum((abs.(ps[:, index, index])).^2)*(bohrtoangstrom/ħ)^2
                ##Note that to stay in keeping with JDFTX conventions, we reconverted to atomic units
            end
        end
    end
    averaged_fermivelocity = sqrt(vFsquared/numintersections)
    isnan(averaged_fermivelocity)  ? println("Got NaN- which typically means you need more sampling points. Try increasing mesh.") : nothing
    return sqrt(vFsquared/numintersections)
end

function convertdos(dos::Real)
    ##Conventions of this package are that the dos will be in 1/angstrom^3*1/eV units, to convert to jdftx units, 
    ##we must do the following:
    return dos*bohrtoangstrom^3*1/eV
end