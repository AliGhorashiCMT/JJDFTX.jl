"""
$(TYPEDSIGNATURES)
Returns the electron phonon coupling, h(ϵ) as defined in Ab initio phonon coupling and optical response of hot electrons in plasmonic metals
"""
function ephcoupling(μ::Real, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, forcematrix::Array{Float64, 3}, 
    cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer; 
    mesh1::Integer=1, mesh2::Integer=1, histwidth1::Real=10, offset::Real=3, energyrange::Real=20, verbose::Bool=true)
# Note: Our densities of states will not include 1/Volume factors 
# Therefore, the six dimensional integral is just a 1/(Nk*Nk') times an average over the Brillouin zone. 
# We use a histogramming method for the product of the two delta functions. 
# We also have an extra factor of 4 compared to the cited paper due to the spin degeneracy that we take into account here
# However, this factor of 4 is cancelled out by the fac 
#2*gfermi/gϵ^2 
    deltaarray = zeros(histwidth1*energyrange, histwidth1*energyrange)
    ϵdosarray = zeros(histwidth1*energyrange)
    for _ in 1:mesh1^3
        k = rand(3)
        ϵks = wannier_bands(HWannier, cellmap, k, nbands) 
        for ϵk in ϵks 
            ϵdosarray[round(Int, (ϵk+offset)*histwidth1)+1] += histwidth1/mesh1^3
        end
    end
    gμ = ϵdosarray[round(Int, (μ+offset)*histwidth1)+1] #Density of states at fermi energy
    verbose && println("Density of States at Fermi Energy: ", gμ)
    prefactor = (2*gμ)*histwidth1*histwidth1*(1/mesh2^6)
    for _ in 1:mesh2^3
        k = rand(3) # Monte Carlo Sampling
        ϵks = wannier_bands(HWannier, cellmap, k, nbands) 
        for _ in 1:mesh2^3
            kprime = rand(3) # Monte Carlo Sampling
            q = kprime - k ## Phonon Wavevector
            phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
            ϵkprimes = wannier_bands(HWannier, cellmap, kprime, nbands) 
            ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
            for (b, ϵk) in enumerate(ϵks)
                for (bprime, ϵkprime) in enumerate(ϵkprimes)
                    for (α, phononomega) in enumerate(phononomegas)
                        #Note that energy scale of phonons is much smaller so we may set the phonon energy to be zero in the delta function 
                        deltaarray[round(Int, histwidth1*(ϵk+offset))+1, round(Int, histwidth1*(ϵkprime+offset))+1] += prefactor*phononomega*(1/ϵdosarray[round(Int, (ϵk+offset)*histwidth1)+1])*(1/ϵdosarray[round(Int, (ϵkprime+offset)*histwidth1)+1])*abs(ephmatrixelements[α, b, bprime])^2
                    end
                end
            end
        end
    end
    return deltaarray, np.einsum("ii->i", deltaarray)
end