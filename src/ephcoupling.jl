"""
$(TYPEDSIGNATURES)
Returns the electron phonon coupling, h(ϵ) as defined in Ab initio phonon coupling and optical response of hot electrons in plasmonic metals
"""
function ephcoupling(ϵ::Real, μ::Real, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, forcematrix::Array{Float64, 3}, 
    cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer; 
    mesh1::Integer=100, mesh2::Integer=100, histwidth1::Real=10, histwidth2::Real=100)
    hϵ = 0 
# Note: Our densities of states will not include 1/Volume factors 
# Therefore, the six dimensional integral is just a 1/(Nk*Nk') times an average over the Brillouin zone. 
# We use a histogramming method for the product of the two delta functions. 
# We also have an extra factor of 4 compared to the cited paper due to the spin degeneracy that we take into account here
# However, this factor of 4 is cancelled out by the fac 
#2*gfermi/gϵ^2 
    relevantks, subsamplingfraction = returnfermikpoint(HWannier, cellmap, nbands, ϵ, histwidth1, mesh=mesh1^3) ##Only look at ks where at least one band is close to ϵ
    nrelevantks = length(relevantks)
    println("Nks: ", nrelevantks )
    println("Subsampling Fraction:", subsamplingfraction)
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    gμ = dosatmu(HWannier, cellmap, nbands, μ, mesh=mesh1, histogram_width=histwidth1) 
    gϵ = dosatmu(HWannier, cellmap, nbands, ϵ, mesh=mesh1, histogram_width=histwidth1) 
    prefactor = (4*gμ/(gϵ)^2)*histwidth1*histwidth2*(1/mesh2^3)*(1/mesh2^3)*(subsamplingfraction)^2
    for _ in 1:mesh2^3
        k = relevantks[rand(1:nrelevantks)] # Monte Carlo Sampling
        ϵks = wannier_bands(HWannier, cellmap, k, nbands) 
        for _ in 1:mesh2^3
            kprime = relevantks[rand(1:nrelevantks)] # Monte Carlo Sampling
            q = kprime - k ## Phonon Wavevector
            phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
            ϵkprimes = wannier_bands(HWannier, cellmap, kprime, nbands) 
            ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
            for (b, ϵk) in enumerate(ϵks)
                for (bprime, ϵkprime) in enumerate(ϵkprimes)
                    for (α, phononomega) in enumerate(phononomegas)
                        (abs(ϵ-ϵk)*histwidth1 < 1  &&  abs(ϵkprime-ϵk-phononomega)*histwidth2 < 1) && (hϵ += prefactor*phononomega*abs(ephmatrixelements[α, b, bprime])^2)
                    end
                end
            end
        end
    end
    return hϵ
end