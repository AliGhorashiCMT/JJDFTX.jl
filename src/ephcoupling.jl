"""
$(TYPEDSIGNATURES)
Returns the electron phonon coupling, h(ϵ) as defined in Ab initio phonon coupling and optical response of hot electrons in plasmonic metals
"""
function ephcoupling(ϵ::Real, μ::Real, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, forcematrix::Array{Float64, 3}, 
    cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer; 
    mesh1::Integer=100, mesh2::Integer=100, histwidth1::Real=10)
    hϵ = 0 
# Note: Our densities of states will not include 1/Volume factors 
# Therefore, the six dimensional integral is just a 1/(Nk*Nk') times an average over the Brillouin zone. 
# We use a histogramming method for the product of the two delta functions. 
# We also have an extra factor of 4 compared to the cited paper due to the spin degeneracy that we take into account here
# However, this factor of 4 is cancelled out by the fac 
#2*gfermi/gϵ^2 
    relevantks = Vector{Vector{Real}}()
    Nkfermi = 0 
    gμ = 0 
    gϵ = 0 
    for _ in 1:mesh1^3
        k=rand(3)
        energies = wannier_bands(HWannier, cellmap, k, nbands)
        atFermi = false
        for energy in energies
            abs(energy-ϵ)*histwidth1<1 && (atFermi=true)
            abs(energy-ϵ)*histwidth1<1 && (gϵ += histwidth1/mesh1^3)
            abs(energy-μ)*histwidth1<1 && (gμ += histwidth1/mesh1^3)
        end
        if atFermi==true
            push!(relevantks, k)
            Nkfermi += 1
        end
    end
    subsamplingfraction = Nkfermi/mesh1^3
    nrelevantks = length(relevantks)
    println("Nks: ", nrelevantks )
    println("Subsampling Fraction:", subsamplingfraction)
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    prefactor = (2*gμ/(gϵ)^2)*histwidth1*histwidth1*(1/mesh2^3)*(1/mesh2^3)*(subsamplingfraction)^2
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
                        #Note that energy scale of phonons is much smaller so we may set the phonon energy to be zero in the delta function 
                        (abs(ϵ-ϵk)*histwidth1 < 1  &&  abs(ϵkprime-ϵk)*histwidth1 < 1) && (hϵ += prefactor*phononomega*abs(ephmatrixelements[α, b, bprime])^2)
                        #(abs(ϵ-ϵk)*histwidth1 < 1  &&  abs(ϵkprime-ϵk-phononomega)*histwidth2 < 1) && (hϵ += prefactor*phononomega*abs(ephmatrixelements[α, b, bprime])^2)
                    end
                end
            end
        end
    end
    return hϵ
end


function ephcoupling(ϵs::Vector{<:Real}, μ::Real, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, forcematrix::Array{Float64, 3}, 
    cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer; 
    mesh1::Integer=100, mesh2::Integer=100, histwidth1::Real=10)
    hϵ = zeros(length(ϵs)) 
# Note: Our densities of states will not include 1/Volume factors 
# Therefore, the six dimensional integral is just a 1/(Nk*Nk') times an average over the Brillouin zone. 
# We use a histogramming method for the product of the two delta functions. 
# We also have an extra factor of 4 compared to the cited paper due to the spin degeneracy that we take into account here
# However, this factor of 4 is cancelled out by the fac 
#2*gfermi/gϵ^2  d
    relevantks = Dict{Real, Vector{Vector{Real}}}() #Dictionary of all k points that intersect the key energy
    NkFermis = Dict{Real, Real}()  # Dictionary of number of kpoints at the key energy
    Doses = Dict{Real, Real}()
    for ϵ in ϵs
        push!(relevantks, ϵ => [Real[]])
        push!(NkFermis, ϵ => 0)
        push!(Doses, ϵ => 0 )
    end
    gμ = 0 
    for _ in 1:mesh1^3
        k=rand(3)
        energies = wannier_bands(HWannier, cellmap, k, nbands)
        for energy in energies
            atEnergies = filter(ϵ->abs(energy-ϵ)*histwidth1<1, ϵs)
            for atEnergy in atEnergies
                push!(relevantks[atEnergy], k)
                Doses[atEnergy] += histwidth1/mesh1^3
            end
            abs(energy-μ)*histwidth1<1 && (gμ += histwidth1/mesh1^3)
        end
    end
    for ϵ in ϵs
        filter!(!isempty, unique!(relevantks[ϵ])) # Get rid of redundant k vectors and only consider ones that are not empty
        NkFermis[ϵ] = length(relevantks)/mesh1^3
    end
    println("Subsamplings : ", NkFermis)
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    for (index, ϵ) in enumerate(ϵs)
        gϵ = Doses[ϵ]
        subsamplingfraction = NkFermis[ϵ]
        relevantk = relevantks[ϵ]
        nrelevantk = length(relevantk)
        prefactor = (2*gμ/(gϵ)^2)*histwidth1*histwidth1*(1/mesh2^3)*(1/mesh2^3)*(subsamplingfraction)^2
        for _ in 1:mesh2^3
            k = relevantk[rand(1:nrelevantk)] # Monte Carlo Sampling
            ϵks = wannier_bands(HWannier, cellmap, k, nbands) 
            for _ in 1:mesh2^3
                kprime = relevantk[rand(1:nrelevantk)] # Monte Carlo Sampling
                q = kprime - k ## Phonon Wavevector
                phononomegas = phonon_dispersion(forcematrix, cellmapph, q)
                ϵkprimes = wannier_bands(HWannier, cellmap, kprime, nbands) 
                ephmatrixelements = eph_matrix_elements(heph, cellmapeph, forcematrix, cellmapph, HWannier, cellmap, k, kprime, nbands)
                for (b, ϵk) in enumerate(ϵks)
                    for (bprime, ϵkprime) in enumerate(ϵkprimes)
                        for (α, phononomega) in enumerate(phononomegas)
                            #Note that energy scale of phonons is much smaller so we may set the phonon energy to be zero in the delta function 
                            (abs(ϵ-ϵk)*histwidth1 < 1  &&  abs(ϵkprime-ϵk)*histwidth1 < 1) && (hϵ[index] += prefactor*phononomega*abs(ephmatrixelements[α, b, bprime])^2)
                            #(abs(ϵ-ϵk)*histwidth1 < 1  &&  abs(ϵkprime-ϵk-phononomega)*histwidth2 < 1) && (hϵ += prefactor*phononomega*abs(ephmatrixelements[α, b, bprime])^2)
                        end
                    end
                end
            end
        end
    end
    return hϵ
end

"""
$(TYPEDSIGNATURES)
Returns the electron phonon coupling, h(ϵ) as defined in Ab initio phonon coupling and optical response of hot electrons in plasmonic metals
"""
function ephcoupling2(ϵ::Real, μ::Real, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, forcematrix::Array{Float64, 3}, 
    cellmapph::Array{Float64, 2}, heph::Array{Float64, 5}, cellmapeph::Array{<:Real, 2}, nbands::Integer; 
    mesh1::Integer=1, mesh2::Integer=1, histwidth1::Real=10, offset::Real=3, energyrange::Real=20)
# Note: Our densities of states will not include 1/Volume factors 
# Therefore, the six dimensional integral is just a 1/(Nk*Nk') times an average over the Brillouin zone. 
# We use a histogramming method for the product of the two delta functions. 
# We also have an extra factor of 4 compared to the cited paper due to the spin degeneracy that we take into account here
# However, this factor of 4 is cancelled out by the fac 
#2*gfermi/gϵ^2 
    ϵarray = zeros(histwidth1*energyrange)
    ϵdosarray = zeros(histwidth1*energyrange)
    nphononmodes = length(phonon_dispersion(forcematrix, cellmapph, [0, 0, 0]))
    for _ in 1:mesh1^3
        k = rand(3)
        ϵks = wannier_bands(HWannier, cellmap, k, nbands) 
        for ϵk in ϵks 
            ϵdosarray[round(Int, (ϵk+offset)*histwidth1)+1]+= histwidth1/mesh1^3
        end
    end
    gμ = ϵdosarray[round(Int, (μ+offset)*histwidth1)]
    prefactor = (2*gμ)*histwidth1*histwidth1*(1/mesh2^3)*(1/mesh2^3)
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
                        (abs(ϵkprime-ϵk)*histwidth1 < 1) && (ϵarray[round(Int, histwidth1*(ϵ+offset))+1] += (prefactor*ϵdosarray[round(Int, histwidth1*(ϵ+offset))+1]*phononomega*abs(ephmatrixelements[α, b, bprime])^2))
                        #(abs(ϵ-ϵk)*histwidth1 < 1  &&  abs(ϵkprime-ϵk-phononomega)*histwidth2 < 1) && (hϵ += prefactor*phononomega*abs(ephmatrixelements[α, b, bprime])^2)
                    end
                end
            end
        end
    end
    return ϵarray
end