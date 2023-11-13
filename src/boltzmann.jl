
#The methods here are designed after "Plasmonics in Argentene by Shankar"

#Below, we try to calculate 1/τ(ω) as defined in cited paper above. 


"""
$(TYPEDSIGNATURES)

Returns the Fermi surface averaged velocity squared. Returns a one and a 3 dimensional array of which the first
index corresponds to the chemical potential, and the second and third indices correspond to cartesian directions.

This corresponds to the coefficient of 1/ω (ω in eV units) in the Drude formula in units of e²/(4ħ)

"""
function drude_conductivity(lattice_vectors::Vector{<:Vector{Float64}}, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    Pwannier::Array{Float64, 4}, ::Val{D}=Val(2); mesh::Integer =10, num_blocks::Integer=10, histogram_width::Integer=10, degeneracy::Integer=1) where D
    V = D == 2 ? unit_cell_area(lattice_vectors) : unit_cell_volume(lattice_vectors) 
    σ =  zeros(200*histogram_width, 3, 3)
    gs = degeneracy
    μarray = collect(range(-100, 100, length=200*histogram_width))
    μarray = np.reshape(μarray, (-1, 1, 1))
    μarray = np.repeat(np.repeat(μarray, mesh^D, axis=1), size(Hwannier)[2], axis=2)
    for i in 1:num_blocks
        println("Block: $i"); flush(stdout)
        ks = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        Eks, Uks= wannier_bands(Hwannier, cell_map, ks);
        Eks = np.repeat(np.reshape(Eks, (1, mesh^D, size(Hwannier)[2])), 200*histogram_width, axis=0)
        vks = np.einsum("kpaa -> kpa ", imag.(momentum_matrix_elements(Uks, cell_map, Pwannier, ks)))  
        σ += np.einsum("ekn, knml -> eml", abs.(Eks .- μarray)*histogram_width .< 0.5, np.einsum("kmn, kln->knml", vks/mₑ, vks/mₑ))*1/mesh^D*gs*(4*ħ^2)*histogram_width*(1/V)
    end
    return collect(range(-100, 100, length=200*histogram_width)), σ/num_blocks
end

btomega(ω::Real, fracroom::Real) = ω/(1-exp(-ω/(kB*295*fracroom))) 

fermi(ω::Real, fracroom::Real) = 1/(exp(ω/(kB*295*fracroom)) + 1) 

bose(ω::Real, fracroom::Real) = 1/(exp(ω/(kB*295*fracroom)) - 1)

function τ(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2},
    Heph::Array{Float64, 5}, celleph_map::Array{<:Real, 2}, ωs::Vector{<:Real}, μ::Real, 
    weight::Union{Val{:gaussian}, Val{:lorentzian}, Val{:histogram}} = Val(:histogram), ::Val{D}=Val(2); 
    histogram_width::Integer=10, esmearing::Real=10, supplysampling::Union{Tuple{<:Matrix{<:Real}, Real}, Nothing}=nothing, 
    supplydos::Union{Nothing, Real}=nothing, mesh::Integer=10, num_blocks::Integer=10, dosmesh::Integer = 10, dos_num_blocks=3, fracroom::Real=1, intraband::Union{Nothing, Integer}=nothing) where D
    
    tauinv = zeros(length(ωs))

    gμ = 
    if !isnothing(supplydos)
            supplydos 
    else
        Es, DOS = density_of_states(Hwannier, cell_map, Val(D); monte_carlo=true, mesh=dosmesh,
        num_blocks=dos_num_blocks, histogram_width=histogram_width) 
        DOS[argmin(abs.(Es .- μ))] #Density of states at fermi energy
    end
    println("DOS at Fermi Energy is: ", gμ)

    relevantks, subsamplingfraction = isnothing(supplysampling) ? returnfermikpoint(Hwannier, cell_map, μ, Val(D); 
    histogram_width = histogram_width, mesh=dosmesh) : supplysampling

    nrelevantks = size(relevantks)[2]

    println("Sampling from ", nrelevantks, "kpoints")

    for i in 1:num_blocks
        println("Block: $i"); flush(stdout)

        ks = relevantks[:, rand(1:nrelevantks, mesh)]
        kprimes = relevantks[:, rand(1:nrelevantks, mesh)]

        Eks, Uks = wannier_bands(Hwannier, cell_map, ks)
        Ekprimes, Ukprimes = wannier_bands(Hwannier, cell_map, kprimes)

        omegaphsquareds, Uphs = diagonalize_phonon(force_matrix, cellph_map, ks, kprimes)
        omegaphs = sqrt.(abs.(omegaphsquareds))

        ephmatrixelements = eph_matrix_elements(Heph, celleph_map, Uks, Ukprimes, omegaphs, Uphs, ks, kprimes)
        ephmatrixelements = (abs.(ephmatrixelements)).^2
        nmodes = size(ephmatrixelements)[3]
        numbands = size(Eks)[2]

        vks = np.einsum("kpaa -> kpa ", imag.(momentum_matrix_elements(Hwannier, cell_map, Pwannier, ks))) ##Find momentum matrix elements for k 
        vkprimes = np.einsum("kpaa -> kpa ", imag.(momentum_matrix_elements(Hwannier, cell_map, Pwannier, kprimes)))
        
        vksabs = np.repeat(np.reshape(sqrt.(np.einsum("ijk -> ik", vks .^2 )), (mesh, 1, -1)), 3, axis=1)
        vks = vks ./ vksabs    
        vkprimesabs = np.repeat(np.reshape(sqrt.(np.einsum("ijk -> ik", vkprimes .^2 )), (mesh, 1, -1)), 3, axis=1)
        vkprimes = vkprimes ./ vkprimesabs    

        velocity_term = 1 .- np.einsum("kin, qim -> kqnm", vks, vkprimes)
        velocity_term = np.repeat(np.reshape(velocity_term, (mesh, mesh, 1, numbands, numbands)), nmodes, axis=2)

        weights = ephmatrixelements .* velocity_term 

        if weight == Val(:histogram)
            weights = np.einsum("kqlnm, kn, qm -> kqlnm", histogram_width^2*weights, abs.((Eks .- μ)*histogram_width) .< 0.5, abs.((Ekprimes .- μ)*histogram_width) .< 0.5)
        elseif weight == Val(:gaussian)
            weights = np.einsum("kqlnm, kn, qm -> kqlnm", 1/(2*pi*esmearing^2)*weights, exp.(-1/2*(Eks .- μ).^2 ./esmearing^2), exp.(-1/2*(Ekprimes .- μ).^2 ./esmearing^2))
        elseif weight == Val(:lorentzian)
            weights = np.einsum("kqlnm, kn, qm -> kqlnm", ((1/π)^2)*weights, imag.(1 ./ (Eks .- μ .+esmearing*1im)), imag.(1 ./ (Ekprimes .-μ .+esmearing*1im)))
        end

        omegaphs = np.repeat(np.repeat(np.reshape(omegaphs, (mesh, mesh, nmodes, 1, 1)), numbands, axis=3), numbands, axis=4) ./ eV
        omegaphs = np.reshape(omegaphs, (1, mesh, mesh, nmodes, numbands, numbands))
        omegaphs = np.repeat(omegaphs, length(ωs), axis=0)

        omegas = np.reshape(ωs, (-1, 1, 1, 1, 1, 1))
        omegas = np.repeat(np.repeat(np.repeat(np.repeat(np.repeat(omegas, mesh, axis=1), mesh, axis=2), nmodes, axis=3), numbands, axis=4), numbands, axis=5)

        weights = np.einsum("wkqlnm, kqlnm -> wkqlnm", (btomega.(omegas .+ omegaphs, fracroom).* bose.(max.(omegaphs, 3e-5), fracroom)  .- btomega.(omegas .- omegaphs, fracroom).*  bose.(-max.(omegaphs, 3e-5), fracroom)) ./ btomega.(omegas, fracroom), weights ) 
        weights = isnothing(intraband) ? weights : np.einsum("wkqlnm, n, m -> wkqlnm", weights, [x == intraband ? 1 : 0 for x in 1:numbands], [x == intraband ? 1 : 0 for x in 1:numbands])
        weights *= 1/mesh^2
        weights = np.einsum("wkqlnm -> w", weights)
        
        tauinv += weights / num_blocks
    end
    return 1e15 ./ ((2π/(ħ*gμ))*tauinv*subsamplingfraction^2)
end


"""
Decay time as derived from P. Allen, Physical Review B 3, 305 (1971).
"""
function τ_allen(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2},
    Heph::Array{Float64, 5}, celleph_map::Array{<:Real, 2}, ωs::Vector{<:Real}, μ::Real, ::Val{D}=Val(2); 
    histogram_width::Integer=10, supplydos::Union{Nothing, Real}=nothing, mesh::Integer=10, num_blocks::Integer=10, dosmesh::Integer = 10, dos_num_blocks=3, fracroom::Real=1) where D

    tauinv = zeros(length(ωs))

    gμ = 
    if !isnothing(supplydos)
            supplydos 
    else
        Es, DOS = density_of_states(Hwannier, cell_map, Val(D); monte_carlo=true, mesh=dosmesh,
        num_blocks=dos_num_blocks, histogram_width=histogram_width) 
        DOS[argmin(abs.(Es .- μ))] #Density of states at fermi energy
    end
    println("DOS at Fermi Energy is: ", gμ)

    for i in 1:num_blocks
        println("Block: $i"); flush(stdout)

        ks = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        kprimes = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))

        Eks, Uks = wannier_bands(Hwannier, cell_map, ks)
        Ekprimes, Ukprimes = wannier_bands(Hwannier, cell_map, kprimes)

        omegaphsquareds, Uphs = diagonalize_phonon(force_matrix, cellph_map, ks, kprimes)
        omegaphs = sqrt.(abs.(omegaphsquareds))

        ephmatrixelements = eph_matrix_elements(Heph, celleph_map, Uks, Ukprimes, omegaphs, Uphs, ks, kprimes)
        ephmatrixelements = (abs.(ephmatrixelements)).^2
        nmodes = size(ephmatrixelements)[3]
        numbands = size(Eks)[2]

        vks = np.einsum("kpaa -> kpa ", imag.(momentum_matrix_elements(Hwannier, cell_map, Pwannier, ks))) ##Find momentum matrix elements for k 
        vkprimes = np.einsum("kpaa -> kpa ", imag.(momentum_matrix_elements(Hwannier, cell_map, Pwannier, kprimes)))
        
        vksabs = np.repeat(np.reshape(sqrt.(np.einsum("ijk -> ik", vks .^2 )), (mesh^D, 1, -1)), 3, axis=1)
        vks = vks ./ vksabs    
        vkprimesabs = np.repeat(np.reshape(sqrt.(np.einsum("ijk -> ik", vkprimes .^2 )), (mesh^D, 1, -1)), 3, axis=1)
        vkprimes = vkprimes ./ vkprimesabs    

        velocity_term = 1 .- np.einsum("kin, qim -> kqnm", vks, vkprimes)
        velocity_term = np.repeat(np.reshape(velocity_term, (mesh^D, mesh^D, 1, numbands, numbands)), nmodes, axis=2)

        weights = ephmatrixelements .* velocity_term 

        omegaphs = np.repeat(np.repeat(np.reshape(omegaphs, (mesh^D, mesh^D, nmodes, 1, 1)), numbands, axis=3), numbands, axis=4) ./ eV
        omegaphs = np.reshape(omegaphs, (1, mesh^D, mesh^D, nmodes, numbands, numbands))
        omegaphs = np.repeat(omegaphs, length(ωs), axis=0)

        omegas = np.reshape(ωs, (-1, 1, 1, 1, 1, 1))
        omegas = np.repeat(np.repeat(np.repeat(np.repeat(np.repeat(omegas, mesh^D, axis=1), mesh^D, axis=2), nmodes, axis=3), numbands, axis=4), numbands, axis=5)
        
        Eks = np.reshape(Eks, (1, mesh^D, 1, 1, numbands, 1))
        Eks = np.repeat(Eks, length(ωs), axis=0)
        Eks = np.repeat(Eks, mesh^D, axis=2)
        Eks = np.repeat(Eks, nmodes, axis=3)
        Eks = np.repeat(Eks, numbands, axis=5)
        
        Ekprimes = np.reshape(Ekprimes, (1, 1, mesh^D, 1, 1, numbands))
        Ekprimes = np.repeat(Ekprimes, length(ωs), axis=0)
        Ekprimes = np.repeat(Ekprimes, mesh^D, axis=1)
        Ekprimes = np.repeat(Ekprimes, nmodes, axis=3)
        Ekprimes = np.repeat(Ekprimes, numbands, axis=4)

        energy_conservation_emission = abs.(Eks - Ekprimes + omegas - omegaphs)*histogram_width .< 0.5
        energy_conservation_emission *= histogram_width

        energy_conservation_absorption = abs.(Eks - Ekprimes + omegas + omegaphs)*histogram_width .< 0.5
        energy_conservation_absorption *= histogram_width

        fks = fermi.(Eks .- μ , fracroom) 
        fkprimes = fermi.(Ekprimes .- μ, fracroom)
        Nqs = bose.(omegaphs, fracroom)
        
        energy_conservation_total = energy_conservation_absorption .* (Nqs .* fks - (1 .+ Nqs).* fkprimes .+ (fks .* fkprimes)) .+
                                    energy_conservation_emission .* ((Nqs .+ 1) .* fks - Nqs.* fkprimes .- (fks .* fkprimes)) 
                                    # Note that if phonon distribution is 0 (low temperature) only phonon emission contributes and we get fk(1-fkprime)
        weights = np.einsum("wkqlnm, kqlnm -> w", energy_conservation_total, weights) / mesh^(2*D)
        weights .*= 1 ./ ωs
        tauinv += weights / num_blocks
    end
    return 1e15 ./ ((2*π/(ħ*gμ))*tauinv) # Returns result in femtoseconds
end

"""
$(TYPEDSIGNATURES)
Interband conductivity as defined in the paper Plasmonics in argentene by Shankar. Instead of using the bloch eigenvectors, this
calculates the local conductivity by interpolating the momentum matrix elements. 

"""
function interbandsigma(lattice_vectors::Vector{<:Vector{Float64}}, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    Pwannier::Array{Float64, 4}, μ::Real, ::Val{D}=Val(2); mesh::Integer =10, num_blocks::Integer = 10,
    histogram_width::Integer=10, energy_range::Real=20, degeneracy::Integer=1) where D

    V = D == 2 ? unit_cell_area(lattice_vectors) : unit_cell_volume(lattice_vectors)
    prefactor = degeneracy*4*π*ħ^2/V
    conds = zeros(ComplexF64, (3, 3, round(Int, histogram_width*energy_range)))
    for i in 1:num_blocks
        println("Block: $i"); flush(stdout)
        ks = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        Eks, Uks= wannier_bands(Hwannier, cell_map, ks);
        numbands = size(Eks)[2]
        Eksn = np.repeat(np.reshape(Eks, (-1, numbands, 1)), numbands, axis=2)
        Eksm = np.repeat(np.reshape(Eks, (-1, 1, numbands)), numbands, axis=1)
        vnks = momentum_matrix_elements(Uks, cell_map, Pwannier, ks)/mₑ
        omeganm = Eksn - Eksm
        weights = prefactor*histogram_width/mesh^D*np.einsum("knm, kpnm, krnm, knm -> kprnm", -np.heaviside(μ .- Eksn, 0.5)  + np.heaviside(μ .- Eksm, 0.5),
        conj.(vnks), vnks, 1 ./ omeganm )
        for (a, b) in Tuple.(CartesianIndices(rand(3, 3)))
            y, _ = np.histogram(omeganm, weights=weights[:, a, b, :, :], bins=round(Int, energy_range*histogram_width), range=(0, energy_range))
            conds[a, b, :] += y
        end
    end
    return conds/num_blocks
end

function Σ(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2},
    Heph::Array{Float64, 5}, celleph_map::Array{<:Real, 2}, ϵ_range::Vector{<:Real}, μ::Real, ::Val{D}=Val(2); 
    histogram_width::Integer=10, mesh::Integer=10, num_blocks::Integer=10, fracroom::Real=1) where D
    
    num_ϵ = length(ϵ_range)
    Σ_kk_ϵ = zeros(num_ϵ)
    nums_at_epsilon = zeros(num_ϵ)
    diff_ϵ = 1/first(diff(ϵ_range))

    for i in 1:num_blocks
        println("Block: $i"); flush(stdout)
        ks = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        kprimes = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
    
        Eks, Uks = wannier_bands(Hwannier, cell_map, ks)
        Ekprimes, Ukprimes = wannier_bands(Hwannier, cell_map, kprimes);
        nbands = size(Eks)[2]
        Eks_ϵ_ϵk = np.repeat(np.reshape(Eks, (1, mesh^D, nbands)), num_ϵ, axis=0)
        ϵs_ϵ_ϵk =  np.repeat(np.repeat(np.reshape(ϵ_range, (-1, 1, 1)), mesh^D, axis=1), nbands, axis=2)
        
        δ_argument = abs.(Eks_ϵ_ϵk .- ϵs_ϵ_ϵk)
        δ_ϵ_ϵk = Int.(δ_argument*diff_ϵ .<  0.5)
        omegaphsquareds, Uphs = diagonalize_phonon(force_matrix, cellph_map, ks, kprimes)
        omegaphs = sqrt.(abs.(omegaphsquareds))
        
        ephmatrixelements = eph_matrix_elements(Heph, celleph_map, Uks, Ukprimes, omegaphs, Uphs, ks, kprimes)
        ephmatrixelements_sqrd  = (abs.(ephmatrixelements)).^2
        omegaphs = omegaphs ./ eV;
        nmodes = size(ephmatrixelements)[3]
    
        Eks = np.repeat(np.repeat(np.reshape(Eks, (mesh^D, 1, 1, nbands, 1)), mesh^D, axis=1), nmodes, axis=2)
        Eks = np.repeat(Eks, nbands, axis=4)
        
        Ekprimes = np.repeat(np.repeat(np.reshape(Ekprimes, (1, mesh^D, 1, 1, nbands)), mesh^D, axis=0), nmodes, axis=2)
        Ekprimes = np.repeat(Ekprimes, nbands, axis=3)
        
        omegaphs = np.repeat(np.repeat(np.reshape(omegaphs, (mesh^D, mesh^D, nmodes, 1, 1)), nbands, axis=4), nbands, axis=3)
        delta_emission = abs.(Eks - Ekprimes - omegaphs)*histogram_width .< 0.5
        delta_absorption = abs.(Eks - Ekprimes + omegaphs)*histogram_width .< 0.5
    
        Nqs = bose.(max.(omegaphs, 3e-5), fracroom)
        fkprimes = fermi.(Ekprimes .- μ, fracroom)
        
        emission_term =  (Nqs .+ 1 - fkprimes) .* delta_emission
        absorption_term = (Nqs + fkprimes) .* delta_absorption
        
        y = np.einsum("kqanm, kqanm, wkn -> w", emission_term + absorption_term, ephmatrixelements_sqrd, δ_ϵ_ϵk)
        nums_at_epsilon += np.einsum("wkn -> w", δ_ϵ_ϵk) / num_blocks
        Σ_kk_ϵ += π*y*histogram_width/(num_blocks*mesh^D)
    end
    return nums_at_epsilon, Σ_kk_ϵ 
end
