
#The methods here are designed after "Plasmonics in Argentene by Shankar"

#Below, we try to calculate 1/τ(ω) as defined in cited paper above. 


"""
$(TYPEDSIGNATURES)

Returns the Fermi surface averaged velocity squared. Note, the current assumption is that the velocity is isotropic. 

This corresponds to the coefficient of 1/ω (ω in eV units) in the Drude formula in units of e²/(4ħ)

"""
function drude_conductivity(lattice_vectors::Vector{<:Vector{Float64}}, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    Pwannier::Array{Float64, 4}, μ::Real, ::Val{D}=Val(2); mesh::Integer =10, num_blocks::Integer=10, histogram_width::Integer=10, degeneracy::Integer=1) where D
    V = D == 2 ? unit_cell_area(lattice_vectors) : unit_cell_volume(lattice_vectors) 
    σ =  zeros(3, 3)
    gs = degeneracy
    for i in 1:num_blocks
        println("Block: $i"); flush(stdout)
        ks = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        Eks, Uks= wannier_bands(Hwannier, cell_map, ks);
        vks = np.einsum("kpaa -> kpa ", imag.(momentum_matrix_elements(Uks, cell_map, Pwannier, ks)))  
        σ += np.einsum("kn, knml -> ml", abs.(Eks .- μ)*histogram_width .< 0.5, np.einsum("kmn, kln->knml", vks/mₑ, vks/mₑ))*1/mesh^D*gs*(4*ħ^2)*histogram_width*(1/V)
    end
    return σ/num_blocks
end

function btomega(ω::Real, fracroom::Real)
#fracroom is fraction of room finite_temperature_chemical_potential
    room = 11606/298
    return ω/(1-exp(-ω*room/fracroom)) #Note that 1/40 eV is room temperature so β = 40 at room temperature
end

function τ(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2},
    Heph::Array{Float64, 5}, celleph_map::Array{<:Real, 2}, ωs::Vector{<:Real}, μ::Real, ::Val{D}=Val(2); histogram_width::Integer=10,
    supplysampling::Union{Tuple{<:Matrix{<:Real}, Real}, Nothing}=nothing, supplydos::Union{Nothing, Real}=nothing, mesh::Integer=10,
    num_blocks::Integer=10, dosmesh::Integer = 10, dos_num_blocks=3, fracroom::Real=1) where D
    
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
        weights = np.einsum("kqlnm, kn, qm -> kqlnm", histogram_width^2*weights, abs.((Eks .- μ)*histogram_width) .< 0.5, abs.((Ekprimes .- μ)*histogram_width) .< 0.5)

        omegaphs = np.repeat(np.repeat(np.reshape(omegaphs, (mesh, mesh, nmodes, 1, 1)), numbands, axis=3), numbands, axis=4) ./ eV
        omegaphs = np.reshape(omegaphs, (1, mesh, mesh, nmodes, numbands, numbands))
        omegaphs = np.repeat(omegaphs, length(ωs), axis=0)

        omegas = np.reshape(ωs, (-1, 1, 1, 1, 1, 1))
        omegas = np.repeat(np.repeat(np.repeat(np.repeat(np.repeat(omegas, mesh, axis=1), mesh, axis=2), nmodes, axis=3), numbands, axis=4), numbands, axis=5)

        temp_weight = (btomega.(omegas .- omegaphs, fracroom) ./ btomega.(omegas, fracroom))
        weights = np.einsum("wkqlnm, kqlnm -> w", temp_weight, weights) / mesh^2
        tauinv += weights / num_blocks
    end

    return 1e15 ./ ((4π^2/(ħ*gμ))*tauinv*subsamplingfraction^2)
end

"""
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
