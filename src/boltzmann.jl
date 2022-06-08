
#The methods here are designed after "Plasmonics in Argentene by Shankar"

#Below, we try to calculate 1/τ(ω) as defined in cited paper above. 


"""
$(TYPEDSIGNATURES)

Returns the Fermi surface averaged velocity squared. Note, the current assumption is that the velocity is isotropic. 

This corresponds to the coefficient of 1/ω (ω in eV units) in the Drude formula in units of e²/(4ħ)

"""
function drude_conductivity(lattice::Vector{<:Vector{Float64}}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, 
    PWannier::Array{Float64, 4}, nbands::Integer, μ::Real; mesh =10, histogram_width=10)
    area = unit_cell_area(lattice);
    σ =  0
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

    Es, DOS = density_of_states(Hwannier, cell_map, Val(D); monte_carlo=true, mesh=dosmesh,
    num_blocks=dos_num_blocks, histogram_width=histogram_width) 
    gμ = DOS[argmin(abs.(Es .- μ))] #Density of states at fermi energy

    println("DOS at Fermi Energy is: ", gμ)

    relevantks, subsamplingfraction = isnothing(supplysampling) ? returnfermikpoint(Hwannier, cell_map, μ, Val(D); 
    histogram_width = histogram_width, mesh=10) : supplysampling

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

    return (4π^2/(ħ*gμ))*tauinv*subsamplingfraction^2
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
