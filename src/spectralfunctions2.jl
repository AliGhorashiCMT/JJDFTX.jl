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
function eliashberg(lattice_vectors::Vector{<:Vector{<:Real}}, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, 
    force_matrix::Array{Float64, 3}, cellph_map::Array{Float64, 2}, Heph::Array{Float64, 5}, celleph_map::Array{<:Real, 2}, μ::Real, ::Val{D} = Val(2), 
    weight::Union{Val{:gaussian}, Val{:lorentzian}, Val{:histogram}} = Val(:histogram); mesh::Integer=10, num_blocks::Integer=10, esmearing::Real=10, histogram_width::Real=10, histogram_width2::Real=3, energyrange::Real=1,
    subsampling::Union{Tuple{<:Matrix{<:Real}, Real}, Nothing}=nothing, dosmu::Union{Real, Nothing}=nothing) where D

    relevantks, subsamplingfraction = isnothing(subsampling) ? returnfermikpoint(Hwannier, cell_map, μ, Val(D); 
    histogram_width = histogram_width2, mesh=10) : subsampling

    nrelevantks = size(relevantks)[2]
    V = unit_cell_volume(lattice_vectors)
    α²F = zeros(Int, energyrange*histogram_width2)
        
    for _ in 1:num_blocks
        ks = relevantks[:, rand(1:nrelevantks, mesh)] # Monte Carlo sampling. Choose an index of a k point at the Fermi level 
        kprimes = relevantks[:, rand(1:nrelevantks, mesh)] # Monte Carlo sampling
        
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
            weights = np.einsum("kqlnm, kn, qm -> kqlnm", (1/π)^2*weights, imag.(1 ./ (Eks .- μ .+esmearing*1im)), imag.(1 ./ (Ekprimes .-μ .+esmearing*1im)))
        end

        weights = V^2*(histogram_width2*weights) ./(dosmu)^2 / (mesh^2)

        omegaphs = np.repeat(np.repeat(np.reshape(omegaphs, (mesh, mesh, nmodes, 1, 1)), numbands, axis=3), numbands, axis=4) ./ eV
        println(maximum(omegaphs))
        y, _ = np.histogram(omegaphs, bins = round(Int, energyrange*histogram_width2), range=(0, energyrange), weights = weights)
        α²F += y 
    end

    return (α²F / num_blocks)*subsamplingfraction*subsamplingfraction #Because we only looked at Fermi kvectors and not arbitrary kvectors   
end
