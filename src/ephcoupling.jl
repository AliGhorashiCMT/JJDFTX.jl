"""
$(TYPEDSIGNATURES)
Returns the electron phonon coupling, h(ϵ) as defined in Ab initio phonon coupling and optical response of hot electrons in plasmonic metals
"""
function ephcoupling(μ::Real, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, force_matrix::Array{Float64, 3}, 
    cellph_map::Array{Float64, 2}, Heph::Array{Float64, 5}, celleph_map::Array{<:Real, 2}, ::Val{D} = Val(2), 
    weight::Union{Val{:gaussian}, Val{:lorentzian}, Val{:histogram}} = Val(:histogram); dosmesh::Integer=1, mesh::Integer=1,
    num_blocks::Integer = 10, dos_num_blocks::Integer =10, histogram_width::Integer=10, esmearing::Real=1, verbose::Bool=true) where D
    
    Es, DOS = density_of_states(Hwannier, cell_map, Val(D); monte_carlo=true, mesh=dosmesh,
    num_blocks=dos_num_blocks, histogram_width=histogram_width) 

    hϵ = zeros(length(Es))
    gμ = DOS[argmin(abs.(Es .- μ))] #Density of states at fermi energy

    verbose && println("Density of States at Fermi Energy: ", gμ)
    for i in 1:num_blocks
        verbose && println("Block: $i"); flush(stdout)
        Energies = Es
        ks = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        kprimes = vcat(rand(D, mesh^D), zeros(3-D, mesh^D))

        Eks, Uks = wannier_bands(Hwannier, cell_map, ks) 
        Ekprimes, Ukprimes = wannier_bands(Hwannier, cell_map, kprimes)

        numbands = size(Eks)[2]

        omegaphsquareds, Uphs = diagonalize_phonon(force_matrix, cellph_map, ks, kprimes)
        omegaphs = sqrt.(abs.(omegaphsquareds))

        ephmatrixelements = eph_matrix_elements(Heph, celleph_map, Uks, Ukprimes, omegaphs, Uphs, ks, kprimes)
        ephmatrixelements = (abs.(ephmatrixelements)).^2
        nmodes = size(ephmatrixelements)[3]

        Energies = np.reshape(Energies, (-1, 1, 1))

        Energies = np.repeat(Energies, mesh^D, axis=1)
        Energies = np.repeat(Energies, numbands, axis=2)

        Eks = np.reshape(Eks, (1, -1, numbands))
        Eks = np.repeat(Eks, size(Energies)[1], axis=0)

        Ekprimes = np.reshape(Ekprimes, (1, -1, numbands))
        Ekprimes = np.repeat(Ekprimes, size(Energies)[1], axis=0)

        dks, dkprimes = 
        if weight == Val(:histogram)
            histogram_width*(abs.(Energies .- Eks)*histogram_width .< 0.5),
            histogram_width*(abs.(Energies .- Ekprimes)*histogram_width .< 0.5)
        elseif weight == Val(:gaussian)
            1/sqrt(2π*esmearing^2)*exp.(-1/2*(Energies .- Eks).^2 ./esmearing^2), 
            1/sqrt(2π*esmearing^2)*exp.(-1/2*(Energies .- Ekprimes).^2 ./esmearing^2)
        elseif weight == Val(:lorentzian)
            1/π*imag.(1 ./ (Eks .- Energies .+esmearing*1im)),
            1/π*imag.(1 ./ (Ekprimes .- Energies .+esmearing*1im))
        end
        weights = np.einsum("e, kqw, kqwnm, ekn, eqm -> ekqwnm",  1 ./ (DOS).^2, omegaphs ./ eV, ephmatrixelements, 
        dks, dkprimes)

        weights *= 2*gμ/(mesh^(2*D))

        Energies = np.reshape(Energies, (-1, mesh^D, 1, 1, numbands, 1))
        Energies = np.repeat(np.repeat(np.repeat(Energies, mesh^D, axis=2), nmodes, axis=3), numbands, axis=5)

        hϵ += first(np.histogram(Energies, bins = length(hϵ), 
            range=(minimum(Es), maximum(Es)), weights = weights))
    end
    return Es, hϵ ./ num_blocks
end