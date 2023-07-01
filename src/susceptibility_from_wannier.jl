"""
$(TYPEDSIGNATURES)

"""
function ImΠ(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real},
    μ::Real, dim::Val{D}=Val(2), wannier_centers::Array{<:Real, 2}=zeros(size(Hwannier)[2], 3); 
    degeneracy::Integer=1, mesh::Integer=100, num_blocks::Integer=10, histogram_width::Integer=100, monte_carlo::Bool = false, verbose::Bool=true, normalized::Bool=true) where D

    verbose && println(q)
    Polarization_Array=zeros(histogram_width*100)
    V = 
        if D == 2 
            unit_cell_area(lattice_vectors) 
        elseif D ==3
            unit_cell_volume(lattice_vectors)
        end

    qnormalized = normalized ? q : normalize_kvector(lattice_vectors, q)

    numbands = size(Hwannier)[2]

    for _ in 1:num_blocks
        kpoints = !monte_carlo ? transpose(make_mesh(mesh, dim)) : vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        kphases = np.exp(-2π*1im*wannier_centers * kpoints)
        
        kplusqpoints = reshape(repeat(qnormalized, mesh^D), (3, mesh^D)) + kpoints
        kplusqphases = np.exp(-2π*1im*wannier_centers * kplusqpoints)

        Eks, Uks = wannier_bands(Hwannier, cell_map, kpoints)
        Ekqs, Ukqs = wannier_bands(Hwannier, cell_map, kplusqpoints)

        Uks = np.einsum("il, lij-> lij ", kphases, Uks)
        Ukqs = np.einsum("il, lij-> lij ", kplusqphases, Ukqs)
        
        overlaps = np.einsum("lji, ljk -> lik", np.conj(Uks), Ukqs) # l indexes the k point, i and k index the band indices
        overlaps = overlaps .* np.conj(overlaps) # lij component is |<i, k_l| j, k_l+q>|^2

        Ekqs_reshaped = np.repeat(np.reshape(Ekqs, (mesh^D, 1, numbands)), numbands, axis=1)
        Eks_reshaped = np.repeat(np.reshape(Eks, (mesh^D, numbands, 1)), numbands, axis=2)
        omegas = Ekqs_reshaped - Eks_reshaped # lij component is E(k_l + q)_j - E(k_l)_i

        f2 = np.heaviside(μ .- Ekqs_reshaped, 0.5)
        f1 = np.heaviside(μ .- Eks_reshaped, 0.5)

        summand = (f2 - f1) .* overlaps

        Polarization_Array += first(np.histogram(omegas, bins=round(Int, histogram_width*100), weights=summand, range=(0, 100)))*π/V*(1/mesh)^D*histogram_width*degeneracy*(1/num_blocks)
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)
Note that this gives the 2d conductivity in units of the universal conductivity, e²/4ħ

"""
function σ(q::Vector{<:Real}, lattice_vectors::Vector{<:Vector{<:Real}}, ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real}, 
    ; normalized::Bool=true, δ::Real=0.01) 
    qabs = normalized ? norm(unnormalize_kvector(lattice_vectors, q)) : norm(q)
    return 4im*ω/(qabs^2)*(real(kramers_kronig(ω, energies, polarizations; δ)) + 1im*polarizations[argmin(abs.(energies .- ω))])
end

"""
$(TYPEDSIGNATURES)
Returns the non-local, non-static dielectric function for dimensions 2 and 3
"""
function ϵ(q::Vector{<:Real}, lattice_vectors::Vector{<:Vector{<:Real}}, ω::Real, energies::Vector{<:Real}, imaginary_polarizations::Vector{<:Real}, 
    ::Val{D}=Val(2), algorithm::Union{Val{:default}, Val{:scipy}, Val{:quadgk}}=Val(:default); normalized::Bool=true, δ::Real=0.01, win_len::Integer=10, kwargs...) where D
    qabs = normalized ? norm(unnormalize_kvector(lattice_vectors, q)) : norm(q)

    real_polarizations = 
        if algorithm == Val(:default)
            kramers_kronig(ω, energies, imaginary_polarizations; δ)
        elseif algorithm == Val(:scipy)
            kramers_kronig_scipy(ω, energies, imaginary_polarizations; kwargs...)
        elseif algorithm == Val(:quadgk)
            first(kramers_kronig_quadgk(ω, energies, imaginary_polarizations; δ, kwargs...))
        end
    return 1-e²ϵ/((4-D)*qabs^(D-1))*(real_polarizations + 1im*smooth(imaginary_polarizations, win_len=win_len)[argmin(abs.(energies .- ω))])
end

"Returns the plasmon confinement factor "
function confinement(lattice_vectors::Vector{<:Vector{<:Real}}, qs::Vector{<:Vector{<:Real}}, omegas::Vector{<:Real}; normalized::Bool = true) 
    qabsv = normalized ? [norm(q) for q in unnormalize_kvector.(Ref(lattice_vectors), qs )] : [norm(q) for q in qs]
    confinements = c*ħ*qabsv ./ omegas
    return qabsv, confinements
end