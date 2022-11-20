"""
$(TYPEDSIGNATURES)

"""
function ImΠ(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real, dim::Val{D}=Val(2); 
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

        kplusqpoints = reshape(repeat(qnormalized, mesh^D), (3, mesh^D)) + kpoints

        Eks, Uks = wannier_bands(Hwannier, cell_map, kpoints)
        Ekqs, Ukqs = wannier_bands(Hwannier, cell_map, kplusqpoints)
        
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

Given a one band Wannier tight binding hamiltonian, this returns the integrand corresponding to the RPA polarization function. 
Note that the extra factor of 2 is due to the fact that time reversal symmetry has been used to condense the integrand into a function 
of just the occupation and not a difference of occupations. This is important since typically numerical algorithms will give incorrect answers
for integrands that are 0 in a large region of phase space.

By default, spin degeneracy is not taken into account, but this can be changed by altering the value of the keyword argument spin. 
"""
function epsilon_integrand(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k₁::Real, k₂::Real, q::Vector{<:Real},
    μ::Real, ω::Real, δ::Real; degeneracy::Integer=1)

    kvector=[k₁, k₂, 0]
    ϵ₁ = first(first(wannier_bands(Hwannier, cell_map, kvector)))
    ϵ₂ = first(first(wannier_bands(Hwannier, cell_map, kvector+q)))
    f = ϵ₁ < μ ? 1 : 0
    integrand = 1/(2π)^2*degeneracy*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*δ)^2)
    return integrand
end

"""
$(TYPEDSIGNATURES)
For calculating ϵ(q, ω) without doing Kramers-Kronig. Due to numerical algorithm limitations, this should only be used 
for intraband (one defect band) calculations.
"""
function ϵ(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}},
    q::Vector{<:Real}, ω::Real, μ::Real, algorithm::Union{Val{:default}, Val{:cubature}}; degeneracy::Integer = 1, δ::Real = 0.01, normalized::Bool=true, kwargs...) 
    kwargsdict=Dict()
    for kwarg in kwargs
        push!(kwargsdict, kwarg.first => kwarg.second)
    end
    qnormalized =  normalized ? q : normalize_kvector(lattice_vectors, q) 
    qabs = normalized ? norm(unnormalize_kvector(lattice_vectors, q)) : norm(q)
    brillouin_area = brillouin_zone_area(lattice_vectors) 
    polarization =  
        if algorithm == Val(:default)
            brillouin_area*pyintegrate.nquad((k₁, k₂) -> real(epsilon_integrand(Hwannier, cell_map, k₁, k₂, 
            qnormalized, μ, ω, δ, degeneracy=degeneracy)), [[0, 1], [0, 1]], opts=kwargsdict)[1]
        elseif algorithm == Val(:cubature)
            brillouin_area*hcubature((k) -> real(epsilon_integrand(Hwannier, cell_map, k[1], k[2], 
            qnormalized, μ, ω, δ, degeneracy=degeneracy)), [0, 0], [1, 1]; kwargs...)[1]
        end
    return 1-e²ϵ/(2qabs)*polarization
end

function im_polarization_cubature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, ω::Real, μ::Real; spin::Integer = 1, ϵ::Real=0.01, kwargs...) 
    qnormalized = normalize_kvector(lattice_vectors, q)
    brillouin_area=brillouin_zone_area(lattice_vectors) 
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand_imaginary(HWannier, cell_map, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]
    return polarization
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
    ::Val{D}=Val(2), algorithm::Union{Val{:default}, Val{:scipy}, Val{:quadgk}}=Val(:default); normalized::Bool=true, δ::Real=0.01, kwargs...) where D
    qabs = normalized ? norm(unnormalize_kvector(lattice_vectors, q)) : norm(q)

    real_polarizations = 
        if algorithm == Val(:default)
            kramers_kronig(ω, energies, imaginary_polarizations; δ)
        elseif algorithm == Val(:scipy)
            kramers_kronig_scipy(ω, energies, imaginary_polarizations; kwargs...)
        elseif algorithm == Val(:quadgk)
            first(kramers_kronig_quadgk(ω, energies, imaginary_polarizations; δ, kwargs...))
        end
    return 1-e²ϵ/abs((4-D)*qabs)*(real_polarizations + 1im*imaginary_polarizations[argmin(abs.(energies .- ω))])
end

"Returns the plasmon confinement factor "
function confinement(lattice_vectors::Vector{<:Vector{<:Real}}, qs::Vector{<:Vector{<:Real}}, omegas::Vector{<:Real}; normalized::Bool = true) 
    qabsv = normalized ? [norm(q) for q in unnormalize_kvector.(Ref(lattice_vectors), qs )] : [norm(q) for q in qs]
    confinements = c*ħ*qabsv ./ omegas
    return qabsv, confinements
end