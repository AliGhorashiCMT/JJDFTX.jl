"""
$(TYPEDSIGNATURES)
Returns the imaginary value of the polarization at frequency omega (eV) and wavevector q (inverse angstrom).
Several methods are provided. Wannier and cell map data may be given either through file names or through passing in 
HWannier and cell-map as dim 3 and dim 2 arrays of floats, respectively.
"""
function im_polarization(wannier_file::AbstractString, cell_map_file::AbstractString, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real;
    spin::Integer = 1, mesh::Integer = 100, histogram_width::Real = 100, normalized::Bool=false) 

    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalized ? q : normalize_kvector(lattice_vectors, q)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/mesh, j/mesh, 0]
        E1=wannier_bands(wannier_file, cell_map_file, kvector  )
        E2=wannier_bands(wannier_file, cell_map_file, kvector+qnormalized  )
        f1=np.heaviside( μ-E1, 0.5)
        f2=np.heaviside( μ-E2, 0.5)
        DeltaE=E2-E1
        DeltaE < 0 && continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)
"""
function im_polarization(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real;
    spin::Integer=1, mesh::Integer=100, histogram_width::Real=100, normalized::Bool=false) 
    
    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalized ? q : normalize_kvector(lattice_vectors, q)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/mesh, j/mesh, 0]
        E1 = wannier_bands(HWannier, cell_map, kvector  )
        E2 = wannier_bands(HWannier, cell_map, kvector+qnormalized  )
        f1 = np.heaviside( μ-E1, 0.5)
        f2 = np.heaviside( μ-E2, 0.5)
        DeltaE = E2-E1
        DeltaE < 0 && continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)

"""
function im_polarization_mc(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; spin::Int=1, mesh::Int=100, histogram_width::Real=100, normalized::Bool=false) 
    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalized ? q : normalize_kvector(lattice_vectors, q)
    xmesh = rand(mesh)
    ymesh = rand(mesh)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i, j, 0]
        E1 = wannier_bands(HWannier, cell_map, kvector  )
        E2 = wannier_bands(HWannier, cell_map, kvector+qnormalized  )
        f1=np.heaviside( μ-E1, 0.5)
        f2=np.heaviside( μ-E2, 0.5)
        DeltaE=E2-E1
        DeltaE < 0 && continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)
"""
function im_polarization(wannier_file::AbstractString, cell_map_file::AbstractString, lattvectors::lattice, q::Array{<:Real, 1}, μ::Real; spin::Integer = 1, mesh::Integer = 100, histogram_width::Real = 100)
    Polarization_Array=zeros(histogram_width*100)
    lattice_vectors = [lattvectors.rvectors[:, 1]*bohrtoangstrom, lattvectors.rvectors[:, 2]*bohrtoangstrom, lattvectors.rvectors[:, 3]*bohrtoangstrom]
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/mesh, j/mesh, 0]
        E1=wannier_bands(wannier_file, cell_map_file, kvector)
        E2=wannier_bands(wannier_file, cell_map_file, kvector+qnormalized)
        f1=np.heaviside( μ-E1, 0.5)
        f2=np.heaviside( μ-E2, 0.5)
        DeltaE=E2-E1
        DeltaE>0 || continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)
"""
function im_polarization(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattvectors::lattice, 
    q::Vector{<:Real}, μ::Real; spin::Integer=1, mesh::Integer=100, histogram_width::Integer=100, normalized::Bool=false) 
    Polarization_Array=zeros(histogram_width*100)
    lattice_vectors = [lattvectors.rvectors[:, 1]*bohrtoangstrom, lattvectors.rvectors[:, 2]*bohrtoangstrom, lattvectors.rvectors[:, 3]*bohrtoangstrom]
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalized ? q : normalize_kvector(lattice_vectors, q)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/mesh, j/mesh, 0]
        E1=wannier_bands(HWannier, cell_map, kvector)
        E2=wannier_bands(Hwannier, cell_map, kvector+qnormalized)
        f1=np.heaviside( μ-E1, 0.5)
        f2=np.heaviside( μ-E2, 0.5)
        DeltaE=E2-E1
        DeltaE>0 || continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)
Returns a Vector of impolarizations given a certain filling fraction
"""
function im_polarizationatfilling(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattvectors::Vector{<:Vector{<:Real}}, 
    filling::Real; spin::Integer=1, mesh::Integer=100, histogram_width::Integer=100, interpolate::Integer=1, 
    kpointsfile::AbstractString="bandstruct.kpoints", offset::Real=2, energy_range::Real=3)

    kpoints = bandstructkpoints2q(filename=kpointsfile, interpolate=interpolate)
    nks = length(kpoints)
    xenergies, yoccupations = find_chemical_potential(HWannier, cell_map; 
    mesh=100, histogram_width=100, energy_range=energy_range, offset=offset, plotoccupations=false)
    μ = xenergies[argmin(abs.(yoccupations .- filling))]
    println("μ is: ", μ)
    impols = Array{Float64, 2}(undef, (nks, histogram_width*100))
    for (idx, q) in enumerate(kpoints)
        println(q)
        impols[idx, :] = im_polarization(HWannier, cell_map, lattvectors, q, μ; spin=spin, mesh=mesh, histogram_width=histogram_width, normalized=true) 
    end
    return impols
end

"""
$(TYPEDSIGNATURES)
Returns a Vector of impolarizations given a certain filling fraction
"""
function im_polarizationatfilling(HWannierDefect::Array{Float64, 3}, HWannierUp::Array{Float64, 3}, HWannierDn::Array{Float64, 3}, 
    cellmapDefect::Array{Float64, 2},cellmapUp::Array{Float64, 2}, cellmapDn::Array{Float64, 2}, lattvectors::Vector{<:Vector{<:Real}}, 
    filling::Real; nbands::Integer=72, spin::Integer=1, mesh::Integer=100, histogram_width::Integer=100, interpolate::Integer=1, 
    kpointsfile::String="bandstruct.kpoints", offset::Real=2, energy_range::Real=3)

    kpoints = bandstructkpoints2q(filename=kpointsfile, interpolate=interpolate)
    nks = length(kpoints)
    xenergies, yoccupations = find_chemical_potential(HWannierDefect, cellmapDefect; 
    mesh=100, histogram_width=100, energy_range=energy_range, offset=offset, plotoccupations=false)
    μ = xenergies[argmin(abs.(yoccupations .- filling))]
    println("μ is: ", μ)
    impols = Array{Float64, 2}(undef, (nks, histogram_width*100))
    for (idx, q) in enumerate(kpoints)
        println(q)
        impols[idx, :] = im_polarization_mixedmesh(HWannierUp, HWannierDn, HWannierDefect, cellmapUp, cellmapDn, cellmapDefect, 
                                nbands, 36, 36, lattvectors, q, μ; spin=spin, intraband_mesh=mesh, interband_mesh=6, exclude_bands_up=[37], exclude_bands_dn=[], histogram_width=histogram_width, normalized=true) 
    end
    return impols
end


"""
$(TYPEDSIGNATURES)
"""
function im_polarization(wannier_file::AbstractString, cell_map_file::AbstractString, nbands::Integer, valence_bands::Integer, lattice_vectors::Array{<:Array{<:Real, 1},1}, 
    q::Array{<:Real, 1}, μ::Real; spin::Integer=1, mesh::Integer=100, histogram_width::Integer=100) 

    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/mesh, j/mesh, 0]
        E1=wannier_bands(wannier_file, cell_map_file, kvector, nbands)
        E2=wannier_bands(wannier_file, cell_map_file, kvector+qnormalized, nbands)
        V1=wannier_vectors(wannier_file, cell_map_file, kvector, nbands)
        V2=wannier_vectors(wannier_file, cell_map_file, kvector+qnormalized, nbands)
        for lower in 1:valence_bands+1
            for upper in valence_bands+1:nbands
                Elower = E1[lower]
                Eupper = E2[upper]
                overlap=(np.abs(np.dot(V1[:, lower], np.conj(V2[:, upper]))))^2;
                f1=np.heaviside( μ-Elower, 0.5)
                f2=np.heaviside( μ-Eupper, 0.5)
                DeltaE=Eupper-Elower
                DeltaE>0 || continue
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*overlap*(1/mesh)^2*histogram_width*spin
            end
        end
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)

### Keyword Arguments
"""
function im_polarization(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, valence_bands::Integer, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; 
    exclude_bands=Int[], spin::Integer=1, mesh::Integer=100, histogram_width::Integer=100, subset::Integer=1, Koffset::Vector{<:Real}=[0, 0, 0], verbose::Bool=true, normalized::Bool=false) 
    verbose && println(q)
    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalized ? q : normalize_kvector(lattice_vectors, q)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/(subset*mesh), j/(subset*mesh), 0] + Koffset
        E1=wannier_bands(HWannier, cell_map, kvector, nbands  )
        E2=wannier_bands(HWannier, cell_map, kvector+qnormalized, nbands  )
        V1=wannier_vectors(HWannier, cell_map, kvector)
        V2=wannier_vectors(HWannier, cell_map, kvector+qnormalized )
        for lower in 1:valence_bands+1
            for upper in valence_bands+1:nbands
                (lower ∉ exclude_bands || upper ∉ exclude_bands) || continue
                Elower = E1[lower]
                Eupper = E2[upper]
                overlap=(np.abs(np.dot(V1[:, lower], np.conj(V2[:, upper]))))^2;
                f1=np.heaviside( μ-Elower, 0.5)
                f2=np.heaviside( μ-Eupper, 0.5)
                DeltaE=Eupper-Elower
                DeltaE>0 || continue
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*overlap*(1/mesh)^2*histogram_width*spin*(1/subset^2)
            end
        end
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)
"""
function im_polarization_mc(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, valence_bands::Integer, 
    lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; exclude_bands::Array{Int, 1}=Int[], spin::Int=1, mesh::Int=100, histogram_width::Int=100) 
    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)
    xmesh = rand(mesh)
    ymesh = rand(mesh)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i, j, 0]
        E1=wannier_bands(HWannier, cell_map, kvector, nbands  )
        E2=wannier_bands(HWannier, cell_map, kvector+qnormalized, nbands  )
        V1=wannier_vectors(HWannier, cell_map, kvector)
        V2=wannier_vectors(HWannier, cell_map, kvector+qnormalized )
        for lower in 1:valence_bands+1
            for upper in valence_bands+1:nbands
                (lower ∉ exclude_bands || upper ∉ exclude_bands) || continue
                Elower = E1[lower]
                Eupper = E2[upper]
                overlap=(np.abs(np.dot(V1[:, lower], np.conj(V2[:, upper]))))^2;
                f1 = np.heaviside( μ-Elower, 0.5)
                f2 = np.heaviside( μ-Eupper, 0.5)
                DeltaE = Eupper-Elower
                DeltaE > 0 || continue
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*overlap*(1/mesh)^2*histogram_width*spin
            end
        end
    end
    return Polarization_Array
end

"""
For susceptibility calculations at finite temperature. Temperature is assumed to be provided in Kelvin. 
"""
function im_polarization_finite_temperature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real, T::Real; spin::Int=1, mesh::Int=100, histogram_width::Real=100) 
    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i/mesh, j/mesh, 0]
        E1 = wannier_bands(HWannier, cell_map, kvector  )
        E2 = wannier_bands(HWannier, cell_map, kvector+qnormalized  )
        f1=1/(1+exp((E1-μ)/(kB*T)))
        f2=1/(1+exp((E2-μ)/(kB*T)))
        DeltaE=E2-E1
        DeltaE>0 || continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end

"""
$(TYPEDSIGNATURES)
"""
function im_polarization_finite_temperature_mc(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real, T::Real; spin::Int=1, mesh::Int=100, histogram_width::Real=100) 
    Polarization_Array=zeros(histogram_width*100)
    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)
    xmesh = rand(mesh)
    ymesh = rand(mesh)
    for (i, j) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kvector=[i, j, 0]
        E1 = wannier_bands(HWannier, cell_map, kvector  )
        E2 = wannier_bands(HWannier, cell_map, kvector+qnormalized  )
        f1=1/(1+exp((E1-μ)/(kB*T)))
        f2=1/(1+exp((E2-μ)/(kB*T)))
        DeltaE=E2-E1
        DeltaE>0 || continue
        Polarization_Array[round(Int, histogram_width*DeltaE+1)] += π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
    end
    return Polarization_Array
end

function im_polarization(wannier_file_up::AbstractString, wannier_file_dn::AbstractString,  cell_map_file_up::AbstractString, 
    cell_map_file_dn::AbstractString, nbands::Integer, valence_bands_up::Integer, valence_bands_dn::Integer, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; kwargs...) 

    #Here we add the independent polarizations from different spin channels 
    spin_up_pol = im_polarization(wannier_file_up, cell_map_file_up, nbands, valence_bands_up, lattice_vectors, q, μ; kwargs... )
    spin_dn_pol = im_polarization(wannier_file_dn, cell_map_file_dn, nbands, valence_bands_dn, lattice_vectors, q, μ; kwargs... )
    return (spin_up_pol + spin_dn_pol)
end

function im_polarization(HWannierup::Array{Float64, 3}, HWannierdn::Array{Float64, 3},  cell_map_up::Array{Float64, 2}, cell_map_dn::Array{Float64, 2}, 
    nbands::Integer, valence_bands_up::Integer, valence_bands_dn::Integer, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; kwargs...)
    #Here we add the independent polarizations from different spin channels 
    spin_up_pol = im_polarization(HWannierup, cell_map_up, nbands, valence_bands_up, lattice_vectors, q, μ; kwargs... )
    spin_dn_pol = im_polarization(HWannierdn, cell_map_dn, nbands, valence_bands_dn, lattice_vectors, q, μ; kwargs... )
    return (spin_up_pol + spin_dn_pol)
end

"""
$(TYPEDSIGNATURES)

"""
function im_polarization_mixedmesh(HWannierup::Array{Float64, 3}, HWannierdn::Array{Float64, 3}, HWannierdefect::Array{Float64, 3},  cell_map_up::Array{Float64, 2}, cell_map_dn::Array{Float64, 2}, cell_map_defect::Array{Float64, 2}, 
    nbands::Integer, valence_bands_up::Integer, valence_bands_dn::Integer, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; interband_mesh::Int=10, intraband_mesh::Int=100, 
    win_len=50, exclude_bands_up = Int[], exclude_bands_dn = Int[], normalized::Bool=true, kwargs...)

    #Here we add the independent polarizations from different spin channels 
    spin_up_pol = im_polarization(HWannierup, cell_map_up, nbands, valence_bands_up, lattice_vectors, q, μ, mesh = interband_mesh, exclude_bands = exclude_bands_up, normalized=normalized; kwargs... )
    spin_dn_pol = im_polarization(HWannierdn, cell_map_dn, nbands, valence_bands_dn, lattice_vectors, q, μ, mesh = interband_mesh, exclude_bands = exclude_bands_dn, normalized=normalized; kwargs... )
    spin_defect_pol = im_polarization(HWannierdefect, cell_map_defect, lattice_vectors, q, μ, mesh = intraband_mesh, normalized=normalized; kwargs... )
    return (smooth(spin_up_pol + spin_dn_pol, win_len=win_len)+spin_defect_pol)
end

#=
For direct momentum matrix elements calculations 
=#
function im_polarization_mixedmesh(filebase::AbstractString, HWannierdefect::Array{Float64, 3},  cell_map_defect::Array{Float64, 2}, 
    nbands::Integer, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; intraband_mesh::Integer=100, win_len=50, kwargs...)

    #Here we add the independent polarizations from different spin channels 
    mixed_pol = nonwannierimpol(filebase, lattice_vectors, q, nbands, μ, Val(2), kwargs...)
    spin_defect_pol = im_polarization(HWannierdefect, cell_map_defect, lattice_vectors, q, μ, mesh = intraband_mesh; kwargs... )
    return (smooth(spin_up_pol + spin_dn_pol, win_len=win_len)+spin_defect_pol)
end

function im_polarization_mixedmesh_mc(HWannierup::Array{Float64, 3}, HWannierdn::Array{Float64, 3}, HWannierdefect::Array{Float64, 3},  cell_map_up::Array{Float64, 2}, cell_map_dn::Array{Float64, 2}, cell_map_defect::Array{Float64, 2}, nbands::Int, valence_bands_up::Int, valence_bands_dn::Int, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; interband_mesh::Int=10, intraband_mesh::Int=100, win_len=50, exclude_bands_up = Int[], exclude_bands_dn = Int[], kwargs...)
    #Here we add the independent polarizations from different spin channels 
    spin_up_pol = im_polarization_mc(HWannierup, cell_map_up, nbands, valence_bands_up, lattice_vectors, q, μ, mesh = interband_mesh, exclude_bands = exclude_bands_up; kwargs... )
    spin_dn_pol = im_polarization_mc(HWannierdn, cell_map_dn, nbands, valence_bands_dn, lattice_vectors, q, μ, mesh = interband_mesh, exclude_bands = exclude_bands_dn; kwargs... )
    spin_defect_pol = im_polarization_mc(HWannierdefect, cell_map_defect, lattice_vectors, q, μ, mesh = intraband_mesh; kwargs... )
    return (smooth(spin_up_pol + spin_dn_pol, win_len=win_len)+spin_defect_pol)
end

function epsilon_integrand(wannier_file::AbstractString, cell_map_file::AbstractString, k₁::Real, k₂::Real, q::Vector{<:Real}, 
    μ::Real, ω::Real, ϵ::Real; spin::Integer=1)

    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(wannier_file, cell_map_file, kvector,  )
    ϵ₂ =wannier_bands(wannier_file, cell_map_file, kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    real(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

function epsilon_integrand_imaginary(wannier_file::AbstractString, cell_map_file::AbstractString, k₁::Real, k₂::Real, 
    q::Vector{<:Real}, μ::Real, ω::Real, ϵ::Real; spin::Integer=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(wannier_file, cell_map_file, kvector  )
    ϵ₂ =wannier_bands(wannier_file, cell_map_file, kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    imag(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

"""
$(TYPEDSIGNATURES)
"""
function epsilon_integrand(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k₁::Real, k₂::Real, q::Vector{<:Real},
    μ::Real, ω::Real, ϵ::Real; spin::Integer=1)

    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(HWannier, cell_map, kvector  )
    ϵ₂ =wannier_bands(HWannier, cell_map, kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    real(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

"""
$(TYPEDSIGNATURES)
"""
function epsilon_integrand_imaginary(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k₁::Real, k₂::Real, q::Vector{<:Real},
    μ::Real, ω::Real, ϵ::Real; spin::Int=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(HWannier, cell_map,  kvector  )
    ϵ₂ =wannier_bands(HWannier, cell_map,  kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    imag(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

"""
$(TYPEDSIGNATURES)
"""
function direct_epsilon(wannier_file::AbstractString, cell_map_file::AbstractString, lattice_vectors::Vector{<:Vector{<:Real}},
     q::Vector{<:Real}, ω::Real, μ::Real; spin::Integer = 1, ϵ::Real = 0.01, normalized::Bool=true, kwargs...) 
    kwargsdict=Dict()
    for kwarg in kwargs
        push!(kwargsdict, kwarg.first => kwarg.second)
    end
    qnormalized =  normalized ? q : normalize_kvector(lattice_vectors, q) 
    qabs = normalized ? sqrt(sum(unnormalize_kvector(lattice_vectors, q).^2)) : sqrt(sum(lattice_vectors, q).^2)
    brillouin_area=brillouin_zone_area(lattice_vectors)   
    polarization=brillouin_area*pyintegrate.nquad((k₁, k₂) -> epsilon_integrand(wannier_file, cell_map_file, k₁, k₂, qnormalized, μ, ω, ϵ, spin=spin), [[0, 1], [0, 1]], opts=kwargsdict)[1]
    1-e²ϵ/(2qabs)*polarization
end

"""
$(TYPEDSIGNATURES)
For calculating ϵ(q, ω) without doing Kramers-Kronig. Due to numerical algorithm limitations, this should only be used 
for intraband (one defect band) calculations.
"""
function direct_epsilon(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}},
    q::Vector{<:Real}, ω::Real, μ::Real; spin::Integer = 1, ϵ::Real = 0.01, normalized::Bool=true, kwargs...) 
    kwargsdict=Dict()
    for kwarg in kwargs
        push!(kwargsdict, kwarg.first => kwarg.second)
    end
    qnormalized =  normalized ? q : normalize_kvector(lattice_vectors, q) 
    qabs = normalized ? sqrt(sum(unnormalize_kvector(lattice_vectors, q).^2)) : sqrt(sum(lattice_vectors, q).^2)
    brillouin_area=brillouin_zone_area(lattice_vectors) 
    polarization=brillouin_area*pyintegrate.nquad((k₁, k₂) -> epsilon_integrand(HWannier, cell_map, k₁, k₂, qnormalized, μ, ω, ϵ, spin=spin), [[0, 1], [0, 1]], opts=kwargsdict)[1]
    1-e²ϵ/(2qabs)*polarization
end

"""
$(TYPEDSIGNATURES)
"""
function direct_plasmon(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}},
    ωs::AbstractRange{<:Float64}, μ::Real; spin::Integer = 1, ϵ::Real = 0.01, normalized::Bool=true, interpolate::Integer=10, 
    kpointsfile::String="bandstruct.kpoints", kwargs...) 
    kpoints = bandstructkpoints2q(filename=kpointsfile, interpolate=interpolate)
    nks = length(kpoints)
    nωs = length(ωs)    
    plasmon = Array{Float64, 2}(undef, (nωs, nks ))
    for (qidx, q) in enumerate(kpoints)
        println(q)
        for (ωidx, ω) in enumerate(ωs)
            plasmon[ωidx, qidx] = direct_epsilon(HWannier, cell_map,lattice_vectors, q, ω, μ, spin = 1, ϵ = 0.1, normalized=true; kwargs...) 
        end
    end
    return plasmon
end 

"""
$(TYPEDSIGNATURES)
Direct 2D integration for Epsilon with HCubature
"""
function direct_epsilon_cubature(wannier_file::AbstractString, cell_map_file::AbstractString, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, 
    ω::Real, μ::Real; spin::Integer = 1, ϵ::Real = 0.01, kwargs...)

    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))
    brillouin_area=brillouin_zone_area(lattice_vectors) 
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand(wannier_file, cell_map_file, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]
    1-e²ϵ/(2qabs)*polarization
end

"""
$(TYPEDSIGNATURES)
"""
function direct_epsilon_cubature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, ω::Real, μ::Real; spin::Integer = 1, ϵ::Real = 0.01, kwargs...)
    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))
    brillouin_area=brillouin_zone_area(lattice_vectors) 
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand(HWannier, cell_map, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]
    1-e²ϵ/(2qabs)*polarization
end

"""
Find the imaginary value of polarization through hcubature 
"""
function im_polarization_cubature(wannier_file::AbstractString, cell_map_file::AbstractString, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, ω::Real, μ::Real; spin::Integer = 1, ϵ::Real=0.01, kwargs...) 
    qnormalized = normalize_kvector(lattice_vectors, q)
    brillouin_area=brillouin_zone_area(lattice_vectors) 
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand_imaginary(wannier_file, cell_map_file, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]
    return polarization
end

function im_polarization_cubature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, ω::Real, μ::Real; spin::Integer = 1, ϵ::Real=0.01, kwargs...) 
    qnormalized = normalize_kvector(lattice_vectors, q)
    brillouin_area=brillouin_zone_area(lattice_vectors) 
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand_imaginary(HWannier, cell_map, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]
    return polarization
end

"""
$(TYPEDSIGNATURES)
returns the non-local, non-static dielectric function
"""
function return_2d_epsilon(q::Vector{<:Real}, lat::Vector{<:Vector{<:Real}},  ω::Real, im_pol::Vector{<:Real}, 
    max_energy::Real, histogram_width::Real, normalized::Bool=true) 

    qabs = normalized ? sqrt(sum(unnormalize_kvector(lat, q).^2)) : sqrt(sum((q.^2)))
    return 1-e²ϵ/abs(2*qabs)*kramers_kronig(ω, im_pol, max_energy, histogram_width)
end

"""
$(TYPEDSIGNATURES)
returns the non-local, non-static dielectric function
"""
function return_2d_epsilons(ωs::AbstractRange{<:Real}, im_pols::Array{<:Real, 2}, lattice::Vector{<:Vector{<:Real}};
    max_energy::Real, histogram_width::Real, kpointsfile::String="bandstruct.kpoints", nlayers::Integer = 1, d::Real=3, interpolate::Integer=1, plotmap::Bool=true) 
    ϵs = zeros(size(im_pols)[1], length(ωs))
    kpoints = bandstructkpoints2q(filename=kpointsfile, interpolate=interpolate)
    for (qidx, qnorm) in enumerate(kpoints)
        print(qnorm)
        impol = im_pols[qidx, :]
        q = sqrt(sum((unnormalize_kvector(lattice, qnorm)).^2))
        for (ωidx, ω) in enumerate(ωs)
            #ϵs[qidx, ωidx] = real(1-e²ϵ/abs(2*q)*kramers_kronig(ω, impol, max_energy, histogram_width))
            ϵs[qidx, ωidx] = return_2d_epsilon(q, ω, impol, max_energy, histogram_width, d, nlayers) 
        end
    end
    plotmap && display(heatmap(transpose(log.(abs.(ϵs))), yticks=(collect(0:(length(ωs))/5:length(ωs)), 
    collect(minimum(ωs):(maximum(ωs)-minimum(ωs))/5:maximum(ωs))), xticks=[], size=(1000, 500)))
    return ϵs
end

"""
$(TYPEDSIGNATURES)
"""
function return_2d_epsilon(q::Real, ω::Real, im_pol::Vector{<:Real}, max_energy::Real, histogram_width::Real, d::Real, num_layers::Integer) 
    epsilon_mat = Array{Float64, 2}(undef, (num_layers, num_layers))
    polarization = real(kramers_kronig(ω, im_pol, max_energy, histogram_width))
    for i in 1:num_layers
        for j in 1:num_layers
            epsilon_mat[i, j] = (i == j ? 1-e²ϵ/abs(2q)*polarization : -exp(-q*d)*e²ϵ/abs(2q)*polarization)
        end
    end
    return det(epsilon_mat)
end

"""
$(TYPEDSIGNATURES)
returns the non-local, non-static dielectric function using scipy functionality
"""
function return_2d_epsilon_scipy(q::Real, ω::Real, im_pol::Vector{<:Real}, max_energy::Real, histogram_width::Real, max_energy_integration::Real) 
    return 1-e²ϵ/abs(2q)*kramers_kronig_scipy(ω, im_pol, max_energy, histogram_width, max_energy_integration)
end

"""
$(TYPEDSIGNATURES)
"""
function return_2d_epsilon_quadgk(q::Real, ω::Real, im_pol::Vector{<:Real}, max_energy::Real, histogram_width::Real, max_energy_integration::Real; δ::Real = 0.1, kwargs... )
    return 1-e²ϵ/(2abs(q))*kramers_kronig_quadgk(ω, im_pol, max_energy, histogram_width, max_energy_integration; δ, kwargs...)  
end