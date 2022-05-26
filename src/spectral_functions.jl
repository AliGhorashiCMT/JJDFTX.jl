"""
$(TYPEDSIGNATURES)
"""
function subsampling(HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, nbands::Integer, μ::Real, histogram_width::Real; mesh=1000)
    Nkfermi = 0 
    for _ in 1:mesh
        k = rand(3)
        eks = wannier_bands(HWannier, cellmap, k, nbands)
        for ek in eks
            if abs(ek-μ)*histogram_width < 0.5
                Nkfermi +=1
            end
        end
    end
    return mesh/Nkfermi
end

"""
$(TYPEDSIGNATURES)
Reproduces the subsampling as defined in Shankar's online notes
"""
function subsampling_gaussian(HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, nbands::Integer, μ::Real, esigma::Real; mesh=1000)
    Nkfermi = 0 
    for _ in 1:mesh
        k = rand(3)
        eks = wannier_bands(HWannier, cellmap, k, nbands)
        for ek in eks
            weight = 1/(esigma*sqrt(2π))*exp(-0.5*((ek-μ)/esigma)^2)
            (abs(weight) > .001/esigma) || continue
            Nkfermi += 1
        end
    end
    return mesh/Nkfermi
end

"""
$(TYPEDSIGNATURES)
For use by the Eliashberg spectral function method above. Returns the density of states per unit volume, without taking spin degeneracy into account.
"""
function dosatmu(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice::Vector{<:Vector{<:Real}}, nbands::Integer, μ::Real; mesh::Integer = 10, histogram_width::Real=3)
    volume = unit_cell_volume(lattice)
    dos = 0 
    for _ in 1:mesh^3
        ϵs = wannier_bands(Hwannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            abs(μ-ϵ)*histogram_width < 0.5 || continue
            dos += histogram_width*(1/mesh)^3*(1/volume)
        end
    end
    return dos
end

function dosatmu(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, μ::Real; mesh::Integer = 10, histogram_width::Real=3)
    #Using multiple dispatch to have a second method that gives dos at mu as 1/ev instead of 1/(ev*volume)
    dos = 0 
    for _ in 1:mesh^3
        ϵs = wannier_bands(Hwannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            (abs(μ-ϵ)*histogram_width) < 0.5 || continue
            dos += histogram_width*(1/mesh)^3
        end
    end
    return dos
end

"""
$(TYPEDSIGNATURES)
"""
function dosatmulorentzian(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice::Vector{<:Vector{<:Real}}, nbands::Integer, μ::Real; mesh::Integer = 10, esmearing::Real=1)
    volume = unit_cell_volume(lattice)
    dos = 0 
    for _ in 1:mesh^3
        ϵs = wannier_bands(Hwannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            dos += abs(1/π*imag(1/(ϵ-μ+1im*esmearing)))*(1/mesh)^3*(1/volume)
        end
    end
    return dos
end

"""
$(TYPEDSIGNATURES)
"""
function dosatmugaussian(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice::Vector{<:Vector{<:Real}}, nbands::Integer, μ::Real; mesh::Integer = 10, esmearing::Real=1)
    volume = unit_cell_volume(lattice)
    dos = 0 
    for _ in 1:mesh^3
        ϵs = wannier_bands(Hwannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            dos += 1/(esmearing*sqrt(2*π))*exp(-0.5*((ϵ-μ)/esmearing)^2)*(1/mesh)^3*(1/volume)
        end
    end
    return dos
end

"""
$(TYPEDSIGNATURES)
"""
function vFsquaredatmu(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, nbands::Integer, μ::Real; mesh::Integer=10, histogram_width::Real=3)
    vFsquared = 0 
    numintersections = 0 
    for _ in 1:mesh^3
        arandk = rand(3)
        ϵs = wannier_bands(Hwannier, cell_map, arandk, nbands)
        ps = momentum_matrix_elements(Hwannier, cell_map, Pwannier, arandk)
        for (index, ϵ) in enumerate(ϵs)
            abs(μ-ϵ)*histogram_width < 0.5 || continue
            numintersections +=1
            vFsquared = vFsquared + sum((abs.(ps[:, index, index])).^2)*(bohrtoangstrom/ħ)^2
                ##Note that to stay in keeping with JDFTX conventions, we reconverted to atomic units
        end
    end
    averaged_fermivelocity = sqrt(vFsquared/numintersections)
    isnan(averaged_fermivelocity)  ? println("Got NaN- which typically means you need more sampling points. Try increasing mesh.") : nothing
    println("NkFermi= ", numintersections)
    println("Subsampling= ", numintersections/(mesh^3))
    println("Inverse Subsampling = ", mesh^3/numintersections)
    return sqrt(vFsquared/numintersections)
end

# Conventions of this package are that the dos will be in 1/angstrom^3*1/eV units, to convert to jdftx units, 
# we must do the following:
convertdos(dos::Real) = dos*bohrtoangstrom^3*1/eV

#Todo Fix units in resistivity
"""
$(TYPEDSIGNATURES)
Return the resistivity from the eliashberg spectral function when given as an array    
"""
function eliashbergresistivity(eliashbergarray::Vector{<:Real}, maxenergy::Real, T::Real, vF::Real, Volume::Real)
    ##First find all relevant quantities in prefactor
    conductancequantum = 7.75*1e-5 ##In Siemens
    electronmass = 9.109*1e-31 ##In kilograms
    vF *= 2e-24*1/electronmass
    vF *= 1e10 ##In angstroms per second
    TempEV = kB*T
    prefactor = 6/(ħ^2*conductancequantum*vF^2)*1/Volume
    println("prefactor is:", prefactor)
    xarray = collect(maxenergy/length(eliashbergarray):maxenergy/length(eliashbergarray):maxenergy)
    interpolated_eliashberg = interpol.interp1d(xarray, eliashbergarray)
    temperature_weights = (xarray ./ TempEV).*exp.(xarray ./ TempEV) ./ (exp.(xarray ./ TempEV) .- 1).^2 ##Temperature weights 
    temperature_weighted_eliashberg = interpol.interp1d(xarray, eliashbergarray.*temperature_weights)
    return prefactor/10*pyintegrate.quad(temperature_weighted_eliashberg, xarray[1], xarray[end])[1]
end

"""
$(TYPEDSIGNATURES)
"""
function eliashbergresistivities(eliashbergarray::Vector{<:Real}, maxenergy::Real, Ts::Vector{<:Real}, vF::Real, Volume::Real)
    ρs = Float64[]
    for T in Ts
        push!(ρs, eliashbergresistivity(eliashbergarray, maxenergy, T, vF, Volume))
    end
    return ρs
end