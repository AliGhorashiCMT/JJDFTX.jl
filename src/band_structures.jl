"""
$(TYPEDSIGNATURES)
Returns an array of k points in the basis of reciprocal lattice vectors, with optional interpolation given by 
keyword argument interpolate. The k points are by default read from a file in the current directory given by 
bandstruct.kpoints, in keeping with JDFTX conventions. This can be changed by passing the keyword argument filename.
"""
function bandstructkpoints2q(; kpointsfile::AbstractString="bandstruct.kpoints", interpolate::Integer=1)
    kpointlist = np.loadtxt(kpointsfile, skiprows=2, usecols=[1, 2, 3])
    numpoints = np.shape(kpointlist)[1]
    kpointsreshaped = Vector{Float64}[]
    for k in 1:numpoints-1
        for interpolatedk in 0:interpolate-1
            push!(kpointsreshaped, kpointlist[k, :] .+ (kpointlist[k+1, :].-kpointlist[k, :])./interpolate.*interpolatedk)
        end
    end
    push!(kpointsreshaped, kpointlist[numpoints, :])
    return kpointsreshaped
end

"""
$(TYPEDSIGNATURES)
Label high symmetry k points on a band structure plot
"""
function label_plots(kticksfile::AbstractString="bandstruct.kpoints.in", kpointsfile::AbstractString="bandstruct.kpoints",
    to_greek::Bool=true)

    kpointscoords=Vector{Float64}[]
    kpointslabels=String[]
    if isfile(kticksfile)
        for line in readlines(kticksfile)
            try
                line_info = split(line)
                push!(kpointslabels, line_info[end])
                push!(kpointscoords, parse.(Float64, line_info[2:4]))
            catch
            end
        end
        xtickindices=Integer[]
        xticklabels=String[]
        for (tick, line) in enumerate(readlines(kpointsfile)[3:end])
            for (kplabel, kpcoord) in zip(kpointslabels, kpointscoords)
                kpointcoord=parse.(Float64, split(line)[2:4])
                isapprox(kpointcoord, kpcoord) || continue
                push!(xtickindices, tick-1) # Convention of -1 since plots start at x = 0 in PyPlot
                push!(xticklabels, kplabel)
                break
            end
        end
        to_greek && replace!(xticklabels, "Gamma" => "Γ")
        xticks(xtickindices, xticklabels)
        xlim(0, xtickindices[end])
    end
end

function load_bandeigs_data(band_file::AbstractString, numpoints::Integer, numbands::Integer, spin::Integer=1)
    energies = np.reshape(np.fromfile(band_file), (numpoints*spin, numbands))*1/eV
    return energies
end

function load_bands_points(band_file::AbstractString, kpointsfile::AbstractString="bandstruct.kpoints", spin::Integer=1)
    numpoints = countlines(kpointsfile) - 2  
    numeigenvals = length(np.fromfile(band_file))
    numbands = convert(Integer, numeigenvals/(numpoints*spin))
    return numpoints, numbands
end

function bandstruct_properties(band_file::AbstractString, numbands::Integer, numpoints::Integer; 
    spin::Integer=1)
    energies = load_bandeigs_data(band_file, numpoints, numbands, spin)
    if spin == 1
        energies_min_max = [(emin, emax) for (emin, emax) in zip(vec(minimum(energies, dims=1)), vec(maximum(energies, dims=1)))]
        return energies_min_max
    elseif spin ==2 
        energies_up = energies[1:numpoints, :]
        energies_dn = energies[numpoints+1:end, :]
        energies_up_min_max = [(emin, emax) for (emin, emax) in zip(vec(minimum(energies_up, dims=1)), vec(maximum(energies_up, dims=1)))]
        energies_dn_min_max = [(emin, emax) for (emin, emax) in zip(vec(minimum(energies_dn, dims=1)), vec(maximum(energies_dn, dims=1)))]
        return energies_up_min_max, energies_dn_min_max
    end
end

function bandstruct_properties(band_file::AbstractString; kpointsfile::AbstractString="bandstruct.kpoints",
    spin::Integer=1)
    
    numpoints, numbands = load_bands_points(band_file, kpointsfile, spin)
    bandstruct_properties(band_file, numbands, numpoints, spin=spin;)
end

"""
$(TYPEDSIGNATURES)

Plots the bands from a non self consistent calculation. First argument must be the file with 
the corresponding band eigenvalues. num_bands is the number of bands included in the calculation. Note
that the spin degeneracy in jdftx is included in the number of k points- not the number of bands. Therefore, 
the k points from 1:num_points will be for one spin species and those from num_points+1 to 2*npoints
correspond to the other spin species.
"""
function plot_bands(band_file::AbstractString, numbands::Integer, numpoints::Integer; 
    whichbands::Union{Nothing, Vector{<:Integer}}=nothing, spin::Integer=1, color_up::AbstractString = "blue", 
    color_dn::AbstractString = "red", color_nospin::AbstractString = "black", kwargs...)
    
    energies = load_bandeigs_data(band_file, numpoints, numbands, spin)
    if spin == 1
        isnothing(whichbands) ? plot(energies, color=color_nospin, label="", linewidth=2; kwargs...) : plot(energies[:, whichbands], color=color_nospin, label="", linewidth=2; kwargs...)
    elseif spin ==2 
        energies_up = energies[1:numpoints, :]
        energies_dn = energies[numpoints+1:end, :]
        isnothing(whichbands) ? plot(energies_up, color=color_up, label="", linewidth=5; kwargs...) : plot(energies_up[:, whichbands], color=color_up, label="", linewidth=5; kwargs...)
        isnothing(whichbands) ? plot(energies_dn, color=color_dn, label="", linewidth=5; kwargs...) : plot(energies_dn[:, whichbands], color=color_dn, label="", linewidth=5; kwargs...)
    end
    ylabel("Energy (eV)")
    xlabel("Wavevector")
end

function plot_bands(band_file::AbstractString; kpointsfile::AbstractString="bandstruct.kpoints",
    kticksfile="bandstruct.kpoints.in", spin::Integer=1, whichbands::Union{Nothing, Vector{<:Integer}}=nothing, 
    to_greek::Bool=true, color_up::AbstractString = "blue", color_dn::AbstractString = "red", 
    color_nospin::AbstractString = "black", kwargs...)

    numpoints, numbands = load_bands_points(band_file, kpointsfile, spin)
    plot_bands(band_file, numbands, numpoints, whichbands=whichbands, color_up = color_up, color_dn = color_dn, color_nospin = color_nospin, 
    spin=spin; kwargs...)
    label_plots(kticksfile, kpointsfile, to_greek)
end

"""
$(TYPEDSIGNATURES)

Plot the Wannier band structure along a kpoints path provided through a file written in JDFTX bandstruct.kpoints conventions
"""
function plot_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}; 
    whichbands::Union{Vector{<:Integer}, Nothing}=nothing, kpointsfile::AbstractString="bandstruct.kpoints",
    kwargs...)

    kpoints = bandstructkpoints2q(kpointsfile=kpointsfile) 
    energies = first(wannier_bands(Hwannier, cell_map, hcat(kpoints...)))
    isnothing(whichbands) ? plot(energies; kwargs...) : plot(energies[:, whichbands]; kwargs...)
end

"""
$(TYPEDSIGNATURES)
"""
function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Vector{<:Real}) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); Es, Us = np.linalg.eigh(H);
    return Es./eV, Us 
end

"""
$(TYPEDSIGNATURES)
"""
function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, kpoints::AbstractArray{<:Real, 2}) 
    Es, Us = np.linalg.eigh(np.einsum("kij, kl -> lij", Hwannier, np.exp(2im*np.pi*cell_map*kpoints)))
    return Es./eV, Us
end

"""
$(TYPEDSIGNATURES)
"""
function hwannier(wannier_file::AbstractString, cell_map_file::AbstractString) 
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = np.loadtxt(wannier_file)
    numbands = Int(sqrt(size(Hwannier)[2]))
    Hwannier = permutedims(reshape(Hwannier, (cell_map_numlines, numbands, numbands)), [1, 3, 2])
    return Hwannier
end

hwannier(filebase::AbstractString) = hwannier("$filebase.txt", "$filebase.map.txt")


