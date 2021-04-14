using Plots
using PyCall
using LinearAlgebra 
using Distances

"""
Plots the bands from a non self consistent calculation. First argument must be the file with 
the corresponding band eigenvalues. num_bands is the number of bands included in the calculation. Note
that the spin degeneracy in jdftx is included in the number of k points- not the number of bands. Therefore, 
the k points from 1:num_points will be for one spin species and those from num_points+1 to 2*npoints
correspond to the other spin species.
"""
function plot_bands(band_file::String, num_bands::Int, num_points::Int; spin::Int=1, kwargs...)
    if spin == 1
        reshaped = reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points )),(num_bands, num_points));
        exactenergies = permutedims(reshaped, [2, 1])*1/eV;
        plot(exactenergies, color="black", label="", linewidth=2; kwargs...)
    elseif spin ==2 
        reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points*2 )),(num_bands, num_points*2));
        exactenergiesup=permutedims(reshaped, [2, 1])[1:num_points, :]*1/eV;
        exactenergiesdown=permutedims(reshaped, [2, 1])[num_points+1:2*num_points, :]*1/eV;

        ##Note that Plots.jl automatically plots 2d arrays columnwise- which is why the band indices now correspond to column indices
        plot(exactenergiesdown, color="black", label="", linewidth=2; kwargs...)
        plot!(exactenergiesup, color="purple", label="", linewidth=2; kwargs...)
    end
end

function wannier_bands(wannier_file::String, cell_map_file::String, k::Array{<:Real, 1}) 
    cell_map=np.loadtxt(cell_map_file)
    cell_map_numlines=countlines(cell_map_file)
    #=
    Note that Julia's reshape method is different from that employed by numerical python. Hence, we must permute the indices of the 
    bands in order to recover the same data that was stored in the text files
    =#
    Hwannier=permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, 1, 1)), [1, 3, 2])
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E[1]/eV 
end

function wannier_bands(wannier_file::String, cell_map_file::String, k::Array{<:Real, 1}, nbands::Int64) 
    cell_map=np.loadtxt(cell_map_file)
    cell_map_numlines=countlines(cell_map_file)
    Hwannier=permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, nbands, nbands)), [1, 3, 2])
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E/eV 
end

function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Array{<:Real, 1}) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E[1]./eV 
end

function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Array{<:Real, 1}, nbands::Int64) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E./eV 
end

function wannier_vectors(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Array{<:Real, 1}) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return U
end

function wannier_vectors(wannier_file::String, cell_map_file::String, k::Array{<:Real, 1}, nbands::Int64) 
    cell_map = np.loadtxt(cell_map_file)
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, nbands, nbands)), [1, 3, 2])
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return U
end

function hwannier(wannier_file::String, cell_map_file::String, nbands::Int64) 
    cell_map = np.loadtxt(cell_map_file)
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, nbands, nbands)), [1, 3, 2])
    return Hwannier
end

function hwannier(wannier_file::String, cell_map_file::String) 
    cell_map = np.loadtxt(cell_map_file)
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, 1, 1)), [1, 3, 2])
    return Hwannier
end