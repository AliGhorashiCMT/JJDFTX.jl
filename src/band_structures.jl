
function label_plots(kticksfile::AbstractString = "bandstruct.kpoints.in", kpointsfile::AbstractString="bandstruct.kpoints", to_greek::Bool=true)
    kpointscoords=Vector{Vector{Float64}}()
    kpointslabels=Vector{String}()
    if isfile(kticksfile)
        for line in readlines(kticksfile)
            try
                push!(kpointslabels, split(line)[end])
                push!(kpointscoords, parse.(Float64, split(line)[2:4]))
            catch
            end
        end
        xtickindices=Vector{Integer}()
        xticklabels=Vector{String}()
        for (tick, line) in enumerate(readlines(kpointsfile)[3:end])
            for (kplabel, kpcoord) in zip(kpointslabels, kpointscoords)
                kpointcoord=parse.(Float64, split(line)[2:4])
                isapprox(kpointcoord, kpcoord) || continue
                push!(xtickindices, tick-1)
                push!(xticklabels, kplabel)
                break
            end
        end
        to_greek && replace!(xticklabels, "Gamma" => "Γ")
        xticks(xtickindices, xticklabels)
    end
end


function bandstruct_properties(band_file::AbstractString, num_bands::Integer, num_points::Integer; 
    spin::Integer=1, kwargs...)
    if spin == 1
        energies = np.reshape(np.fromfile(band_file), (num_points, num_bands))*1/eV;
        energies_min_max = [(min, max) for (min, max) in zip(vec(minimum(energies, dims=1)), vec(maximum(energies, dims=1)))]
        return energies_min_max
    elseif spin ==2 
        energies = np.reshape(np.fromfile(band_file), (num_points*2, num_bands))*1/eV;
        energies_up = energies[1:num_points, :]
        energies_dn = energies[num_points+1:end, :]
        energies_up_min_max = [(min, max) for (min, max) in zip(vec(minimum(energies_up, dims=1)), vec(maximum(energies_up, dims=1)))]
        energies_dn_min_max = [(min, max) for (min, max) in zip(vec(minimum(energies_dn, dims=1)), vec(maximum(energies_dn, dims=1)))]
        return energies_up_min_max, energies_dn_min_max
    end
end

function bandstruct_properties(band_file::AbstractString; kpointsfile::AbstractString="bandstruct.kpoints",
    kticksfile="bandstruct.kpoints.in", spin::Integer=1, whichbands::Union{Nothing, Vector{<:Integer}}=nothing, to_greek::Bool=false, kwargs...)
    numpoints = countlines(kpointsfile) - 2  
    numeigenvals = length(np.fromfile(band_file))
    numbands = convert(Integer, numeigenvals/(numpoints*spin))
    bandstruct_properties(band_file, numbands, numpoints, spin=spin; kwargs...)
end


"""
$(TYPEDSIGNATURES)

Plots the bands from a non self consistent calculation. First argument must be the file with 
the corresponding band eigenvalues. num_bands is the number of bands included in the calculation. Note
that the spin degeneracy in jdftx is included in the number of k points- not the number of bands. Therefore, 
the k points from 1:num_points will be for one spin species and those from num_points+1 to 2*npoints
correspond to the other spin species.
"""
function plot_bands(band_file::AbstractString, num_bands::Integer, num_points::Integer; 
    whichbands::Union{Nothing, Vector{<:Integer}}=nothing, spin::Integer=1, kwargs...)
    if spin == 1
        energies = np.reshape(np.fromfile(band_file), (num_points, num_bands))*1/eV;
        isnothing(whichbands) ? plot(energies, color="red", label="", linewidth=2; kwargs...) : plot(energies[:, whichbands], color="red", label="", linewidth=2; kwargs...)
    elseif spin ==2 
        energies = np.reshape(np.fromfile(band_file), (num_points*2, num_bands))*1/eV;
        energies_up = energies[1:num_points, :]
        energies_dn = energies[num_points+1:end, :]
        isnothing(whichbands) ? plot(energies_up, color="black", label="", linewidth=5; kwargs...) : plot(energies_up[:, whichbands], color="black", label="", linewidth=5; kwargs...)
        isnothing(whichbands) ? plot(energies_dn, color="red", label="", linewidth=5; kwargs...) : plot(energies_dn[:, whichbands], color="red", label="", linewidth=5; kwargs...)
    end
end

function plot_bands(band_file::AbstractString; kpointsfile::AbstractString="bandstruct.kpoints",
    kticksfile="bandstruct.kpoints.in", spin::Integer=1, whichbands::Union{Nothing, Vector{<:Integer}}=nothing, to_greek::Bool=false, kwargs...)
    numpoints = countlines(kpointsfile) - 2  
    numeigenvals = length(np.fromfile(band_file))
    numbands = convert(Integer, numeigenvals/(numpoints*spin))
    plot_bands(band_file, numbands, numpoints, whichbands=whichbands, spin=spin; kwargs...)
    ylabel("Energy (eV)")
    xlabel("Wavevector")
    kpointscoords=Vector{Vector{Float64}}()
    kpointslabels=Vector{String}()
    if isfile(kticksfile)
        for line in readlines(kticksfile)
            try
                push!(kpointslabels, split(line)[end])
                push!(kpointscoords, parse.(Float64, split(line)[2:4]))
            catch
            end
        end
        xtickindices=Vector{Integer}()
        xticklabels=Vector{String}()
        for (tick, line) in enumerate(readlines(kpointsfile)[3:end])
            for (kplabel, kpcoord) in zip(kpointslabels, kpointscoords)
                kpointcoord=parse.(Float64, split(line)[2:4])
                isapprox(kpointcoord, kpcoord) || continue
                push!(xtickindices, tick-1)
                push!(xticklabels, kplabel)
                break
            end
        end
        to_greek && replace!(xticklabels, "Gamma" => "Γ")
        xticks(xtickindices, xticklabels)
    end
end



"""
$(TYPEDSIGNATURES)
Plot several band structures overlayed on one another (assuming they take the same path through kspace)
"""
function plotmanybands(kpoints::AbstractString, bandfiles::Vector{<:AbstractString}, spin::Val{2}; shifts::Union{Vector{<:Real}, Nothing}=nothing, 
    μs::Union{Vector{<:Real}, Nothing}=nothing, whichbands::Vector{<:Integer}=Int[], kwargs...)
    plotly()
    numkpoints = size(np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3]))[1] ##Get number of kpoints at which bands are evaluated
    numbandfiles = length(bandfiles)
    newshifts = shifts isa Nothing ? zeros(numbandfiles) : shifts ##Take into account possiblity of no shifts
    colors = collect(1:numbandfiles)
    numbandseach = Int.(first.(np.shape.(np.fromfile.(bandfiles))) ./(2*numkpoints)) ##Find number of bands for each file
    if isempty(whichbands)
        for (index, (numbands, bandfile, shift)) in enumerate(zip(numbandseach, bandfiles, newshifts))
            reshaped=reshape(read!(bandfile, Array{Float64}(undef, numbands*numkpoints*2)),(numbands, numkpoints*2)) .+ shift*eV;
            exactenergiesup=permutedims(reshaped, [2, 1])[1:numkpoints, :]*1/eV;
            exactenergiesdown=permutedims(reshaped, [2, 1])[numkpoints+1:2*numkpoints, :]*1/eV;
            index == 1 ? display(plot(exactenergiesup, color = colors[index], size=(1000, 500), legend=false; kwargs...)) : display(plot!(exactenergiesup, color = colors[index], size=(1000, 500), legend=false; kwargs...) )
            display(plot!(exactenergiesdown, color = colors[index], size=(1000, 500), legend=false; kwargs...)) 
        end
    else
        for (index, (numbands, bandfile, shift)) in enumerate(zip(numbandseach, bandfiles, newshifts))
            reshaped=reshape(read!(bandfile, Array{Float64}(undef, numbands*numkpoints*2)),(numbands, numkpoints*2)) .+ shift*eV;
            exactenergiesup=permutedims(reshaped, [2, 1])[1:numkpoints, whichbands]*1/eV;
            exactenergiesdown=permutedims(reshaped, [2, 1])[numkpoints+1:2*numkpoints, whichbands]*1/eV;
            index == 1 ? display(plot(exactenergiesup, color = colors[index], size=(1000, 500), legend=false; kwargs...)) : display(plot!(exactenergiesup, color = colors[index], size=(1000, 500), legend=false; kwargs...) )
            display(plot!(exactenergiesdown, color = colors[index], size=(1000, 500), legend=false; kwargs...)) 
        end
    end
    ylabel!("Energy (eV)", yguidefontsize=20)
    yticks!(round.(collect(ylims()[1]:(ylims()[2]-ylims()[1])/10:ylims()[2]), digits=2);ytickfontsize=20)
    xticks!(Float64[])
    if μs isa Vector{<:Real}
        for i in 1:numbandfiles
            display(hline!([μs[i]+newshifts[i]], linewidth=2, color=colors[i]))
        end
    end
end

function plotmanybands(kpoints::AbstractString, bandfiles::Vector{<:AbstractString}, spin::Val{1}; shifts::Union{Vector{<:Real}, Nothing}=nothing, 
    μs::Union{Vector{<:Real}, Nothing}=nothing, whichbands::Vector{<:Integer}=Int[], kwargs...)
    plotly()
    numkpoints = size(np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3]))[1] ##Get number of kpoints at which bands are evaluated
    numbandfiles = length(bandfiles)
    newshifts = shifts isa Nothing ? zeros(numbandfiles) : shifts ##Take into account possiblity of no shifts
    colors = collect(1:numbandfiles)
    numbandseach = Int.(first.(np.shape.(np.fromfile.(bandfiles))) ./numkpoints) ##Find number of bands for each file
    if isempty(whichbands)
        for (index, (numbands, bandfile, shift)) in enumerate(zip(numbandseach, bandfiles, newshifts))
            index == 1 ? display(plot(np.reshape(np.fromfile(bandfile)*1/eV .+ shift, (numkpoints, numbands)), color = colors[index], size=(1000, 500), legend=false; kwargs...)) : display(plot!(np.reshape(np.fromfile(bandfile)*1/eV .+ shift, (numkpoints, numbands)), color = colors[index], size=(1000, 500), legend=false; kwargs...) )
        end
    else
        for (index, (numbands, bandfile, shift)) in enumerate(zip(numbandseach, bandfiles, newshifts))
            index == 1 ? display(plot(np.reshape(np.fromfile(bandfile)*1/eV .+ shift, (numkpoints, numbands))[:, whichbands], color = colors[index], size=(1000, 500), legend=false; kwargs...)) : display(plot!(np.reshape(np.fromfile(bandfile)*1/eV .+ shift, (numkpoints, numbands))[:, whichbands], color = colors[index], size=(1000, 500), legend=false;kwargs...))
        end
    end
    ylabel!("Energy (eV)", yguidefontsize=20)
    yticks!(round.(collect(ylims()[1]:(ylims()[2]-ylims()[1])/10:ylims()[2]), digits=2);ytickfontsize=20)
    xticks!(Float64[])
    if μs isa Vector{<:Real}
        for i in 1:numbandfiles
            display(hline!([μs[i]+newshifts[i]], linewidth=2, color=colors[i]))
        end
    end
end


"""
$(TYPEDSIGNATURES)

Plot the Wannier band structure along a kpoints path provided through a file written in JDFTX bandstruct.kpoints conventions
"""
function plotwannierbands(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer; whichbands::Union{Vector{<:Integer}, Nothing}=nothing,
    kpoints::AbstractString="bandstruct.kpoints", overlay::Bool=false, kwargs...)
    kpointlist = np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    energiesatkpoints = Array{Float64, 2}(undef, (num_kpoints, nbands))
    for k in 1:num_kpoints
        energiesatkpoints[k, :] = wannier_bands(HWannier, cell_map, kpointlist[k, :], nbands)
    end
    if overlay
        isnothing(whichbands) ? plot(energiesatkpoints; kwargs...) : plot(energiesatkpoints[:, whichbands]; kwargs...)
    else
        isnothing(whichbands) ? plot(energiesatkpoints; kwargs...) : plot(energiesatkpoints[:, whichbands]; kwargs...)
    end
end

"""
$(TYPEDSIGNATURES)

Plot the wannier bands along the supplied k point path. 
"""
function plotwannierbands(HWannierUp::Array{Float64, 3}, HWannierDn::Array{Float64, 3}, cellmapUp::Array{Float64, 2}, 
    cellmapDn::Array{Float64, 2}, nbands::Integer; whichbands::Union{Vector{<:Integer}, Nothing}=nothing, kpoints::AbstractString="bandstruct.kpoints", kwargs...)
    kpointlist = np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    energiesatkpointsUp = Array{Float64, 2}(undef, (num_kpoints, nbands))
    energiesatkpointsDn = Array{Float64, 2}(undef, (num_kpoints, nbands))
    for k in 1:num_kpoints
        energiesatkpointsUp[k, :] = wannier_bands(HWannierUp, cellmapUp, kpointlist[k, :], nbands)
        energiesatkpointsDn[k, :] = wannier_bands(HWannierDn, cellmapDn, kpointlist[k, :], nbands)
    end
    isnothing(whichbands) ? plot(energiesatkpointsUp,  color="green"; kwargs...) : plot(energiesatkpointsUp[:, whichbands], color="green"; kwargs...)
    isnothing(whichbands) ? plot(energiesatkpointsDn,  color="orange"; kwargs...) : plot(energiesatkpointsDn[:, whichbands], color="orange"; kwargs...)
end


"""
$(TYPEDSIGNATURES)
"""
function plotbandsoverlayedwannier(band_file::AbstractString, ntotalbands::Integer, HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    nwannierbands::Integer, numpoints::Integer; spin::Integer=1, kpoints::AbstractString="bandstruct.kpoints", kwargs...)
    plot_bands(band_file, ntotalbands, numpoints, spin=spin; kwargs...)
    plotwannierbands(HWannier, cell_map, nwannierbands, kpoints=kpoints; linestyle = :dashdot, kwargs... )
end

"""
$(TYPEDSIGNATURES)
"""
function plotbandsoverlayedwannier(band_file::AbstractString, HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    nwannierbands::Integer=2; spin::Integer=1, kpointsfile::AbstractString="bandstruct.kpoints", kwargs...)
    plot_bands(band_file, kpointsfile=kpointsfile, spin=spin; linewidth=1, kwargs...)
    plotwannierbands(HWannier, cell_map, nwannierbands, kpoints=kpointsfile; overlay=true,  linewidth=5, linestyle = :dashdot, kwargs... )
end

"""
$(TYPEDSIGNATURES)
"""
function wannier_bands(wannier_file::AbstractString, cell_map_file::AbstractString, k::Vector{<:Real}) 
    cell_map=np.loadtxt(cell_map_file)
    cell_map_numlines=countlines(cell_map_file)
    #=
    Note that Julia's reshape method is different from that employed by numerical python. Hence, we must permute the indices of the 
    bands in order to recover the same data that was stored in the text files
    =#
    Hwannier=permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, 1, 1)), [1, 3, 2])
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, _ =np.linalg.eigh(H);
    return E[1]/eV 
end

"""
$(TYPEDSIGNATURES)
"""
function wannier_bands(wannier_file::AbstractString, cell_map_file::AbstractString, k::Vector{<:Real}, nbands::Integer) 
    cell_map=np.loadtxt(cell_map_file)
    cell_map_numlines=countlines(cell_map_file)
    Hwannier=permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, nbands, nbands)), [1, 3, 2])
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E/eV 
end

"""
$(TYPEDSIGNATURES)

Returns the wannier energy dispersion at the supplied k point in eV. Note that JDFTX provides energies in Hartree so an explicit
conversion takes place. 
"""
function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Vector{<:Real}) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, _=np.linalg.eigh(H);
    return E[1]./eV 
end

"""
$(TYPEDSIGNATURES)
"""
function wannier_bands(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Vector{<:Real}, nbands::Integer) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E./eV 
end

function wannier_vectors(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k::Vector{<:Real}) 
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return U
end

function wannier_vectors(wannier_file::AbstractString, cell_map_file::AbstractString, k::Vector{<:Real}, nbands::Integer) 
    cell_map = np.loadtxt(cell_map_file)
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, nbands, nbands)), [1, 3, 2])
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return U
end

"""
$(TYPEDSIGNATURES)
"""
function hwannier(wannier_file::AbstractString, cell_map_file::AbstractString, nbands::Integer) 
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, nbands, nbands)), [1, 3, 2])
    return Hwannier
end

"""
$(TYPEDSIGNATURES)
"""
function hwannier(wannier_file::AbstractString, cell_map_file::AbstractString) 
    cell_map_numlines = countlines(cell_map_file)
    Hwannier = permutedims(reshape(np.loadtxt(wannier_file), (cell_map_numlines, 1, 1)), [1, 3, 2])
    return Hwannier
end


## Band Properties: Widths, Max, cell_map_numlines

