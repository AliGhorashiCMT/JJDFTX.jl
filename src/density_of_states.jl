"""
$(TYPEDSIGNATURES)
Overlay the bandstructure with the density of states.
"""
function bandsoverlayedDOS(dosfile::AbstractString, band_file::AbstractString, num_bands::Integer, num_points::Integer, 
    energy_range::Tuple{<:Real, <:Real}; spin::Integer=1)
    dosdata = try 
        np.loadtxt(dosfile)
    catch 
        np.loadtxt(dosfile, skiprows=1)
    end
    if spin == 2
        reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points*2 )),(num_bands, num_points*2));
        exactenergiesup=permutedims(reshaped, [2, 1])[1:num_points, :]*1/eV;
        exactenergiesdown=permutedims(reshaped, [2, 1])[num_points+1:2*num_points, :]*1/eV;
        plot(exactenergiesdown, color="black", label="", linewidth=2, ylims = collect(energy_range))
        B = plot!(exactenergiesup, color="purple", label="", linewidth=2, ylabel = "Energy (eV)", ylims = collect(energy_range), xticks=false)    
        @assert isapprox(num_bands, sum(diff(dosdata[:, 1]).*dosdata[2:end, 2]), atol=1e-1) "DOS not propertly normalized. Make sure files are correct"
    elseif spin == 1
        reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points)),(num_bands, num_points));
        exactenergies=permutedims(reshaped, [2, 1])[1:num_points, :]*1/eV;
        B = plot(exactenergies, color="purple", label="", linewidth=2, ylabel = "Energy (eV)", ylims = collect(energy_range), xticks=false)    
        @assert isapprox(num_bands*2, sum(diff(dosdata[:, 1]).*dosdata[2:end, 2]), atol=1e-1) "DOS not propertly normalized. Make sure files are correct"
    end
    lowerDOS = argmin(abs.(dosdata[:, 1]*1/eV .- energy_range[1]))
    upperDOS = argmin(abs.(dosdata[:, 1]*1/eV .- energy_range[2]))
    C = plot(dosdata[:, 2]*eV, dosdata[:, 1]*1/eV, linewidth=2, ylims = collect(energy_range), xlims = [0, maximum((dosdata[:, 2]*eV)[lowerDOS:upperDOS])], xlabel = "DOS (1/eV)")
    plot(B, C, size = (1000, 500), legend = false)
end

function bandsoverlayedDOS(dosfile::AbstractString, bandfile::AbstractString; energyrange::Tuple{<:Real, <:Real}=(-10, 10), kpointsfile::AbstractString="bandstruct.kpoints", spin::Integer=1 )
    numpoints = countlines(kpointsfile) - 2  
    numeigenvals = length(np.fromfile(bandfile))
    numbands = convert(Integer, numeigenvals/(spin*numpoints))
    bandsoverlayedDOS(dosfile, bandfile, numbands, numpoints, energyrange, spin=spin)
end

"""
$(TYPEDSIGNATURES)

For z-spin calculations where the density of states is outputed as spinUp and spinDn. The number of kpoints and the number of bands must 
be explicitly given. If not given, the method looks for bandstruct.kpoints to determine the number of points. This is multiplied by 2 due to the 
spin polarization and the number of bands is taken to be the size of the eigenvalues file divided by the number of kpoints. 

## Args

`dosfile1` : spin up DOS file (file extension dosUp)

`dosfile2` : spin dn DOS file (file extension dosDn)

`bandfile` : eigenvalues binary file (file extension .eigenvals)

`energy_range` : range of energies to be plotted


"""
function bandsoverlayedDOS2(dosfile1::AbstractString, dosfile2::AbstractString, band_file::AbstractString, num_bands::Integer, 
    num_points::Integer, energy_range::Tuple{<:Real, <:Real}; color_up="blue", color_dn="red", label_plot::Bool=true,
    kticksfile::AbstractString="bandstruct.kpoints.in", kpointsfile::AbstractString="bandstruct.kpoints", return_tot::Bool=false,
    band_subplot::Vector{<:Int}=[1, 2, 1], dos_subplot::Vector{<:Int}=[1, 2, 2], dos_yticks::Bool=true, kwargs...)

    energies = np.reshape(np.fromfile(band_file), (num_points*2, num_bands))*1/eV
    energies_up = energies[1:num_points, :]
    energies_dn = energies[num_points+1:end, :]
    subplot(band_subplot...)
    plot(energies_up, color=color_up, label="", linewidth=2; kwargs...)
    ylim(collect(energy_range))
    plot(energies_dn, color=color_dn, label="", linewidth=2; kwargs...)
    ylabel("Energy (eV)")

    label_plot && label_plots(kticksfile, kpointsfile)
    dosdata1 = try 
        np.loadtxt(dosfile1)
    catch 
        np.loadtxt(dosfile1, skiprows=1)
    end
    dosdata2 = try 
         np.loadtxt(dosfile2)
    catch 
        np.loadtxt(dosfile2, skiprows=1)
    end

    lowerDOS1 = argmin(abs.(dosdata1[:, 1]*1/eV .- energy_range[1]))
    upperDOS1 = argmin(abs.(dosdata1[:, 1]*1/eV .- energy_range[2]))
    lowerDOS2 = argmin(abs.(dosdata2[:, 1]*1/eV .- energy_range[1]))
    upperDOS2 = argmin(abs.(dosdata2[:, 1]*1/eV .- energy_range[2]))
    max1 = maximum((dosdata1[:, 2]*eV)[lowerDOS1:upperDOS1])
    max2 = maximum((dosdata2[:, 2]*eV)[lowerDOS2:upperDOS2])
    max = maximum([max1, max2])

    subplot(dos_subplot...)
    plot(dosdata1[:, 2]*eV, dosdata1[:, 1]*1/eV, linewidth=2, color=color_up; kwargs...)
    plot(dosdata2[:, 2]*eV, dosdata2[:, 1]*1/eV, linewidth=2, color=color_dn; kwargs...)
    !dos_yticks && yticks(Float64[])

    return_tot  && println("Total number of electrons in range: ", sum(diff(dosdata1[:, 1])[lowerDOS1:upperDOS1] .* (dosdata1[:, 2])[lowerDOS1:upperDOS1])+sum(diff(dosdata2[:, 1])[lowerDOS2:upperDOS2] .* (dosdata2[:, 2])[lowerDOS2:upperDOS2] ))
    xlabel("DOS (1/eV)")
    xlim(0, max)
    ylim(collect(energy_range))

end

function bandsoverlayedDOS2(dosfile1::AbstractString, dosfile2::AbstractString, band_file::AbstractString, energy_range::Tuple{<:Real, <:Real}=(-100, 100); kwargs...)
    numpoints = countlines("bandstruct.kpoints") - 2  
    numeigenvals = length(np.fromfile(band_file))
    numbands = convert(Integer, numeigenvals/(numpoints*2))
    bandsoverlayedDOS2(dosfile1, dosfile2, band_file, numbands, numpoints, energy_range; kwargs...)
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states(dosfile_1::AbstractString, dosfile_2::AbstractString; color_up::AbstractString="blue", 
    color_dn::AbstractString="red", kwargs... )
    parseddos1 = try
        np.loadtxt(dosfile_1)
    catch
        np.loadtxt(dosfile_1, skiprows=1)
    end
    parseddos2 = try
        np.loadtxt(dosfile_2)
    catch
        np.loadtxt(dosfile_2, skiprows=1)
    end
    plot(parseddos1[:, 1]*1/eV, parseddos1[:, 2]*eV, linewidth=4, color=color_up, label="Spin Up"; kwargs...)
    plot(parseddos2[:, 1]*1/eV, parseddos2[:, 2]*eV, linewidth=4, color=color_dn, label="Spin Down"; kwargs...)
    ylabel("Density of States (1/eV/Cell)")
    xlabel("Energy (eV)")
end

function find_chemical_potential(dosfile::AbstractString)
    parseddos = try
        np.loadtxt(dosfile)
    catch
        np.loadtxt(dosfile, skiprows=1)
    end
    x, y = parseddos[:, 1]*1/eV, parseddos[:, 2]*eV
    return x, [sum(y[1:idx-1] .* diff(x[1:idx])) for idx in eachindex(x)[2:end]]
end

"""
$(TYPEDSIGNATURES)


Returns relevant properties of the density of states. Since this is basically only bandgaps, this function returns the energies at which 
the density of states vanishes.
"""
function dos_properties(dosfile_1::AbstractString)
    parseddos = try
        np.loadtxt(dosfile_1)
    catch
        np.loadtxt(dosfile_1, skiprows=1)
    end
    energies, dos = parseddos[:, 1]*1/eV, parseddos[:, 2]*eV
    num_data = length(dos)
    zero_indices = Int[]
    for (idx, density) in enumerate(dos)
        idx == 1 && continue
        idx == num_data && continue
        ((iszero(round(density, digits=3)) && !iszero(round(dos[idx-1], digits=3))) || (iszero(round(density, digits=3)) && !iszero(round(dos[idx+1], digits=3))))  || continue
        push!(zero_indices, idx)
    end
    return energies[zero_indices]
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states(dosfile::AbstractString; returntot::Bool=false, kwargs...)
    parseddos = try
        np.loadtxt(dosfile)
    catch
        np.loadtxt(dosfile, skiprows=1)
    end
    plot(parseddos[:, 1]*1/eV, parseddos[:, 2]*eV, linewidth=4, label="Spin Unpolarized"; kwargs...)
    returntot && println("Total number of electrons is: ", sum(diff(parseddos[:, 1]).*parseddos[2:end, 2]))
    ylabel("Density of States (1/eV/Cell)")
    xlabel("Energy (eV)")
end

"""
The typical density of states outputed by JDFTX is per unit cell. However, sometimes it is more relevant to know the 
density of states per unit volume. This is simply equivalent to dividing the conventional DOS by the unit cell 
volume 
"""
function density_of_states_per_area(dosfile_1::AbstractString, lattice_vecs::Vector{<:Vector{<:Real}}; kwargs...)
    ucell_area = unit_cell_area(lattice_vecs)
    try 
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, 
            size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, 
            size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    end
end

function density_of_states_per_area(dosfile_1::AbstractString, lattice_vecs::lattice; kwargs...)
    ucell_area = unit_cell_area(lattice_vecs)
    try
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    end
end

function density_of_states_per_area(dosfile_1::AbstractString, dosfile_2::AbstractString, lattice_vecs::Vector{<:Vector{<:Real}}; kwargs... )
    ucell_area = unit_cell_area(lattice_vecs)
    try 
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    end
    try
        plot!(np.loadtxt(dosfile_2)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    catch
        plot!(np.loadtxt(dosfile_2, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    end
end

function density_of_states_per_area(dosfile_1::AbstractString, dosfile_2::AbstractString, lattice_vecs::lattice; kwargs... )
    ucell_area = unit_cell_area(lattice_vecs)
    try 
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    end
    try
        plot!(np.loadtxt(dosfile_2)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    catch
        plot!(np.loadtxt(dosfile_2, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    end
end

"""
$(TYPEDSIGNATURES)

Calculate the density of states of a one band wannier interpolation model using the HCubature package. This may be helpful for cross checking values
obtained through histogramming. Note that the limits of integration are 0 to 1 in each direction because we are integrating in the reciprocal lattice basis. 
Therefore, the Jacobian factor cancels out with the other 
"""
function density_of_states_wannier_quad(wannier_file::AbstractString, cell_map_file::AbstractString, ϵ::Real; δ=.1, kwargs...) 
    1/π*hcubature(vec->imag(-1/(ϵ-wannier_bands(wannier_file, cell_map_file, [vec[1], vec[2], 0])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]
end

function density_of_states_wannier_quad(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵ::Real; δ::Real = 0.1, kwargs...)
    1/π*hcubature(vec->imag(-1/(ϵ-wannier_bands(HWannier, cell_map, [vec[1], vec[2], 0])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]
end

function density_of_states_wannier_quad(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵs::Vector{<:Real}; δ::Real = 0.1, kwargs...)
    dos_vec = Float64[]
    for ϵ in ϵs
        push!(dos_vec, density_of_states_wannier_quad(HWannier, cell_map, ϵ, δ=δ;  kwargs...))
    end
    dos_vec
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states_wannier_scipy_quad(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵ::Real; δ::Real = 0.1, kwargs...) 
    nquad = pyintegrate.nquad
    optdict = Dict()
    for kwarg in kwargs
        push!(optdict, kwarg[1]=>kwarg[2])
    end
    1/π*nquad((x, y)->imag(-1/(ϵ-wannier_bands(HWannier, cell_map, [x, y, 0])+1im*δ)), [[0, 1], [0, 1]], opts=optdict)[1]
end

function density_of_states_wannier_scipy_quad(wannier_file::AbstractString, cell_map_file::AbstractString, ϵ::Real; δ::Real = 0.1, kwargs...) 
    nquad = pyintegrate.nquad
    optdict=Dict()
    for kwarg in kwargs
        push!(optdict, kwarg[1]=>kwarg[2])
    end
    1/π*nquad((x, y)->imag(-1/(ϵ-wannier_bands(wannier_file, cell_map_file, [x, y, 0])+1im*δ)), [[0, 1], [0, 1]], opts=optdict)[1]
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states_wannier_quad_check(wannier_file::AbstractString, cell_map_file::AbstractString, ϵmin::Real, ϵmax::Real, numpoints::Integer; δ=.1, kwargs...) 
    ϵdif=(ϵmax-ϵmin)/numpoints
    dosarray=[]
    for i in 0:numpoints
        ϵ=ϵmin+ϵdif*i
        push!(dosarray, density_of_states_wannier_quad(wannier_file, cell_map_file, ϵ; δ, kwargs... ))
    end
    return sum(ϵdif*dosarray)
end

function density_of_states_wannier_quad_check(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵmin::Real, ϵmax::Real, numpoints::Integer; δ::Real = 0.1, kwargs...) 
    ϵdif=(ϵmax-ϵmin)/numpoints
    dosarray=[]
    for i in 0:numpoints
        ϵ=ϵmin+ϵdif*i
        push!(dosarray, density_of_states_wannier_quad(HWannier, cell_map, ϵ; δ, kwargs... ))
    end
    return sum(ϵdif*dosarray)
end

function density_of_states_wannier(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}; mesh::Integer = 100,
    histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    WannierDOS=np.zeros(histogram_width*energy_range)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵ = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
        WannierDOS[round(Int, histogram_width*(ϵ+offset))] +=histogram_width*(1/mesh)^2
    end
    return WannierDOS
end


"""
$(TYPEDSIGNATURES)
"""
function wannierbandsoverlayedDOS(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ::Val{2}; 
    kpointsfile::AbstractString="bandstruct.kpoints", kticksfile::AbstractString="bandstruct.kpoints.in", 
    mesh::Integer = 100, histogram_width::Real = 100, spin::Integer=1, band_subplot::Vector{<:Int}=[1, 2, 1], 
    dos_subplot::Vector{<:Int}=[1, 2, 2], kwargs...)

    kpointlist = np.loadtxt(kpointsfile, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    energiesatkpoints = Vector{Float64}()
    for k in 1:num_kpoints
        push!(energiesatkpoints, wannier_bands(HWannier, cell_map, kpointlist[k, :]))
    end
    WannierDOSGather = Float64[]
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵ = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
        push!(WannierDOSGather, ϵ)
    end
    offset = minimum(WannierDOSGather) - 0.2
    energy_range = maximum(WannierDOSGather)-minimum(WannierDOSGather) + 0.4
    WannierDOS = np.zeros(round(Int, histogram_width*energy_range))
    for ϵ in WannierDOSGather
        WannierDOS[round(Int, histogram_width*(ϵ-offset))] += spin*histogram_width*(1/mesh)^2
    end
    subplot(band_subplot...)
    label_plots(kticksfile, kpointsfile)
    ylabel("Energy (eV)")
    plot(energiesatkpoints; kwargs...)
    ylim([offset, offset+energy_range])
    subplot(dos_subplot...)
    plot(WannierDOS, range(offset, offset+energy_range, length=length(WannierDOS)); kwargs...)
    ylim([offset, offset+energy_range])
    yticks(Float64[])
    xlabel("DOS (1/eV)")
    @assert sum(WannierDOS .* 1/histogram_width) ≈ spin "Error in Normalization of DOS"
end

"""
$(TYPEDSIGNATURES)
"""
function wannierbandsoverlayedDOS(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, kpoints::AbstractString, ::Val{3};
    mesh::Integer = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    kpointlist = np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    energiesatkpoints = Vector{Float64}()
    for k in 1:num_kpoints
        push!(energiesatkpoints, wannier_bands(HWannier, cell_map, kpointlist[k, :]))
    end
    WannierDOS = np.zeros(round(Int, histogram_width*energy_range))
    for (xmesh, ymesh, zmesh) in Tuple.(CartesianIndices(rand(mesh, mesh, mesh)))
        ϵ = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, zmesh/mesh])
        WannierDOS[round(Int, histogram_width*(ϵ+offset))] += histogram_width*(1/mesh)^3
    end
    plot(energiesatkpoints, ylims=[-offset, energy_range-offset])
    plot( WannierDOS, collect(1:histogram_width*energy_range))
end

"""
$(TYPEDSIGNATURES)
"""
function wannierbandsoverlayedDOS(HWannierUp::Array{Float64, 3}, cell_mapUp::Array{Float64, 2}, HWannierDn::Array{Float64, 3}, 
    cell_mapDn::Array{Float64, 2}; kpointsfile::AbstractString="bandstruct.kpoints", kticksfile="bandstruct.kpoints.in", 
    mesh::Integer=100, histogram_width::Real=100, band_subplot::Vector{<:Int}=[1, 2, 1], dos_subplot::Vector{<:Int}=[1, 2, 2], 
    color_up="blue", color_dn="red", kwargs...)

    kpointlist = np.loadtxt(kpointsfile, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    energiesatkpointsUp = Vector{Float64}()
    energiesatkpointsDn= Vector{Float64}()
    for k in 1:num_kpoints
        push!(energiesatkpointsUp, wannier_bands(HWannierUp, cell_mapUp, kpointlist[k, :]))
        push!(energiesatkpointsDn, wannier_bands(HWannierDn, cell_mapDn, kpointlist[k, :]))
    end
    WannierDOSGatherUp = Float64[]
    WannierDOSGatherDn = Float64[]
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵ = wannier_bands(HWannierUp, cell_mapUp, [xmesh/mesh, ymesh/mesh, 0])
        push!(WannierDOSGatherUp, ϵ)
        ϵ = wannier_bands(HWannierDn, cell_mapDn, [xmesh/mesh, ymesh/mesh, 0])
        push!(WannierDOSGatherDn, ϵ)
    end
    offset = minimum([WannierDOSGatherDn..., WannierDOSGatherUp...]) - 0.2
    energy_range = maximum([WannierDOSGatherDn..., WannierDOSGatherUp...]) - minimum([WannierDOSGatherDn..., WannierDOSGatherUp...]) + 0.4

    WannierDOSUp = np.zeros(round(Int, histogram_width*energy_range))
    WannierDOSDn = np.zeros(round(Int, histogram_width*energy_range))
    for ϵ in WannierDOSGatherUp
        WannierDOSUp[round(Int, histogram_width*(ϵ-offset))] += histogram_width*(1/mesh)^2
    end
    for ϵ in WannierDOSGatherDn  
        WannierDOSDn[round(Int, histogram_width*(ϵ-offset))] += histogram_width*(1/mesh)^2
    end
    subplot(band_subplot...)
    plot(energiesatkpointsUp; color=color_up, kwargs...)
    plot(energiesatkpointsDn; color=color_dn, kwargs...)
    ylabel("Energy (eV)")
    ylim([offset, offset+energy_range])
    label_plots(kticksfile, kpointsfile)
    subplot(dos_subplot...)
    plot(WannierDOSUp, range(offset, offset+energy_range, length=length(WannierDOSUp)), color=color_up; kwargs...)
    plot(WannierDOSDn, range(offset, offset+energy_range, length=length(WannierDOSDn)), color=color_dn; kwargs...)
    yticks(Float64[])
    xlabel("DOS(1/eV)")
    @assert sum(WannierDOSUp .* 1/histogram_width) ≈ 1 "Error in Normalization of DOS"
    @assert sum(WannierDOSDn .* 1/histogram_width) ≈ 1 "Error in Normalization of DOS"
    ylim([offset, offset+energy_range])
end

"""
$(TYPEDSIGNATURES)
"""
function wannierbandsoverlayedDOS(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer; 
    kpoints::AbstractString = "bandstruct.kpoints", mesh::Integer = 100, histogram_width::Real = 100, spin::Integer=1, kwargs...)

    kpointlist = np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    energiesatkpoints = Array{Float64, 2}(undef, (num_kpoints, nbands))
    for k in 1:num_kpoints
        energiesatkpoints[k, :] = wannier_bands(HWannier, cell_map, kpointlist[k, :], nbands)
    end
    WannierDOSGather = Float64[]
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵs = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0], nbands)
        for ϵ in ϵs
            push!(WannierDOSGather, ϵ)
        end
    end
    offset = minimum(WannierDOSGather) - 1
    energy_range = maximum(WannierDOSGather)-minimum(WannierDOSGather) + 3
    numdospoints = round(Int, histogram_width*energy_range)+1
    Energies = collect(energy_range/numdospoints:energy_range/numdospoints:energy_range) .+ offset
    WannierDOS = np.zeros(numdospoints)
    for ϵ in WannierDOSGather
        WannierDOS[round(Int, histogram_width*(ϵ-offset))+1] += spin*histogram_width*(1/mesh)^2
    end
    println(sum(WannierDOS.*1/histogram_width) )
    @assert sum(WannierDOS.*1/histogram_width) ≈ nbands*spin #Check normalization of density of states per unit cell.
    A = plot(energiesatkpoints, ylims=[offset, energy_range+offset], xticks = false, legend=false, ylabel = "Energy (eV)"; kwargs...)
    B = plot(WannierDOS, Energies, legend=false, xlabel = "DOS (1/eV)", yticks = false; kwargs...)
    display(plot(A, B, size=(1000, 500)))
end

"""
$(TYPEDSIGNATURES)
"""
function wannierbandsoverlayedDOS(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, kpoints::AbstractString, nbands::Integer, ::Val{3}; 
    mesh::Int = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0, kwargs...)

    kpoints_list = bandstructkpoints2q(filename=kpoints)
    energiesatkpoints = Array{Float64, 2}(undef, (length(kpoints_list), nbands))
    for (k, kpoint) in enumerate(kpoints_list)
        energiesatkpoints[k, :] = wannier_bands(HWannier, cell_map, kpoint, nbands)
    end
    WannierDOS = np.zeros(round(Int, histogram_width*energy_range))
    for _ in 1:mesh^3
        ϵs = wannier_bands(HWannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            WannierDOS[round(Int, histogram_width*(ϵ+offset))] += histogram_width*(1/mesh)^3
        end
    end
    figure()
    subplot(1, 2, 1)
    plot(energiesatkpoints, ylims=[-offset, energy_range-offset])
    ylabel("Energy (eV)")
    subplot(1, 2, 2)
    plot(WannierDOS, collect(1:histogram_width*energy_range)./histogram_width .-offset)
    xlabel(xlabel = "DOS (1/eV)")
end

"""
$(TYPEDSIGNATURES)
"""
function bandsoverlayedwannierDOS(band_file::AbstractString, dosfile::AbstractString, spin::Integer, ntotalbands::Integer, HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    kpoints::AbstractString, nbands::Integer, ::Val{3}; mesh::Integer = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0, kwargs...)
    kpointlist = np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    reshaped_exactenergies = reshape(read!(band_file, Array{Float64}(undef, ntotalbands*num_kpoints )),(ntotalbands, num_kpoints));
    exactenergies = permutedims(reshaped_exactenergies, [2, 1])*1/eV;
    energiesatkpoints = Array{Float64, 2}(undef, (num_kpoints, nbands))
    for k in 1:num_kpoints
        energiesatkpoints[k, :] = wannier_bands(HWannier, cell_map, kpointlist[k, :], nbands)
    end
    WannierDOS = np.zeros(round(Int, histogram_width*energy_range))
    for _ in 1:mesh^3
        ϵs = wannier_bands(HWannier, cell_map, rand(3), nbands)
        for ϵ in ϵs
            try
                WannierDOS[round(Int, histogram_width*(ϵ+offset))] += histogram_width*(1/mesh)^3
            catch
            end
        end
    end
    scfdosdata = try 
        np.loadtxt(dosfile)
    catch 
        np.loadtxt(dosfile, skiprows=1)
    end
    A = plot(energiesatkpoints, ylims=[-offset, energy_range-offset], xticks = false, legend=false, ylabel = "Energy (eV)", color="blue", linestyle=:dashdot, linewidth=5)
    A = plot!(exactenergies, ylims=[-offset, energy_range-offset], linewidth=3, color="red")
    B = plot(WannierDOS*spin, collect(1:histogram_width*energy_range)./histogram_width .-offset, legend=false, xlabel = "DOS (1/eV)", yticks = false, linestyle=:dashdot, linewidth=2)
    B = plot!(scfdosdata[:, 2]*eV, scfdosdata[:, 1]*1/eV, ylims=[-offset, energy_range-offset], xlims=[0, 1], linewidth=3)
    plot(A, B, size=(1000, 500); kwargs...)
end

function wannierbandsoverlayedDOS(HWannierUp::Array{Float64, 3}, cell_mapUp::Array{Float64, 2}, HWannierDn::Array{Float64, 3}, cell_mapDn::Array{Float64, 2},
    kpoints::AbstractString, nbands::Integer; mesh::Integer = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0, kwargs...)

    kpointlist = np.loadtxt(kpoints, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    energiesatkpointsUp = Array{Float64, 2}(undef, (num_kpoints, nbands))
    energiesatkpointsDn = Array{Float64, 2}(undef, (num_kpoints, nbands))
    for k in 1:num_kpoints
        energiesatkpointsUp[k, :] = wannier_bands(HWannierUp, cell_mapUp, kpointlist[k, :], nbands)
        energiesatkpointsDn[k, :] = wannier_bands(HWannierDn, cell_mapDn, kpointlist[k, :], nbands)
    end
    WannierDOSUp = np.zeros(round(Int, histogram_width*energy_range))
    WannierDOSDn = np.zeros(round(Int, histogram_width*energy_range))
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵups = wannier_bands(HWannierUp, cell_mapUp, [xmesh/mesh, ymesh/mesh, 0], nbands)
        ϵdns = wannier_bands(HWannierDn, cell_mapDn, [xmesh/mesh, ymesh/mesh, 0], nbands)
        for ϵup in ϵups
            WannierDOSUp[round(Int, histogram_width*(ϵup+offset))] += histogram_width*(1/mesh)^2
        end
        for ϵdn in ϵdns
            WannierDOSDn[round(Int, histogram_width*(ϵdn+offset))] += histogram_width*(1/mesh)^2
        end
    end
    A = plot(energiesatkpointsUp, ylims=[-offset, energy_range-offset], xticks = false, legend=false, ylabel = "Energy (eV)", linewidth=5; kwargs...)
    Aprim = plot!(energiesatkpointsDn, ylims=[-offset, energy_range-offset], xticks = false, legend=false, ylabel = "Energy (eV)", linewidth=5; kwargs...)
    B = plot(WannierDOSUp, collect(1:histogram_width*energy_range)./histogram_width .- offset, legend=false, xlabel = "DOS (1/eV)", yticks = false, linewidth=5; kwargs...)
    C = plot!(WannierDOSDn, collect(1:histogram_width*energy_range)./histogram_width .- offset,legend=false, xlabel = "DOS (1/eV)", yticks = false, linewidth=5; kwargs...)
    plot(A, C, size=(1000, 500))
end

"""
$(TYPEDSIGNATURES)
Returns an array of k points in the basis of reciprocal lattice vectors, with optional interpolation given by 
keyword argument interpolate. The k points are by default read from a file in the current directory given by 
bandstruct.kpoints, in keeping with JDFTX conventions. This can be changed by passing the keyword argument filename.
"""
function bandstructkpoints2q(;kpointsfile::AbstractString="bandstruct.kpoints", interpolate::Integer=1)
    kpointlist = np.loadtxt(kpointsfile, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    kpointsreshaped = Vector{Vector{Float64}}()
    for k in 1:num_kpoints-1
        for interpolatedk in 0:interpolate-1
            push!(kpointsreshaped, kpointlist[k, :] .+ (kpointlist[k+1, :].-kpointlist[k, :])./interpolate.*interpolatedk)
        end
    end
    push!(kpointsreshaped, kpointlist[num_kpoints, :])
    return kpointsreshaped
end

"""
$(TYPEDSIGNATURES)
The standard DOS function using wannier functions returns the density of states per eV per unit cell. 
At times it is more convenient to obtain the DOS per eV per angstrom^2 
"""
function density_of_states_wannier_per_area(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}; 
    mesh::Integer = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    ucell_area = unit_cell_area(lattice_vectors)
    WannierDOS=np.zeros(histogram_width*energy_range)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵ = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
        WannierDOS[round(Int, histogram_width*(ϵ+offset))] += histogram_width*(1/mesh)^2
    end
    return WannierDOS/ucell_area
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states_wannier(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, ::Val{3};
    exclude_bands::Vector{<:Integer} = Int[], mesh::Integer = 100, histogram_width::Integer = 100, energy_range::Real = 10, offset::Real = 0)
    !isempty(exclude_bands) && (0 ∈ exclude_bands && exclude_bands .+= 1) #Check for accidental 0 based indexing
    WannierDOS=np.zeros(histogram_width*energy_range)
    energies = collect(0:1/histogram_width:energy_range-1/histogram_width) .- offset
    for (xmesh, ymesh, zmesh) in Tuple.(CartesianIndices(rand(mesh, mesh, mesh)))
        ϵs =  wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, zmesh/mesh], nbands)
        for (band, ϵ) in enumerate(ϵs)
            (band ∈ exclude_bands) && continue
            WannierDOS[round(Int, histogram_width*(ϵ+offset))] += histogram_width*(1/mesh)^3
        end
    end
    @assert sum(WannierDOS.*1/histogram_width) ≈ nbands - length(exclude_bands) #Verify normalization 
    return energies, WannierDOS
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states_wannier(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, ::Val{2};
    exclude_bands::Vector{<:Integer} = Int[], mesh::Integer = 100, histogram_width::Integer = 100, energy_range::Real = 10, offset::Real = 0)
    !isempty(exclude_bands) && (0 ∈ exclude_bands && exclude_bands .+= 1) #Check for accidental 0 based indexing
    WannierDOS=np.zeros(histogram_width*energy_range)
    energies = collect(0:1/histogram_width:energy_range-1/histogram_width) .- offset
    for (xmesh, ymesh, zmesh) in Tuple.(CartesianIndices(rand(mesh, mesh, mesh)))
        ϵs =  wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, zmesh/mesh], nbands)
        for (band, ϵ) in enumerate(ϵs)
            (band ∈ exclude_bands) && continue
            WannierDOS[round(Int, histogram_width*(ϵ+offset))] += histogram_width*(1/mesh)^3
        end
    end
    @assert sum(WannierDOS.*1/histogram_width) ≈ nbands - length(exclude_bands) #Verify normalization 
    return energies, WannierDOS
end

density_of_states_wannier(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer; kwargs...) = 
density_of_states_wannier(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, Val(2); kwargs...)

function density_of_states_montecarlo(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, ::Val{2}; exclude_bands = Int[], mesh::Integer = 100, histogram_width::Integer = 100, energy_range::Real = 10, offset::Real = 0)
    WannierDOS=np.zeros(histogram_width*energy_range)
    randks = rand(mesh^2, 2)
    energies = collect(0:1/histogram_width:energy_range-1/histogram_width) .- offset
    for randk in eachrow(randks)
        ϵs = wannier_bands(HWannier, cell_map, [collect(randk)..., 0], nbands)
        for (ϵ, band) in enumerate(ϵs)
            (band in exclude_bands) && continue
            WannierDOS[round(Int, histogram_width*(ϵ+offset))] += histogram_width*(1/mesh)^2
        end
    end
    @assert sum(WannierDOS.*1/histogram_width) ≈ nbands - length(exclude_bands) #Verify normalization 
    return energies, WannierDOS
end

function density_of_states_montecarlo(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2},
    nbands::Integer, ::Val{3}; exclude_bands::Vector{<:Integer} = Integer[], mesh::Integer = 100, 
    histogram_width::Integer = 100, energy_range::Real = 10, offset::Real = 0)
    energies = collect(0:1/histogram_width:energy_range-1/histogram_width) .- offset
    WannierDOS=np.zeros(histogram_width*energy_range)
    randks = rand(mesh^3, 3 )
    for randk in eachrow(randks)
        ϵ =  wannier_bands(HWannier, cell_map, collect(randk), nbands)
        for band in 1:nbands
            (band in exclude_bands) && continue
            ϵ_band = ϵ[band]
            WannierDOS[round(Int, histogram_width*(ϵ_band+offset))] += histogram_width*(1/mesh)^3
        end
    end
    @assert sum(WannierDOS.*1/histogram_width) ≈ nbands - length(exclude_bands) #Verify normalization 
    return energies, WannierDOS
end

density_of_states_montecarlo(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer; kwargs...) = 
density_of_states_montecarlo(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, Val(2); kwargs...)

"""
$(TYPEDSIGNATURES)

### Keyword Arguments

`plotoccupations` : Whether to plot filling as a function of energy with vertical lines indicating 1/4, 1/2, 3/4 fillings

`histogram_width` : The histogram binning that is used to compute the density of states
"""
function find_chemical_potential(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}; 
    mesh::Integer=100, histogram_width::Integer=100, energy_range::Real=10, offset::Real=0, plotoccupations::Bool=true)

    doss = density_of_states_wannier(HWannier, cell_map, mesh=mesh, histogram_width=histogram_width, energy_range=energy_range, offset=offset )
    totalstates = []
    for i in 1:length(doss)
        push!(totalstates, [i/histogram_width-offset, sum(doss[1:i]*1/histogram_width)])
    end
    xenergies = []
    yoccupations = []
    for i in 1:length(doss)
        push!(xenergies, totalstates[i][1])
        push!(yoccupations, totalstates[i][2])
    end
    onequarter  = xenergies[argmin(abs.(yoccupations .- 0.25))]
    halffilling = xenergies[argmin(abs.(yoccupations .- 0.5))]
    threequarter = xenergies[argmin(abs.(yoccupations .- 0.75))]
    println("Quarter Filling: ", onequarter)
    println("Half Filling: ", halffilling)
    println("Three Quarters Filling: ", threequarter)
    plotoccupations && plot(xenergies, yoccupations, linewidth=5)
    plotoccupations && vlines([onequarter, halffilling, threequarter], -1, 2*maximum(yoccupations), linewidth=5)
    return xenergies, yoccupations
end

function find_chemical_potential(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer;
    mesh::Real = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)
    doss = density_of_states_wannier(HWannier, cell_map, nbands, mesh=mesh, histogram_width=histogram_width, energy_range=energy_range, offset=offset )
    totalstates = []
    for i in 1:length(doss[2])
        push!(totalstates, [i/histogram_width-offset, sum(doss[2][1:i]*1/histogram_width)])
    end
    xenergies = []
    yoccupations = []
    for i in 1:length(doss[2])
        push!(xenergies, totalstates[i][1])
        push!(yoccupations, totalstates[i][2])
    end
    return xenergies, yoccupations
end

"""
Returns an array of occupations. Method to find chemical potential at finite temperature. 
"""
function finite_temperature_chemical_potential(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, T::Real;
    mesh::Real = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    doss = density_of_states_wannier(HWannier, cell_map, mesh=mesh, histogram_width=histogram_width, energy_range=energy_range, offset=offset )
    occupations_array = Float64[]
    for i in 1:length(doss)
        μ = i/histogram_width-offset
        Fermi = x -> 1/(exp((x-μ)/(kB*T))+1)
        Occupations= Fermi.((1/histogram_width-offset):1/histogram_width:(length(doss)/histogram_width-offset))
        push!(occupations_array, sum(Occupations.*doss)*1/histogram_width)
    end
    return collect((1/histogram_width-offset):1/histogram_width:(length(doss)/histogram_width-offset)), occupations_array
end

function find_num_phonons(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}; mesh::Integer = 100, histogram_width::Integer = 100, energy_range::Real = 2)
    doss = phonon_density_of_states(force_matrix, phonon_cell_map; mesh=mesh, histogram_width=histogram_width, energy_range=energy_range)
    totalstates = []
    for i in 1:length(doss)
        push!(totalstates, [i/histogram_width, sum(doss[1:i]*1/histogram_width)])
    end
    xenergies = []
    yoccupations = []
    for i in 1:length(doss)
        push!(xenergies, totalstates[i][1])
        push!(yoccupations, totalstates[i][2])
    end
    return xenergies, yoccupations
end


"""
$(TYPEDSIGNATURES)
"""
function phonon_density_of_states(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, ::Val{2}; mesh::Integer = 100, 
    histogram_width::Integer = 100, energy_range::Real = 2)
    PhononDOS=np.zeros(histogram_width*energy_range)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ωs=  phonon_dispersion(force_matrix, phonon_cell_map, [xmesh/mesh, ymesh/mesh, 0])
        for ω in ωs
            ω < 0 && continue
            PhononDOS[round(Int, histogram_width*ω)+1] += histogram_width*(1/mesh)^2
        end
    end
    return PhononDOS
end

"""
$(TYPEDSIGNATURES)
"""
function phonon_density_of_states(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, ::Val{3}; mesh::Integer = 100, 
    histogram_width::Integer = 100, energy_range::Real = 2)
    PhononDOS=np.zeros(histogram_width*energy_range)
    for (xmesh, ymesh, zmesh) in Tuple.(CartesianIndices(rand(mesh, mesh, mesh)))
        ωs=  phonon_dispersion(force_matrix, phonon_cell_map, [xmesh/mesh, ymesh/mesh, zmesh/mesh])
        for ω in ωs
            ω < 0 && continue
            PhononDOS[round(Int, histogram_width*ω)+1] += histogram_width*(1/mesh)^3
        end
    end
    return PhononDOS
end

phonon_density_of_states(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}; mesh::Integer = 100, histogram_width::Integer = 100, energy_range::Real = 2) = phonon_density_of_states(force_matrix, phonon_cell_map, Val(2); mesh= mesh, histogram_width= histogram_width, energy_range = energy_range)
    

function phononbandsoverlayedDOS(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}; kpointsfile::AbstractString="bandstruct.kpoints",
    kticksfile::AbstractString="bandstruct.kpoints.in", mesh::Integer=100, histogram_width::Integer=100, 
    band_subplot::Vector{<:Int}=[1, 2, 1], dos_subplot::Vector{<:Int}=[1, 2, 2], return_tot::Bool=false, kwargs...)

    kpointlist = np.loadtxt(kpointsfile, skiprows=2, usecols=[1, 2, 3])
    num_kpoints = np.shape(kpointlist)[1]
    PhononBands = Array{Float64, 2}( undef, (length(phonon_dispersion(force_matrix, phonon_cell_map, [0, 0, 0])), num_kpoints))
    for k in 1:num_kpoints
        PhononBands[:, k] = phonon_dispersion(force_matrix, phonon_cell_map, kpointlist[k, :])
    end
    
    PhononDOSGather = Float64[]

    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ωs =  phonon_dispersion(force_matrix, phonon_cell_map, [xmesh/mesh, ymesh/mesh, 0])
        for ω in ωs
            ω < 0 && continue
            push!(PhononDOSGather, ω)
        end
    end

    energy_range = maximum(PhononDOSGather) + 2/histogram_width
    PhononDOS = np.zeros(round(Int, histogram_width*energy_range))

    for ω in PhononDOSGather
        PhononDOS[round(Int, histogram_width*ω)+1] += histogram_width*(1/mesh)^2
    end

    subplot(band_subplot...)
    plot(transpose(PhononBands); kwargs...)
    ylim(0, energy_range)
    ylabel("Energy (eV)")
    label_plots(kticksfile, kpointsfile)
    subplot(dos_subplot...)
    plot(PhononDOS, range(0, energy_range, length=length(PhononDOS)); kwargs...)
    ylim([0, energy_range])
    yticks(Float64[])
    xlabel("DOS (1/eV)")
    return_tot && println(sum(PhononDOS)/histogram_width)
end

"""
$(TYPEDSIGNATURES)
"""
function phonon_density_of_states_per_area(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, lattice_vecs::Vector{<:Vector{<:Real}}; 
    mesh::Integer = 100, histogram_width::Integer = 100, energy_range::Real = 2)

    PhononDOS = np.zeros(histogram_width*energy_range)
    ucell_area = unit_cell_area(lattice_vecs)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ωs=  phonon_dispersion(force_matrix, phonon_cell_map, [xmesh/mesh, ymesh/mesh, 0])
        for ω in ωs
            ω < 0 && continue
            PhononDOS[round(Int, histogram_width*ω)+1]=PhononDOS[round(Int, histogram_width*ω)+1]+histogram_width*(1/mesh)^2
        end
    end
    return PhononDOS/ucell_area
end

