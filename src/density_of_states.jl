function make_mesh(mesh::Integer, ::Val{2})
    x = y = collect(range(0, 1-1/mesh, length=mesh))

    x = np.repeat(x, mesh)
    y = repeat(y, mesh)
    z = zeros(mesh^2)
    return hcat(x, y, z)
end

function make_mesh(mesh::Integer, ::Val{3})
    x = y = z = collect(range(0, 1-1/mesh, length=mesh))

    x = np.repeat(x, mesh^2)
    y = repeat(np.repeat(y, mesh), mesh)
    z = repeat(z, mesh^2)

    return hcat(x, y, z)
end


"""
Returns the density of states data as outputed by dump End DOS. We convert the energies dimension from Hartree to eV and 
we convert the DOS dimension from 1/Hartree to 1/eV. 
"""
function load_dos_data(dosfile::AbstractString)
    dosdata = try 
        np.loadtxt(dosfile)
    catch 
        np.loadtxt(dosfile, skiprows=1)
    end
    return dosdata[:, 1]*1/eV, dosdata[:, 2]*eV
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states(dosfile_up::AbstractString, dosfile_dn::AbstractString; color_up::AbstractString="blue", 
    color_dn::AbstractString="red", return_tot::Bool=false, kwargs... )

    energies_up, dos_up = load_dos_data(dosfile_up)
    energies_dn, dos_dn = load_dos_data(dosfile_dn)

    plot(energies_up, dos_up, linewidth=4, color=color_up, label="Spin Up"; kwargs...)
    plot(energies_dn, dos_dn, linewidth=4, color=color_dn, label="Spin Down"; kwargs...)
    return_tot && println("Total number of spin up electrons is: ", sum(diff(energies_up).*dos_up[2:end]))
    return_tot && println("Total number of spin down electrons is: ", sum(diff(energies_dn).*dos_dn[2:end]))

    ylabel("Density of States (1/eV/Cell)")
    xlabel("Energy (eV)")
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states(dosfile::AbstractString; return_tot::Bool=false, kwargs...)
    energies, dos = load_dos_data(dosfile)
    plot(energies, dos, linewidth=4; kwargs...)
    return_tot && println("Total number of electrons is: ", sum(diff(energies).*dos[2:end]))
    ylabel("Density of States (1/eV/Cell)")
    xlabel("Energy (eV)")
end

function find_chemical_potential(energies::Vector{<:Real}, dos::Vector{<:Real})
    return energies[2:end], [sum(dos[1:idx-1] .* diff(energies[1:idx])) for idx in eachindex(energies)[2:end]]
end

function find_chemical_potential(dosfile::AbstractString)
    energies, dos = load_dos_data(dosfile)
    find_chemical_potential(energies, dos)
end

"""
$(TYPEDSIGNATURES)
Returns relevant properties of the density of states. Since this is basically only bandgaps, this function returns the energies at which 
the density of states vanishes.
"""
function dos_properties(dosfile::AbstractString)
    energies, dos = load_dos_data(dosfile)
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
Overlay the bandstructure with the density of states.
"""
function bands_overlayed_dos(dos_info::Tuple{Vector{<:Float64}, Vector{<:Real}}, energies::Array{Float64, 2}; 
    energy_range::Tuple{<:Real, <:Real} = (-10, 10), kpointsfile::AbstractString="bandstruct.kpoints", 
    kticksfile::AbstractString = "bandstruct.kpoints.in", label_plot::Bool=true, band_subplot::Vector{<:Integer}=[1, 2, 1],
    dos_subplot::Vector{<:Integer}=[1, 2, 2], dos_yticks::Bool=false, return_tot::Bool=false, kwargs...)
    
    dos_energies, dos = dos_info

    subplot(band_subplot...)
    plot(energies, color="black", linewidth=2; kwargs...) 
    ylim(collect(energy_range))
    ylabel("Energy (eV)")
    
    label_plot && label_plots(kticksfile, kpointsfile)    
    lowerDOS = argmin(abs.(dos_energies .- energy_range[1]))
    upperDOS = argmin(abs.(dos_energies .- energy_range[2]))
    
    subplot(dos_subplot...)
    plot(dos, dos_energies, color="black", linewidth=2; kwargs...)
    ylim(collect(energy_range))
    xlim([0, maximum(dos[lowerDOS:upperDOS])])
    xlabel("DOS (1/eV)")
    !dos_yticks && yticks(Float64[])

    return_tot && println("Total number of electrons in range: ", 
    sum(diff(dos_energies)[lowerDOS:upperDOS] .* dos[lowerDOS:upperDOS]))
end

function bands_overlayed_dos(dos_info_up::Tuple{Vector{<:Float64}, Vector{<:Real}}, dos_info_dn::Tuple{Vector{<:Float64}, Vector{<:Real}}, 
    energies_up::Array{<:Float64, 2}, energies_dn::Array{<:Float64, 2}; energy_range::Tuple{<:Real, <:Real}=(-10, 10),
    color_up="blue", color_dn="red", label_plot::Bool=true, kticksfile::AbstractString="bandstruct.kpoints.in", kpointsfile::AbstractString="bandstruct.kpoints", 
    return_tot::Bool=false, band_subplot::Vector{<:Integer}=[1, 2, 1], dos_subplot::Vector{<:Integer}=[1, 2, 2], dos_yticks::Bool=true, kwargs...)

    dos_energies_up, dos_up = dos_info_up
    dos_energies_dn, dos_dn = dos_info_dn

    subplot(band_subplot...)
    plot(energies_up, color=color_up, label="", linewidth=2; kwargs...)
    ylim(collect(energy_range))
    plot(energies_dn, color=color_dn, label="", linewidth=2; kwargs...)
    ylabel("Energy (eV)")

    label_plot && label_plots(kticksfile, kpointsfile)

    lowerDOS_up = argmin(abs.(dos_energies_up .- energy_range[1]))
    upperDOS_up = argmin(abs.(dos_energies_up .- energy_range[2]))
    lowerDOS_dn = argmin(abs.(dos_energies_dn .- energy_range[1]))
    upperDOS_dn = argmin(abs.(dos_energies_dn .- energy_range[2]))
    max1 = maximum((dos_up)[lowerDOS_up:upperDOS_up])
    max2 = maximum((dos_dn)[lowerDOS_dn:upperDOS_dn])
    max = maximum([max1, max2])

    subplot(dos_subplot...)
    plot(dos_up, dos_energies_up, linewidth=2, color=color_up; kwargs...)
    plot(dos_dn, dos_energies_dn, linewidth=2, color=color_dn; kwargs...)
    !dos_yticks && yticks(Float64[])

    return_tot && println("Total number of electrons in range: ", 
    sum(diff(dos_energies_up)[lowerDOS_up:upperDOS_up] .* (dos_up)[lowerDOS_up:upperDOS_up]) + 
    sum(diff(dos_energies_dn)[lowerDOS_dn:upperDOS_dn] .* (dos_dn)[lowerDOS_dn:upperDOS_dn]))

    xlabel("DOS (1/eV)")
    xlim(0, max)
    ylim(collect(energy_range))
end

function bands_overlayed_dos(dosfile::AbstractString, band_file::AbstractString; kpointsfile::AbstractString="bandstruct.kpoints", 
    kwargs...)

    numpoints, numbands = load_bands_points(band_file, kpointsfile, 1)
    energies = load_bandeigs_data(band_file, numpoints, numbands, 1)
    dos_info = load_dos_data(dosfile)
    dos_energies, dos = dos_info
    @assert isapprox(numbands*2, sum(diff(dos_energies).*dos[2:end]), atol=1e-1) "DOS not propertly normalized. Make sure files are correct"

    bands_overlayed_dos(dos_info, energies; kpointsfile, kwargs...)
end

function bands_overlayed_dos(dosfile_up::AbstractString, dosfile_dn::AbstractString, band_file::AbstractString; 
    kpointsfile::AbstractString="bandstruct.kpoints", kwargs...)

    numpoints, numbands = load_bands_points(band_file, kpointsfile, 2)
    energies = load_bandeigs_data(band_file, numpoints, numbands, 2)

    energies_up = energies[1:numpoints, :]
    energies_dn = energies[numpoints+1:end, :]

    dos_info_up = load_dos_data(dosfile_up)
    dos_info_dn = load_dos_data(dosfile_dn)

    dos_energies_up, dos_up = dos_info_up
    dos_energies_dn, dos_dn = dos_info_dn
    
    @assert isapprox(numbands, sum(diff(dos_energies_up).*dos_up[2:end]), atol=1e-1) "DOS not propertly normalized. Make sure files are correct"
    @assert isapprox(numbands, sum(diff(dos_energies_dn).*dos_dn[2:end]), atol=1e-1) "DOS not propertly normalized. Make sure files are correct"

    bands_overlayed_dos(dos_info_up, dos_info_dn, energies_up, energies_dn; kpointsfile, kwargs...)
end

function dos_cubature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵ::Real; δ::Real = 0.1, kwargs...)
    1/π*hcubature(vec->imag(-1/(ϵ-wannier_bands(HWannier, cell_map, [vec[1], vec[2], 0])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]
end

function dos_cubature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵs::Vector{<:Real}; δ::Real = 0.1, kwargs...)
    dos_vec = Float64[]
    for ϵ in ϵs
        push!(dos_vec, density_of_states_wannier_quad(HWannier, cell_map, ϵ, δ=δ;  kwargs...))
    end
    return dos_vec
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

function density_of_states_wannier_quad_check(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵmin::Real, ϵmax::Real, numpoints::Integer; δ::Real = 0.1, kwargs...) 
    ϵdif=(ϵmax-ϵmin)/numpoints
    dosarray=[]
    for i in 0:numpoints
        ϵ=ϵmin+ϵdif*i
        push!(dosarray, density_of_states_wannier_quad(HWannier, cell_map, ϵ; δ, kwargs... ))
    end
    return sum(ϵdif*dosarray)
end

function collect_dos(DOS_GATHER::Vector{<:Real}; histogram_width::Integer=10)
    offset = minimum(DOS_GATHER) - 1
    energy_range = maximum(DOS_GATHER) - minimum(DOS_GATHER) + 2

    DOS = zeros(round(Int, histogram_width*energy_range))
    energies = collect(range(offset, offset + energy_range, length=length(DOS))) 

    for ϵ in DOS_GATHER 
        DOS[round(Int, histogram_width*(ϵ-offset))] += histogram_width
    end
    
    return energies, DOS
end

"""
$(TYPEDSIGNATURES)
"""
function density_of_states(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, dim::Val{D} = Val(2);
    whichbands::Union{Vector{<:Integer}, Nothing} = nothing, monte_carlo::Bool = false, mesh::Integer = 100,
    num_blocks::Integer=10, histogram_width::Integer = 100) where D

    DOS_GATHER = Float64[]

    for _ in 1:num_blocks
        kpoints = !monte_carlo ? transpose(make_mesh(mesh, dim)) : vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        Es, _ = wannier_bands(Hwannier, cell_map, kpoints)
        Es = isnothing(whichbands) ? Es : Es[:, whichbands]
        DOS_GATHER = [DOS_GATHER..., vec(Es)...]
    end
    energies, dos = collect_dos(DOS_GATHER, histogram_width = histogram_width)
    dos *= (1/mesh^D)*(1/num_blocks)
    !isnothing(whichbands) ? (@assert sum(dos * 1/histogram_width) ≈ length(whichbands)) :
    (@assert sum(dos * 1/histogram_width) ≈ size(Hwannier)[2])
    return energies, dos
end

function find_chemical_potential(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2};
    mesh::Real = 100, histogram_width::Real = 100)
    energies, dos = density_of_states_wannier(HWannier, cell_map, mesh=mesh, histogram_width=histogram_width)
    find_chemical_potential(energies, dos)
end

function find_num_phonons(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}; kwargs...)
    energies, dos = phonon_dos(force_matrix, cellph_map; kwargs...)
    find_chemical_potential(energies, dos)    
end

"""
$(TYPEDSIGNATURES)
"""
function phonon_dos(force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}, dim::Val{D} = Val(2); mesh::Integer = 100,
    num_blocks::Integer = 10, histogram_width::Integer = 100, monte_carlo::Bool = false, kwargs...) where D
    
    DOS_GATHER = Float64[]
    for _ in 1:num_blocks
        kpoints = !monte_carlo ? transpose(make_mesh(mesh, dim)) : vcat(rand(D, mesh^D), zeros(3-D, mesh^D))
        ωs =  phonon_dispersion(force_matrix, cellph_map, kpoints; kwargs...)
        DOS_GATHER = [DOS_GATHER..., ωs...]
    end

    energies, dos = collect_dos(DOS_GATHER, histogram_width = histogram_width)
    dos *= (1/mesh^D)*(1/num_blocks)
    @assert isapprox(sum(dos*1/histogram_width), size(force_matrix)[2], atol=1e-3) "$(sum(dos * 1/histogram_width))"
    return energies, dos
end

