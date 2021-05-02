#========================== First we calculate electronic heat capacities ==========================#
"""
$(TYPEDSIGNATURES)
"""
function fermiderivative(ϵ::Real, μ::Real, T::Real)
    T < 0 && error("Temperature must be a positive real number") 
    β = 1/(kB*T)
    x = exp(β*(ϵ-μ))
    fermid = ((ϵ-μ)/(kB*T^2))*(x/(x+1)^2)
    isnan(fermid) ? println("NaN value- increase temperature") : nothing
    return fermid 
end

"""
$(TYPEDSIGNATURES)
"""
function electron_heatcapacity(μ::Real, T::Real, lat::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer; exclude_bands = Int[], mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)
    WannierDOS=np.zeros(histogram_width*energy_range)
    DOSweightedϵ=np.zeros(round(Int, histogram_width*energy_range))
    vol = unit_cell_volume(lat)
    for (xmesh, ymesh, zmesh) in Tuple.(CartesianIndices(rand(mesh, mesh, mesh)))
        ϵs=wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, zmesh/mesh], nbands)
        for (bandidx, ϵ) in enumerate(ϵs)
            bandidx ∈ exclude_bands && continue 
            WannierDOS[round(Int, histogram_width*(ϵ+offset))]=WannierDOS[round(Int, histogram_width*(ϵ+offset))]+histogram_width*(1/mesh)^3
            DOSweightedϵ[round(Int, histogram_width*(ϵ+offset))]=DOSweightedϵ[round(Int, histogram_width*(ϵ+offset))]+(ϵ-μ)*fermiderivative(ϵ, μ, T)*histogram_width*1/vol*(1/mesh)^3
        end
    end
    @assert sum(WannierDOS ./ histogram_width) ≈ nbands #Check normalization of DOS
    return sum(DOSweightedϵ  ./ histogram_width) 
end

"""
$(TYPEDSIGNATURES)
Heat capacity for 2d systems. 
"""
function electron_heatcapacity(μ::Real, T::Real, lat::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, ::Val{2}; exclude_bands = Int[], mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)
    WannierDOS=np.zeros(histogram_width*energy_range)
    DOSweightedϵ=np.zeros(round(Int, histogram_width*energy_range))
    vol = unit_cell_area(lat)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵs=wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0], nbands)
        for (bandidx, ϵ) in enumerate(ϵs)
            bandidx ∈ exclude_bands && continue 
            WannierDOS[round(Int, histogram_width*(ϵ+offset))]=WannierDOS[round(Int, histogram_width*(ϵ+offset))]+histogram_width*(1/mesh)^2
            DOSweightedϵ[round(Int, histogram_width*(ϵ+offset))]=DOSweightedϵ[round(Int, histogram_width*(ϵ+offset))]+(ϵ-μ)*fermiderivative(ϵ, μ, T)*histogram_width*1/vol*(1/mesh)^2
        end
    end
    @assert sum(WannierDOS ./ histogram_width) ≈ nbands - length(exclude_bands) #Check normalization of DOS
    return sum(DOSweightedϵ  ./ histogram_width) 
end

"""
$(TYPEDSIGNATURES)
An analytic approximation 
"""
function diraccone_heatcapacity(T::Real)
    return 4*kB^2*T/(36*π)
end

"""
$(TYPEDSIGNATURES)

"""
function electron_heatcapacities(μ::Real, Ts::Union{Vector{<:Real}, StepRange{<:Integer, <:Integer}, UnitRange{<:Integer}}, lat::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer; exclude_bands = Int[], mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)
    cs = Float64[]
    for T in Ts
        println("T: ", T)
        push!(cs, electron_heatcapacity(μ, T, lat, HWannier, cell_map, nbands; exclude_bands, mesh, histogram_width, energy_range, offset))
    end
    return cs
end

"""
$(TYPEDSIGNATURES)

"""
function electron_heatcapacities(μ::Real, Ts::Union{Vector{<:Real}, StepRange{<:Integer, <:Integer}, UnitRange{<:Integer}}, lat::Vector{<:Vector{<:Real}}, HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Integer, ::Val{2}; exclude_bands = Int[], mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)
    cs = Float64[]
    for T in Ts
        println("T: ", T)
        push!(cs, electron_heatcapacity(μ, T, lat, HWannier, cell_map, nbands, Val(2); exclude_bands, mesh, histogram_width, energy_range, offset))
    end
    return cs
end


#========================== Next we calculate the lattice heat capacity (corresponding to phonon contributions) ==========================#

#= We must integrate over energies from 0 to infinity. The integrand is ϵD(ϵ)nb'(ϵ)
where the derivative is with respect to temperature
=#

"""
$(TYPEDSIGNATURES)
"""
function bosederivative(ϵ::Real, T::Real)
    T < 0 && error("Temperature must be a positive real number") 
    β = 1/(kB*T)
    x = exp(β*ϵ)
    bosed = (ϵ/(kB*T^2))*(x/(x-1)^2)
    isnan(bosed) ? println("NaN value- increase temperature or decrease phonon energy") : nothing
    return bosed 
end

"""
$(TYPEDSIGNATURES)
Returns the heat capacity in eV/(Kelvin*angstrom^3)
"""
function lattice_heatcapacity(T::Real, lat::Vector{<:Vector{<:Real}}, forcematrix::Array{<:Real, 3}, phononcellmap::Array{<:Real, 2}; mesh::Integer=100, histogram_width::Real=1000, energy_range::Real=0.1)
    PhononDOS=np.zeros(round(Int, histogram_width*energy_range))
    DOSweightedϵ=np.zeros(round(Int, histogram_width*energy_range))
    vol = unit_cell_volume(lat)
    nmodes = length(phonon_dispersion(forcematrix, phononcellmap, [0, 0, 0]))
    for (xmesh, ymesh, zmesh) in Tuple.(CartesianIndices(rand(mesh, mesh, mesh)))
        ωs=  phonon_dispersion(forcematrix, phononcellmap, [xmesh/mesh, ymesh/mesh, zmesh/mesh])
        for ω in ωs
            if ω>0
                PhononDOS[round(Int, histogram_width*ω)+1]=PhononDOS[round(Int, histogram_width*ω)+1]+histogram_width*(1/mesh)^3
                DOSweightedϵ[round(Int, histogram_width*ω)+1]=DOSweightedϵ[round(Int, histogram_width*ω)+1]+ω*bosederivative(ω, T)*histogram_width*1/vol*(1/mesh)^3
            end
        end
    end
    @assert (sum(PhononDOS ./ histogram_width) ≈ nmodes) #Check Normalization
    return (sum(DOSweightedϵ ./ histogram_width))
end
