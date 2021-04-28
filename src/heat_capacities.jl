#========================== First we calculate electronic heat capacities ==========================#
function electron_heatcapacity()

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
