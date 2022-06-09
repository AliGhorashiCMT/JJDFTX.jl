#========================== First we calculate electronic heat capacities ==========================#
"""
$(TYPEDSIGNATURES)
"""
function fermiderivative(ϵ::Real, μ::Real, T::Real)
    T < 0 && error("Temperature must be a positive real number") 
    β = 1/(kB*T)
    x = exp(β*(ϵ-μ))
    fermid = ((ϵ-μ)/(kB*T^2))*(x/(x+1)^2)
    isnan(fermid) && println("NaN value- increase temperature")
    return fermid 
end

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
"""
function electron_heatcapacity(lattice_vectors::Vector{<:Vector{<:Real}}, Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, 
    μ::Real, Ts::Vector{<:Real}, ::Val{D}=Val(2); kwargs...) where D
    V = D == 2 ? unit_cell_area(lattice_vectors) : unit_cell_volume(lattice_vectors)
    Energies, DOS = density_of_states(Hwannier, cell_map, Val(D), monte_carlo=true; kwargs...)   
    return [sum(diff(Energies).*((Energies .- μ).*DOS.*fermiderivative.(Energies, μ, T)*1/V)[1:end-1]) for T in Ts]
end

"""
$(TYPEDSIGNATURES)
Returns the heat capacity in eV/(Kelvin*angstrom^3)
"""
function lattice_heatcapacity(lattice_vectors::Vector{<:Vector{<:Real}}, force_matrix::Array{<:Real, 3}, cellph_map::Array{<:Real, 2}, 
    Ts::Vector{<:Real}, ::Val{D}=Val(2); kwargs...) where D
    V = D == 2 ? unit_cell_area(lattice_vectors) : unit_cell_volume(lattice_vectors)
    Energies, DOS = phonon_dos(force_matrix, cellph_map, Val(D), monte_carlo=true; kwargs...)
    return [sum(diff(Energies).*(Energies.*DOS.*bosederivative.(Energies, T)*1/V)[1:end-1]) for T in Ts]
end

"""
$(TYPEDSIGNATURES)
An analytic approximation 
"""
function diraccone_heatcapacity(μ::Real, T::Real)
    #Factor of 36 comes from denominator of fermi velocity squared, factor of 4 comes from derivative of T^2, spin, valley degeneracies and dividing by (2pi)^2 in the denominator
    return (π^2/6)*4*kB^2*T*μ/(5.75^2*π) #Units of ev^2/Kelvin*1/(eV^2*angstrom^2)
end

"""
$(TYPEDSIGNATURES)
"""
function graphene_phononheatcapacity(T::Real)
    #Take sound velocity to be 20 km/s -> 2*10^14 angstroms/second. Therefore, dispersion is given by 6.6*2*10^(-2)*k 
    # Take a cutoff wavevector of 1 angstroms^(-1). 
    prefactor = (1/(1.3*6.6e-2)^2)*(2/(4*π^2))+(1/(2*6.6e-2)^2)*(1/(4*π^2)) # Different logitudinal and transverse sound velocities. 
    prefactor*sum(bosederivative.(collect(0.0001:0.001:0.16), Ref(T)).*collect(0.0001:0.001:0.16))*0.001
end

