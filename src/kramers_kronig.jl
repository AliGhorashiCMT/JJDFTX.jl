"""
$(TYPEDSIGNATURES)
Applies the kramers-kronig relations onto a 1 dimensional array of numbers consisting of the imaginary value of the polarization to return the real value of polarization
"""
function kramers_kronig(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real}; δ::Real=0.01) 
    2/π*sum(diff(energies).*(polarizations.*energies./(energies.^2 .- (ω+ δ*1im)^2))[1:end-1])
end

function kramers_kronig_reverse(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real}; δ::Real=0.01) 
    -2/π*sum(diff(energies).*(polarizations.*ω./(energies.^2 .- (ω+ δ*1im)^2))[1:end-1])
end

function kramers_kronig(im_pol::Function, ω::Real, max_energy_integration::Real=10; min_energy_integration::Real=0, kwargs...)
    cauchy_inner_function(omegaprime) = 2/pi*im_pol(omegaprime)*omegaprime/(omegaprime+ω)
    ErrorAbs=1e-20
    return pyintegrate.quad(cauchy_inner_function, min_energy_integration, max_energy_integration, weight="cauchy", 
        epsrel=ErrorAbs, epsabs=ErrorAbs, limit=175,  wvar= ω ; kwargs...)[1]
end

function kramers_kronig_reverse(im_pol::Function, ω::Real, max_energy_integration=10; kwargs...) 
    cauchy_inner_function(omegaprime) = -2/pi*im_pol(omegaprime)*ω/(omegaprime+ω)
    ErrorAbs=1e-20
    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy", 
        epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75, wvar=ω; kwargs...)[1]
end

"""
$(TYPEDSIGNATURES)
"""
function kramers_kronig_reverse_scipy(ω::Real, re_pol::Vector{<:Real}, max_energy::Real, domega::Real, max_energy_integration::Real; kwargs...) 
    interpolated_res=interpol.interp1d(0:domega:max_energy, re_pol)
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=-2/pi*interpolated_res(omegaprime)*ω/(omegaprime+ω)
    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
end

"Applies the kramers-kronig relations but with scipy's cauchy weight; kwargs for scipy.integrate.quad supported"
function kramers_kronig_scipy(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real}; kwargs...) 
    interpolated_ims = interpol.interp1d(energies, polarizations)
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime) = 2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
    return pyintegrate.quad(cauchy_inner_function, minimum(energies), maximum(energies), weight="cauchy", epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75, wvar= ω; kwargs...)[1]
end

function kramers_kronig_reverse_quadgk(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real}; δ::Real = 0.1, kwargs...) 
    interpolated_res = interpol.interp1d(energies, polarizations)
    inner_function(omegaprime) = -2/pi*interpolated_res(omegaprime)*ω/(omegaprime^2-(ω+1im*δ)^2)
    return real(quadgk(inner_function, minimum(energies), maximum(energies); kwargs...)[1])
end

function kramers_kronig_quadgk(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real}; δ::Real = 0.1, kwargs...) 
    interpolated_ims = interpol.interp1d(energies, polarizations)
    inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime^2-(ω+1im*δ)^2)
    return real(quadgk(inner_function, minimum(energies), maximum(energies); kwargs...)[1])
end
