"""
$(TYPEDSIGNATURES)
Applies the kramers-kronig relations onto a 1 dimensional array of numbers consisting of the imaginary value of the polarization to return the real value of polarization
"""
function kramers_kronig(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real},
    real_or_imaginary::Union{Val{:real}, Val{:imaginary}} = Val(:imaginary); δ::Real=0.01)
    kramers_kronig_result = 
        if real_or_imaginary == Val(:imaginary)
            2/π*sum(diff(energies).*(polarizations.*energies./(energies.^2 .- (ω+ δ*1im)^2))[1:end-1])
        elseif real_or_imaginary == Val(:real)
            -2/π*sum(diff(energies).*(polarizations.*ω./(energies.^2 .- (ω+ δ*1im)^2))[1:end-1])
        end
        return kramers_kronig_result
end

function kramers_kronig(polarization_function::Function, ω::Real, max_energy_integration::Real=10, 
    real_or_imaginary::Union{Val{:real}, Val{:imaginary}} = Val(:imaginary); min_energy_integration::Real=0, kwargs...)
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime) = 
        if real_or_imaginary == Val(:imaginary)
            2/pi*polarization_function(omegaprime)*omegaprime/(omegaprime+ω)
        elseif real_or_imaginary == Val(:real)
            -2/pi*im_pol(omegaprime)*ω/(omegaprime+ω)
        end
    return pyintegrate.quad(cauchy_inner_function, min_energy_integration, max_energy_integration, weight="cauchy", 
        epsrel=ErrorAbs, epsabs=ErrorAbs, limit=175,  wvar= ω ; kwargs...)[1]
end

"Applies the kramers-kronig relations but with scipy's cauchy weight; kwargs for scipy.integrate.quad supported"
function kramers_kronig_scipy(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real},
    real_or_imaginary::Union{Val{:real}, Val{:imaginary}} = Val(:imaginary); kwargs...) 
    interpolation_function = interpol.interp1d(energies, polarizations)
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime) = 
        if real_or_imaginary == Val(:imaginary)
            2/pi*interpolation_function(omegaprime)*omegaprime/(omegaprime+ω)
        elseif real_or_imaginary == Val(:real)
            -2/pi*interpolation_function(omegaprime)*ω/(omegaprime+ω)
        end
        upper_limit = energies[findall(x -> !isapprox(x, 0), polarizations)[end]]
    return pyintegrate.quad(cauchy_inner_function, minimum(energies), upper_limit, weight="cauchy", epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75, wvar= ω; kwargs...)[1]
end

function kramers_kronig_quadgk(ω::Real, energies::Vector{<:Real}, polarizations::Vector{<:Real}, 
    real_or_imaginary::Union{Val{:real}, Val{:imaginary}} = Val(:imaginary); δ::Real = 0.1, kwargs...) 
    interpolation_function = interpol.interp1d(energies, polarizations)
    inner_function(omegaprime) = 
        if real_or_imaginary == Val(:imaginary)
            2/pi*interpolation_function(omegaprime)*omegaprime/(omegaprime^2-(ω+1im*δ)^2)
        elseif real_or_imaginary == Val(:real)
            -2/pi*interpolation_function(omegaprime)*ω/(omegaprime^2-(ω+1im*δ)^2)
        end
        upper_limit = energies[findall(x -> !isapprox(x, 0), polarizations)[end]]
    return real(first(quadgk(inner_function, minimum(energies), upper_limit; kwargs...)))
end
