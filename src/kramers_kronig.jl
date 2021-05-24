"""
$(TYPEDSIGNATURES)
Applies the kramers-kronig relations onto a 1 dimensional array of numbers consisting of the imaginary value of the polarization to return the real value of polarization
"""
function kramers_kronig(ω::Real, im_pol::Vector{<:Real}, max_energy::Real, histogram_width::Real) 
    omegaprime = 
        try             
            @assert length((0:histogram_width*max_energy)) == length(im_pol)
            collect(0:histogram_width*max_energy).*1/histogram_width
        catch 
            collect(0:histogram_width*max_energy)[1:end-1].*1/histogram_width
        end
    sum(1/histogram_width*2/π*im_pol.*omegaprime./(omegaprime.^2 .- (ω+ω*0.003im)^2))
end

function kramers_kronig(im_pol::Function, ω::Real, max_energy_integration=10; kwargs...)
    cauchy_inner_function(omegaprime) = 2/pi*im_pol(omegaprime)*omegaprime/(omegaprime+ω)
    ErrorAbs=1e-20
    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy", 
        epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
end
#=
function kramers_kronig_nosym(im_pol::Function, ω::Real, max_energy_integration=10; kwargs...)
    cauchy_inner_function(omegaprime) = 1/pi*im_pol(omegaprime)
    ErrorAbs=1e-20
    return pyintegrate.quad(cauchy_inner_function, -max_energy_integration, max_energy_integration, weight="cauchy", 
        epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
end

function kramers_kronig_reverse_nosym(im_pol::Function, ω::Real, max_energy_integration=10; kwargs...)
    cauchy_inner_function(omegaprime) = -1/pi*im_pol(omegaprime)
    ErrorAbs=1e-30
    return pyintegrate.quad(cauchy_inner_function, -max_energy_integration, max_energy_integration, weight="cauchy", 
        epsrel=ErrorAbs, epsabs=ErrorAbs, limit=750, wvar= ω; kwargs...)[1]
end
=#
function kramers_kronig_reverse(im_pol::Function, ω::Real, max_energy_integration=10; kwargs...) 
    cauchy_inner_function(omegaprime) = -2/pi*im_pol(omegaprime)*ω/(omegaprime+ω)
    ErrorAbs=1e-20
    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy", 
        epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75, wvar=ω; kwargs...)[1]
end

"""
$(TYPEDSIGNATURES)
Mostly provided for checking the reciprocity relation between the real and imaginary susceptibilities. Give the real susceptibility to obtain the imaginary susceptibility

"""
function kramers_kronig_reverse(ω::Real, re_pol::Vector{<:Real}, max_energy::Real, domega::Real) 
    omegaprime=collect(0:domega:max_energy)
    sum(-domega*2/π*re_pol.*ω./(omegaprime.^2 .- (ω+0.003im)^2))
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
function kramers_kronig_scipy(ω::Real, im_pol::Vector{<:Real}, max_energy::Real, histogram_width::Real, max_energy_integration::Real; kwargs...) 
    interpolated_ims = 
        try 
            interpol.interp1d(0:1/histogram_width:(max_energy), im_pol)
        catch 
            interpol.interp1d(collect(0:1/histogram_width:(max_energy))[0:end-1], im_pol)
        end

    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)
    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]
end

function kramers_kronig_reverse_quadgk(ω::Real, re_pol::Vector{<:Real}, max_energy::Real, domega::Real, max_energy_integration::Real ; δ::Real = 0.1, kwargs...) 
    interpolated_res = 
        try 
            interpol.interp1d(0:domega:max_energy, re_pol)
        catch
            interpol.interp1d(0:domega:max_energy, re_pol)[1:end-1]
        end
    inner_function(omegaprime)=-2/pi*interpolated_res(omegaprime)*ω/(omegaprime^2-(ω+1im*δ)^2)
    return real(quadgk(inner_function, 0, max_energy_integration; kwargs...)[1])
end

function kramers_kronig_quadgk(ω::Real, im_pol::Vector{<:Real}, max_energy::Real, histogram_width::Real, max_energy_integration::Real; δ::Real = 0.1, kwargs...) 
    interpolated_ims = 
        try 
            interpol.interp1d(0:1/histogram_width:(max_energy), im_pol)
        catch 
            interpol.interp1d(0:1/histogram_width:(max_energy)[1:end-1], im_pol)
        end
    inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime^2-(ω+1im*δ)^2)
    return real(quadgk(inner_function, 0, max_energy_integration; kwargs...)[1])
end
