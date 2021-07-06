
function te_plasmon(σ::Complex, ω::Real, q::Real)
    return 1-1im*ω*σ/(2*137*(ħ*c))*1/(sqrt(q^2-ω^2/(ħ*c)^2))
end

function te_plasmon(σ::Complex, ω::Real, npoints::Integer=1000, frac::Real=1e7)
    qlist = collect(1:(1/(npoints*frac)):(1+1/(frac)))
    plas_dispersion = zeros(length(qlist))
    for (i, q) in enumerate(qlist)
        q = ω/(ħ*c)*(1+i/(frac*npoints))
        plas_dispersion[i] = log(abs(real(te_plasmon(σ, ω, q))))
    end
    return qlist, plas_dispersion
end

function multilayer_te_plasmon(σ::Complex, ω::Real, npoints::Integer=1000, frac::Real=1e7, n::Integer=1)
    qlist = collect(1:(1/(npoints*frac)):(1+1/(frac)))
    plas_dispersion = zeros(length(qlist))
    for (i, q) in enumerate(qlist)
        q = ω/(ħ*c)*(1+i/(frac*npoints))
        plas_dispersion[i] = log(abs(real(   1+n*(te_plasmon(σ, ω, q)-1)    )))
    end
    return qlist, plas_dispersion
end