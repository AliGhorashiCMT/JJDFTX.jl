
function te_plasmon(σ::Complex, ω::Real, q::Real)
    return 1-1im*ω*σ/(2*137*2000)*1/(sqrt(q^2-ω^2/2000^2))
end