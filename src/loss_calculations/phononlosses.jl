"""
$(TYPEDSIGNATURES)
First order phonon loss through electron phonon interactions. 
"""
function forderphononloss(hwannier::Array{<:Real, 3}, cellmap::Array{<:Real, 2}, HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Integer, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2},
    qnorm::Array{<:Real, 1}, phononband::Integer, μ::Real; histogram_width::Real=100, mesh::Integer=100 )
    ω = phonon_dispersion(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, qnorm::Array{<:Real, 1})[phononband]
    inverselifetime = 0 
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        k1 = [xmesh/mesh, ymesh/mesh, 0]
        k2 = k1 + qnorm
        gsquared = (abs(eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, k1, k2)[phononband]))^2
        ϵ1 = wannier_bands(hwannier, cellmap, k1)
        ϵ2 = wannier_bands(hwannier, cellmap, k2) 
        (abs((ϵ2-ϵ1-ω)*histogram_width) < 0.5) && (inverselifetime += 2π/ħ*heaviside(μ-ϵ1)*(1-heaviside(μ-ϵ2))/(mesh)^2*gsquared*histogram_width)
    end
    return inverselifetime
end