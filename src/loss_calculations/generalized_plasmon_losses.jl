"""
$(TYPEDSIGNATURES)

Calculates the damping of plasmons at lowest order in velocity gauge. 
"""
function landau_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, histogram_width::Integer, 
    mesh::Integer, q::Vector{<:Real}, μ::Real, energy_range::Real, nbands::Integer) 
    lossarray = zeros(histogram_width*energy_range)
    ucellarea = unit_cell_area(lattice_vectors)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵ1s = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0], nbands)
        ϵ2s = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized, nbands)
        ωs = [ϵ2 - ϵ1 for ϵ1 in ϵ1s for ϵ2 in ϵ2s]
        ϵtuples = [(ϵ1, ϵ2) for ϵ1 in ϵ1s for ϵ2 in ϵ2s]
        for ((ϵ1, ϵ2), ω) in zip(ϵtuples,  ωs )
            f1 = ϵ1<μ ? 1 : 0
            f2 = ϵ2>μ ? 1 : 0
            (f1>0 && f2>0) ? lossarray[round(Int, ω*histogram_width)+1 ] = lossarray[round(Int, ω*histogram_width)+1] + 2π/ħ*e²ϵ/4*ω/qabs*1/ucellarea*f1*f2*(1/mesh)^2*histogram_width : continue
        end
    end
    return lossarray
end
