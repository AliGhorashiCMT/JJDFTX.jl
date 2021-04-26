function landau_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, histogram_width::Integer, mesh::Integer, q::Vector{<:Real}, μ::Real, energy_range::Real, nbands::Integer) 
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

function first_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real, nbands::Integer; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10, subsampling::Int=1, offset::Vector{<:Real}=zeros(3)) 
    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kinitial = [xmesh/(subsampling*mesh), ymesh/(subsampling*mesh), 0] + offset
        ϵinitials = wannier_bands(HWannier, cell_map, kinitial, nbands)
        ϵmiddles = wannier_bands(HWannier, cell_map, kinitial+qnormalized, nbands)        
        finitials = ϵinitials .< μ 
        fmiddles1 = ϵmiddles .> μ 
        for (xmesh1, ymesh1) in Tuple.(CartesianIndices(rand(mesh, mesh)))
            kphonon = [xmesh1/(subsampling*mesh), ymesh1/(subsampling*mesh), 0]
            phonon_energies = phonon_dispersion(force_matrix, phonon_cell_map, kphonon)
            phonon_mat_elements1= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, HWannier, cell_map, kinitial, kinitial+kphonon,  nbands)
            phonon_mat_elements2= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, HWannier, cell_map, kinitial+qnormalized, kinitial+kphonon+qnormalized, nbands)
            for (g1, g2, ϵphonon) in zip(eachslice(phonon_mat_elements1; dims=1), eachslice(phonon_mat_elements2; dims=1), phonon_energies)
                ϵmiddles2 = wannier_bands(HWannier, cell_map, kinitial+kphonon, nbands)        
                fmiddles2 = ϵmiddles2 .> μ
                ϵfinals = wannier_bands(HWannier, cell_map, kinitial+kphonon+qnormalized, nbands)        
                ffinals = ϵfinals .> μ
                for (findex, (ϵfinal, ffinal)) in enumerate(zip(ϵfinals, ffinals))
                    for (iindex, (ϵinitial, finitial)) in enumerate(zip(ϵinitials, finitials))
                        T=0 
                        ω = ϵfinal-ϵinitial+ϵphonon
                        if ω<0 || finitial*ffinal==0 
                            continue
                        end
                        ω<1 ? println(ω) : nothing
                        for (mindex, (ϵmiddle, fmiddle1, ϵmiddle2, fmiddle2)) in enumerate(zip(ϵmiddles, fmiddles1, ϵmiddles2, fmiddles2))
                            T += fmiddle1*g2[mindex, findex]/(ϵmiddle-ϵinitial-ω*(1+.01im))+ fmiddle2*g1[mindex, iindex]/(ϵmiddle2-ϵinitial+ϵphonon*(1+.01im))
                        end
                        ω > 0 ? lossarray[round(Int, ω*histogram_length+1)] += 1/cell_area*(abs(T)^2)*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^4*histogram_length*(1/subsampling^4) : nothing
                    end
                end
            end
        end
    end
    return lossarray
end
    