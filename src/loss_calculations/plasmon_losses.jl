"""
$(TYPEDSIGNATURES)
"""
function landau_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Vector{<:Vector{<:Real}}, 
    histogram_width::Integer, mesh::Integer, q::Vector{<:Real}, μ::Real, energy_range::Real) 

    lossarray = zeros(histogram_width*energy_range)
    ucellarea = unit_cell_area(lattice_vectors)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵ1 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
        ϵ2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
        f1 = ϵ1 < μ ? 1 : 0
        f2 = ϵ2 > μ ? 1 : 0
        ω = ϵ2-ϵ1
        (ω > 0) && (lossarray[round(Int, ω*histogram_width)+1 ] += 2π/ħ*e²ϵ/4*ω/qabs*1/ucellarea*f1*f2*(1/mesh)^2*histogram_width)
    end
    return lossarray
end


#=
Next, we will do the same as the above, except with a full phonon band calculation (calculate over all phonon branches, taking into account 
the different electron phonon matrix elements )
=#
"""
$(TYPEDSIGNATURES)

Compute the first order damping of plasmons due to the electron-phonon interaction
"""
function first_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, 
    lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; histogram_width::Real=100, mesh::Integer=30, energy_range::Real=10, subsampling::Real=1, offset::Vector{<:Real}=[0, 0, 0]) 
    lossarray = zeros(histogram_width*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        kinitial = [xmesh/(subsampling*mesh), ymesh/(subsampling*mesh), 0] + offset
        ϵinitial = wannier_bands(HWannier, cell_map, kinitial)
        ϵmiddle = wannier_bands(HWannier, cell_map, kinitial + qnormalized )        
        finitial = ϵinitial<μ ? 1 : 0
        fmiddle1 = ϵmiddle>μ ? 1 : 0
        for (xmesh1, ymesh1) in Tuple.(CartesianIndices(rand(mesh, mesh)))
            phonon_energies = phonon_dispersion(force_matrix, phonon_cell_map, [xmesh1/mesh, ymesh1/mesh, 0])
            phonon_mat_elements1= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, kinitial, kinitial + [xmesh1/mesh, ymesh1/mesh, 0])
            phonon_mat_elements2= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, kinitial + qnormalized , kinitial+[xmesh1/mesh, ymesh1/mesh, 0]+qnormalized)
            for phonon in 1:length(phonon_energies)
                g1= phonon_mat_elements1[phonon]
                g2= phonon_mat_elements2[phonon]
                ϵphonon = phonon_energies[phonon]
                ϵmiddle2 = wannier_bands(HWannier, cell_map, kinitial+[xmesh1/mesh, ymesh1/mesh, 0])        
                fmiddle2 = ϵmiddle2>μ ? 1 : 0
                ϵfinal = wannier_bands(HWannier, cell_map, kinitial+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                ffinal = ϵfinal>μ ? 1 : 0
                ω = ϵfinal-ϵinitial+ϵphonon
                ω>0 && (lossarray[round(Int, ω*histogram_width+ 1 )] += 1/cell_area*abs(fmiddle1*g2/(ϵmiddle-ϵinitial-ω)+ fmiddle2*g1/(ϵmiddle2-ϵinitial+ϵphonon))^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^4*histogram_width)
            end
        end
    end
    return lossarray
end

"""
$(TYPEDSIGNATURES)
"""
function second_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, 
    lattice_vectors::Vector{<:Vector{<:Real}}, q::Vector{<:Real}, μ::Real; histogram_width::Real=100, mesh::Integer=30, energy_range::Real=10) 
    #= In all cases, we consider an initial state with an electron at kx, ky = xmesh/mesh, ymesh/mesh and a plasmon with wavevector q 
     Furthermore, the final state is an electron with momentum (xmesh+xmesh1+xmesh2)/mesh, (ymesh+ymesh1+ymesh2)/mesh - qnormalized
    as well as two emitted phonons. We sum over all intermmediate states and square the sum in the integrand. The three possible decay channels are 
    plasmon absorption first, second, or third. Note that thus for each iteration of the sum, we're summing contributions from virtual 
    electronic intermmediate states. Note also that the matrix elements included here contain Fermi occupation functions 
     =#
    lossarray = zeros(histogram_width*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)
    for (xmesh, ymesh) in Tuple.(CartesianIndices(rand(mesh, mesh)))
        ϵinitial = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
            #first middle state (plasmon absorbed first)
        ϵmiddle = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        
        finitial = ϵinitial < μ ? 1 : 0
        finitial == 1 || continue
        println(finitial)
        flush(stdout)
        fmiddle1 = ϵmiddle>μ ? 1 : 0
        for (xmesh1, ymesh1) in Tuple.(CartesianIndices(rand(mesh, mesh)))
            #second middle state (phonon absorbed first)
            ϵmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
            fmiddle2 = ϵmiddle2>μ ? 1 : 0
            #Then consider second middle states from one phonon and plasmon absoption first, (phonon absorbed last)
            ϵsecondmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0] + qnormalized)        
            fsecondmiddle2 = ϵsecondmiddle2>μ ? 1 : 0
            for (xmesh2, ymesh2) in Tuple.(CartesianIndices(rand(mesh, mesh)))
                #first consider second middle states originating from 2 phonons being absorbed first (plasmon absorbed last)
                ϵsecondmiddle1 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0])        
                fsecondmiddle1 = ϵsecondmiddle1>μ ? 1 : 0
                #The final energy will always be that of the electronic state corresponding to the original kvector plus the two phonon kvectors and the plasmon kvector
                ϵfinal = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                ffinal = ϵfinal>μ ? 1 : 0
                ffinal == 1 || continue
                phonon_energies1 = phonon_dispersion(force_matrix, phonon_cell_map, [xmesh1/mesh, ymesh1/mesh, 0])
                phonon_energies2 = phonon_dispersion(force_matrix, phonon_cell_map, [xmesh2/mesh, ymesh2/mesh, 0])
                #=
                    There are four possible matrix elements, corresponding to the three possible phonon assisted transitions
                =#
                phonon_mat_elements1= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map,[xmesh/mesh, ymesh/mesh, 0], [xmesh/mesh+xmesh1/mesh, ymesh/mesh+ymesh1/mesh, 0])
                phonon_mat_elements2= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map,[xmesh/mesh, ymesh/mesh, 0]+qnormalized, [xmesh/mesh+xmesh1/mesh, ymesh/mesh+ymesh1/mesh, 0]+qnormalized)
                phonon_mat_elements3= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map,[xmesh/mesh+xmesh1/mesh, ymesh/mesh+ymesh1/mesh, 0], [xmesh/mesh+xmesh1/mesh+xmesh2/mesh, ymesh/mesh+ymesh1/mesh+ymesh2/mesh, 0])
                phonon_mat_elements4= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map,[xmesh/mesh+xmesh1/mesh, ymesh/mesh+ymesh1/mesh, 0]+qnormalized, [xmesh/mesh+xmesh1/mesh+xmesh2/mesh, ymesh/mesh+ymesh1/mesh+ymesh2/mesh, 0]+qnormalized)
                for first_phonon in 1:length(phonon_energies1)
                    for second_phonon in 1:length(phonon_energies2)
                        ϵ1 = phonon_energies1[first_phonon]
                        ϵ2 = phonon_energies2[second_phonon]
                        ω = ϵfinal-ϵinitial+(ϵ1+ϵ2)
                        g1 = phonon_mat_elements2[first_phonon] #Plasmon has been absorbed first, phonon absorbed second
                        g2 = phonon_mat_elements4[second_phonon] #Plasmon absorbed first, phonon absorbed second, phonon absorbed thrid
                        g3 = phonon_mat_elements1[first_phonon] #phonon absorbed first
                        g4 = phonon_mat_elements3[second_phonon] #phonon absorbed first, phonon absorbed second
                        g5 = g2 #phonon absorbed first, plasmon absorbed second, phonon absorbed third
                            #=
                                First term is plasmon first, phonon second, phonon third
                                Second term is phonon first, phonon second, plasmon third
                                Last term is phonon first, plasmon second, phonon third
                            =#
                        ω > 0 || continue
                        lossarray[round(Int, ω*histogram_width + 1 )] += 1/cell_area*abs(g1/(ϵmiddle-ϵinitial-ω-.01*1im)*g2/(ϵsecondmiddle2-ϵinitial-ω-.01*1im+ϵ1)*fmiddle1*fsecondmiddle2 + fmiddle2*g3/(ϵmiddle2-ϵinitial+ϵ1+.01*1im)*( g4*fsecondmiddle1/(ϵsecondmiddle1-ϵinitial+.01*im+ϵ1+ϵ2) + g5*fsecondmiddle2/(ϵsecondmiddle2-ϵinitial-ω-.01*im+ϵ1) ))^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^6*histogram_width
                    end
                end
            end
        end
    end
    return lossarray
end




