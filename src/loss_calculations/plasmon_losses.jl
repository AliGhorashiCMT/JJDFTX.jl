#We calculate losses at 0th (landau damping), 1st, and 2nd orders in phonon-assisted damping

function landau_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1},1}, histogram_width::Int, mesh::Int, q::Array{<:Real, 1}, μ::Real, energy_range::Real) 
    lossarray = zeros(histogram_width*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            ϵ1 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0])
            ϵ2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
            f1 = ϵ1<μ ? 1 : 0
            f2 = ϵ2>μ ? 1 : 0
            if f1>0 && f2>0
                ω = ϵ2-ϵ1
                lossarray[round(Int, (ω+offset)*histogram_length  )] = lossarray[round(Int, (ω+offset)*histogram_length  )] + 2π/ħ*e²ϵ/4*ω/qabs*f1*f2*(1/mesh)^2*histogram_length
            end
        end
    end
    return lossarray
end

function landau_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, histogram_width::Int, mesh::Int, q::Array{<:Real, 1}, μ::Real, energy_range::Real) 
    lossarray = zeros(histogram_width*energy_range)
    ucellarea = unit_cell_area(lattice_vectors)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            ϵ1 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
            ϵ2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
            f1 = ϵ1<μ ? 1 : 0
            f2 = ϵ2>μ ? 1 : 0
            if f1>0 && f2>0
                ω = ϵ2-ϵ1
                lossarray[round(Int, ω*histogram_width)+1 ] = lossarray[round(Int, ω*histogram_width)+1] + 2π/ħ*e²ϵ/4*ω/qabs*1/ucellarea*f1*f2*(1/mesh)^2*histogram_width
            end
        end
    end
    return lossarray
end

function first_order_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real, ϵphonon::Real, gph::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            ϵinitial = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0])
            ϵmiddle = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0
            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    ϵmiddle2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    ϵfinal = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                    ffinal = ϵfinal>μ ? 1 : 0

                    ω = ϵfinal-ϵinitial+ϵphonon
                    if ω>0
                        lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*(1/(ϵmiddle-ϵinitial-ω)+1/(ϵmiddle2-ϵinitial+ϵphonon))^2*gph^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*fmiddle1*ffinal*(1/mesh)^4*histogram_length
                    end
                end
            end
        end
    end
    return lossarray
end

function first_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real, ϵphonon::Real, gph::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            ϵinitial = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
            ϵmiddle = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0
            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    ϵmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    ϵfinal = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                    ffinal = ϵfinal>μ ? 1 : 0

                    ω = ϵfinal-ϵinitial+ϵphonon
                    if ω>0
                        lossarray[round(Int, ω*histogram_length + 1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*(fmiddle1/(ϵmiddle-ϵinitial-ω)+ fmiddle2/(ϵmiddle2-ϵinitial+ϵphonon))^2*gph^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^4*histogram_length
                    end
                end
            end
        end
    end
    return lossarray
end



#=
Next, we will do the same as the above, except with a full phonon band calculation (calculate over all phonon branches, taking into account 
the different electron phonon matrix elements )
=#

function first_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
lossarray = zeros(histogram_length*energy_range)
qabs = sqrt(sum(q.^2))
qnormalized = normalize_kvector(lattice_vectors, q)
cell_area = unit_cell_area(lattice_vectors)
for xmesh in 1:mesh
    for ymesh in 1:mesh
        ϵinitial = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
        ϵmiddle = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

        finitial = ϵinitial<μ ? 1 : 0
        fmiddle1 = ϵmiddle>μ ? 1 : 0
        for xmesh1 in 1:mesh
            for ymesh1 in 1:mesh

                phonon_energies = phonon_dispersion(force_matrix, phonon_cell_map, [xmesh1/mesh, ymesh1/mesh, 0])

                phonon_mat_elements1= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map,[xmesh/mesh, ymesh/mesh, 0], [xmesh/mesh+xmesh1/mesh, ymesh/mesh+ymesh1/mesh, 0])
                phonon_mat_elements2= eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map,[xmesh/mesh, ymesh/mesh, 0]+qnormalized , [xmesh/mesh+xmesh1/mesh, ymesh/mesh+ymesh1/mesh, 0]+qnormalized)

                for phonon in 1:length(phonon_energies)

                    g1= phonon_mat_elements1[phonon]
                    g2= phonon_mat_elements2[phonon]

                    ϵphonon = phonon_energies[phonon]

                    ϵmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    ϵfinal = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                
                    ffinal = ϵfinal>μ ? 1 : 0

                    ω = ϵfinal-ϵinitial+ϵphonon
                    if ω>0
                        lossarray[round(Int, ω*histogram_length + 1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*abs(fmiddle1*g2/(ϵmiddle-ϵinitial-ω)+ fmiddle2*g1/(ϵmiddle2-ϵinitial+ϵphonon))^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^4*histogram_length
                    end
                end
            end
        end
    end
end
return lossarray
end





function second_order_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real, ϵphonon::Real, gph::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    

    #= In all cases, we consider an initial state with an electron at kx, ky = xmesh/mesh, ymesh/mesh and a plasmon with wavevector q 
     Furthermore, the final state is an electron with momentum (xmesh+xmesh1+xmesh2)/mesh, (ymesh+ymesh1+ymesh2)/mesh - qnormalized
    as well as two emitted phonons. We sum over all intermmediate states and square the sum in the integrand. The three possible decay channels are 
    plasmon absorption first, second, or third. Note that thus for each iteration of the sum, we're summing contributions from virtual 
    electronic intermmediate states. Note also that the matrix elements included here contain Fermi occupation functions 
     =#

    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)

    for xmesh in 1:mesh
        for ymesh in 1:mesh

            ϵinitial = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0])
            #first middle state (plasmon absorbed first)
            ϵmiddle = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0

            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    #second middle state (phonon absorbed first)
                    ϵmiddle2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    #Then consider second middle states from one phonon and plasmon absoption first, (phonon absorbed last)

                    ϵsecondmiddle2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0] + qnormalized)        
                    fsecondmiddle2 = ϵsecondmiddle1>μ ? 1 : 0


                    for xmesh2 in 1:mesh
                        for ymesh2 in 1:mesh

                            #first consider second middle states originating from 2 phonons being absorbed first (plasmon absorbed last)
                            ϵsecondmiddle1 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0])        
                            fsecondmiddle1 = ϵsecondmiddle1>μ ? 1 : 0

                            #The final energy will always be that of the electronic state corresponding to the original kvector plus the two phonon kvectors and the plasmon kvector
                            ϵfinal = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                            ffinal = ϵfinal>μ ? 1 : 0

                            ω = ϵfinal-ϵinitial+ϵphonon
                            if ω>0
                                lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*( 1/(ϵmiddle-ϵinitial-ω)*1/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon)*fmiddle1*fsecondmiddle2 + fmiddle2/(ϵmiddle2-ϵinitial+ϵphonon)*( fsecondmiddle1/(ϵsecondmiddle1-ϵinitial+2*ϵphonon) + fsecondmiddle2/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon) ))^2*gph^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^6*histogram_length
                            end
                        end
                    end
                end
            end
        end
    end
    return lossarray
end


function second_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real, ϵphonon::Real, gph::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    

    #= In all cases, we consider an initial state with an electron at kx, ky = xmesh/mesh, ymesh/mesh and a plasmon with wavevector q 
     Furthermore, the final state is an electron with momentum (xmesh+xmesh1+xmesh2)/mesh, (ymesh+ymesh1+ymesh2)/mesh - qnormalized
    as well as two emitted phonons. We sum over all intermmediate states and square the sum in the integrand. The three possible decay channels are 
    plasmon absorption first, second, or third. Note that thus for each iteration of the sum, we're summing contributions from virtual 
    electronic intermmediate states. Note also that the matrix elements included here contain Fermi occupation functions 
     =#

    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)

    for xmesh in 1:mesh
        for ymesh in 1:mesh

            ϵinitial = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
            #first middle state (plasmon absorbed first)
            ϵmiddle = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0

            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    #second middle state (phonon absorbed first)
                    ϵmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    #Then consider second middle states from one phonon and plasmon absoption first, (phonon absorbed last)

                    ϵsecondmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0] + qnormalized)        
                    fsecondmiddle2 = ϵsecondmiddle2>μ ? 1 : 0


                    for xmesh2 in 1:mesh
                        for ymesh2 in 1:mesh

                            #first consider second middle states originating from 2 phonons being absorbed first (plasmon absorbed last)
                            ϵsecondmiddle1 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0])        
                            fsecondmiddle1 = ϵsecondmiddle1>μ ? 1 : 0

                            #The final energy will always be that of the electronic state corresponding to the original kvector plus the two phonon kvectors and the plasmon kvector
                            ϵfinal = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                            ffinal = ϵfinal>μ ? 1 : 0

                            ω = ϵfinal-ϵinitial+2*ϵphonon
                            if ω>0
                                lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*( 1/(ϵmiddle-ϵinitial-ω)*1/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon)*fmiddle1*fsecondmiddle2 + fmiddle2/(ϵmiddle2-ϵinitial+ϵphonon)*( fsecondmiddle1/(ϵsecondmiddle1-ϵinitial+2*ϵphonon) + fsecondmiddle2/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon) ))^2*gph^4*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^6*histogram_length
                            end
                        end
                    end
                end
            end
        end
    end
    return lossarray
end


function second_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    

    #= In all cases, we consider an initial state with an electron at kx, ky = xmesh/mesh, ymesh/mesh and a plasmon with wavevector q 
     Furthermore, the final state is an electron with momentum (xmesh+xmesh1+xmesh2)/mesh, (ymesh+ymesh1+ymesh2)/mesh - qnormalized
    as well as two emitted phonons. We sum over all intermmediate states and square the sum in the integrand. The three possible decay channels are 
    plasmon absorption first, second, or third. Note that thus for each iteration of the sum, we're summing contributions from virtual 
    electronic intermmediate states. Note also that the matrix elements included here contain Fermi occupation functions 
     =#

    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)

    for xmesh in 1:mesh
        for ymesh in 1:mesh

            ϵinitial = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
            #first middle state (plasmon absorbed first)
            ϵmiddle = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0

            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    #second middle state (phonon absorbed first)
                    ϵmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    #Then consider second middle states from one phonon and plasmon absoption first, (phonon absorbed last)

                    ϵsecondmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0] + qnormalized)        
                    fsecondmiddle2 = ϵsecondmiddle2>μ ? 1 : 0


                    for xmesh2 in 1:mesh
                        for ymesh2 in 1:mesh

                            #first consider second middle states originating from 2 phonons being absorbed first (plasmon absorbed last)
                            ϵsecondmiddle1 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0])        
                            fsecondmiddle1 = ϵsecondmiddle1>μ ? 1 : 0

                            #The final energy will always be that of the electronic state corresponding to the original kvector plus the two phonon kvectors and the plasmon kvector
                            ϵfinal = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                            ffinal = ϵfinal>μ ? 1 : 0

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

                                    if ω>0
                                        #=
                                            First term is plasmon first, phonon second, phonon third
                                            Second term is phonon first, phonon second, plasmon third
                                            Last term is phonon first, plasmon second, phonon third
                                        =#
                                        lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*abs( g1/(ϵmiddle-ϵinitial-ω)*g2/(ϵsecondmiddle2-ϵinitial-ω+ϵ1)*fmiddle1*fsecondmiddle2 + fmiddle2*g3/(ϵmiddle2-ϵinitial+ϵ1)*( g4*fsecondmiddle1/(ϵsecondmiddle1-ϵinitial+ϵ1+ϵ2) + g5*fsecondmiddle2/(ϵsecondmiddle2-ϵinitial-ω+ϵ1) ))^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^6*histogram_length
                                    end 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return lossarray
end




