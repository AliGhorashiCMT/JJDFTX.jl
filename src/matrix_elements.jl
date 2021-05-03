"""
$(TYPEDSIGNATURES)
"""
function eph_matrix_elements(HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, k1::Array{<:Real, 1}, k2::Array{<:Real, 1})
    omegaPh, Uph = phonon_dispersionmodes(force_matrix, phonon_cell_map, k1-k2)
    #Note that the phonon energies given by phonon dispersionmodes are in eV, so they must be converted 
    omegaPh *= eV
    phase1 = exp.((2im*π )*(cellMapEph*k1))
    phase2 = exp.((2im*π)*(cellMapEph*k2))
    normFac = np.sqrt(0.5 ./ np.maximum(omegaPh,1e-6))
    g= vec(np.einsum("xy, xab-> yab", Uph, #Rotate to phonon eigenbasis
        np.einsum("R,Rxab->xab", phase2, #Fourier transform from r2 -> k2
        np.einsum("r,rRxab->Rxab", conj(phase1), #Fourier transform from r1 -> k1
        HePhWannier)))).*normFac  #Phonon amplitude factor
    return g/eV
end

function eph_matrix_elements(HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, wannier_file::String, cell_map_file::String, k1::Array{<:Real, 1}, k2::Array{<:Real, 1}, nbands::Int)
    omegaPh, Uph = phonon_dispersionmodes(force_matrix, phonon_cell_map, k1-k2)
    ##Note that the phonon energies given by phonon dispersionmodes are in eV, so they must be converted 
    omegaPh *= eV
    phase1 = exp.((2im*π )*(cellMapEph*k1))
    phase2 = exp.((2im*π)*(cellMapEph*k2))
    normFac = np.sqrt(0.5 ./ np.maximum(omegaPh,1e-6))
    U2 = wannier_vectors(wannier_file::String, cell_map_file::String, k2, nbands) 
    U1 = wannier_vectors(wannier_file::String, cell_map_file::String, k1, nbands) 
    g = np.einsum("bd, ycb-> ycd", U2, #Rotate to electron 2 eigenbasis
    np.einsum("ac,yab -> ycb", conj(U1), #Rotate to electron 1 eigenbasis
    np.einsum("xy, xab-> yab", Uph, #Rotate to phonon eigenbasis
    np.einsum("R, Rxab -> xab", phase2, #Fourier transform from r2 -> k2
    np.einsum("r, rRxab -> Rxab", conj(phase1), #Fourier transform from r1 -> k1
    HePhWannier))))).*normFac #Phonon amplitude factor
    return g/eV
end

function eph_matrix_elements(HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, 
    phonon_cell_map::Array{<:Real, 2}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, 
    k1::Array{<:Real, 1}, k2::Array{<:Real, 1}, nbands::Integer)
    omegaPh, Uph = phonon_dispersionmodes(force_matrix, phonon_cell_map, k1-k2)
    ##Note that the phonon energies given by phonon dispersionmodes are in eV, so they must be converted 
    omegaPh *= eV
    phase1 = exp.((2im*π )*(cellMapEph*k1))
    phase2 = exp.((2im*π)*(cellMapEph*k2))
    normFac = np.sqrt(0.5 ./ np.maximum(omegaPh,1e-4))
    U2 = wannier_vectors(HWannier, cellmap, k2) 
    U1 = wannier_vectors(HWannier, cellmap, k1) 
    g = np.einsum("bd, ycb-> ycd", U2, #Rotate to electron 2 eigenbasis
    np.einsum("ac,yab -> ycb", conj(U1), #Rotate to electron 1 eigenbasis
    np.einsum("xy, xab-> yab", Uph, #Rotate to phonon eigenbasis
    np.einsum("R, Rxab -> xab", phase2, #Fourier transform from r2 -> k2
    np.einsum("r, rRxab -> Rxab", conj(phase1), #Fourier transform from r1 -> k1
    HePhWannier))))).*normFac #Phonon amplitude factor
    return g/eV
end

function plot_ephmatrixelements(HePhWannier::Array{<:Real, 5}, cellMapEph::Array{<:Real, 2}, force_matrix::Array{<:Real, 3}, 
    phonon_cell_map::Array{<:Real, 2}, HWannier::Array{Float64, 3}, cellmap::Array{Float64, 2}, 
    k1::Array{<:Real, 1}, k2::Array{<:Real, 1}, nbands::Integer; kpointfile::String = "bandstruct.kpoints", band1::Integer=4, band2::Integer=5,
    phononband::Integer=5)
    ephmatelements = Vector{Float64}()
    for q in eachrow(np.loadtxt(kpointfile, skiprows=2, usecols=[1, 2, 3]))
        push!(ephmatelements, ((abs(eph_matrix_elements(HePhWannier, cellMapEph, force_matrix, phonon_cell_map, HWannier, cellmap, [0, 0, 0], collect(q), nbands)[phononband, band1, band2]))^2))
    end
    plot!(ephmatelements, linewidth=5)
end

"""
$(TYPEDSIGNATURES)
"""
function momentum_matrix_elements(Pwannier::Array{Float64, 4}, cell_map::Array{Float64, 2}, k::Array{<:Real, 1})
    phase = np.exp(2im*π*cell_map*k); 
    Pk = np.tensordot(phase, Pwannier, axes=1); 
    #= 
    JDFTX output is in Atomic units. Therefore, the units of the momentum matrix elements are in ħ/a₀
    a₀ is the Bohr radius, which is approximately 0.529 Angstrom. ħ is 6.6*10-16 eV*seconds. Since all 
    physical calculations perfomed in jdftx_to_plot assume lengths to be given in angstroms and energies 
    to be given in eV, we multiply Pk by ħ/bohrtoangstrom
    =#
    return Pk*ħ/bohrtoangstrom
end

function momentum_matrix_elements(Hwannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, Pwannier::Array{Float64, 4}, k::Array{<:Real, 1})
    phase = np.exp(2im*π*cell_map*k); 
    Pk = np.tensordot(phase, Pwannier, axes=1); 
    Vk = wannier_vectors(Hwannier, cell_map, k) 
    Pk = np.einsum("ba, pbc, cd-> pad",   #Sum using Einstein notation for
    conj(Vk), Pk, Vk)  
    #= 
    JDFTX output is in Atomic units. Therefore, the units of the momentum matrix elements are in ħ/a₀
    a₀ is the Bohr radius, which is approximately 0.529 Angstrom. ħ is 6.6*10-16 eV*seconds. Since all 
    physical calculations perfomed in jdftx_to_plot assume lengths to be given in angstroms and energies 
    to be given in eV, we multiply Pk by ħ/bohrtoangstrom
    =#
    return Pk*ħ/bohrtoangstrom
end

"""
$(TYPEDSIGNATURES)
"""
function pwannier(pwannier_file::String, cell_map_file::String, nbands::Int64) 
    cell_map = np.loadtxt(cell_map_file)
    cell_map_numlines = countlines(cell_map_file)
    Pwannier = np.reshape(np.loadtxt(pwannier_file), (cell_map_numlines, 3, nbands, nbands))
    return Pwannier
end

function pwannier(pwannier_file::String, cell_map_file::String) 
    cell_map = np.loadtxt(cell_map_file);
    cell_map_numlines = countlines(cell_map_file);
    nbands = np.int(np.sqrt(np.size(np.loadtxt(pwannier_file)) //(cell_map_numlines*3)));
    println("The number of bands detected is: ", nbands, "\nIf this is incorrect, something went wrong at some point somewhere")
    Pwannier = np.reshape(np.loadtxt(pwannier_file), (cell_map_numlines, 3, nbands, nbands));
    return Pwannier;
end


