# All utilities for parsing jdftx output

"""
Assuming a planar crystal with some out of plane ion, get the distance of the ion from the plane (convention is that 
the plane is at z=0 and that the z direction is perpendicular to the planar crystal)
"""
function get_d(filename::AbstractString, ion::AbstractString)
    lattice_z = loadlattice(filename)[3][3]
    return parse.(Float64, split(filter(line -> contains(line, "ion $ion"), readlines(filename))[1])[5])*lattice_z
 end

 """
 Get the- possibly fixed- initial magnetization.  
 """
 function get_mag(filename::AbstractString)
    magnetization_line = split(first(filter(line->contains(line, "elec-initial-magnetization"), readlines(filename))))
    get_fixed = magnetization_line[3]
    get_fixed == "no" && @warn "Magnetization not fixed- use with care"
    return parse(Float64, magnetization_line[2])
end

"""
Outputs the result of the listEnergy script of jdftx
"""
function list_energy(filename::AbstractString)
    e = readlines(`listEnergy $filename`)
    e = parse(Float64, first(split(first(e))))
    return e/eV
end