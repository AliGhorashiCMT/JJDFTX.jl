# All utilities for parsing jdftx output

"""
Assuming a planar crystal with some out of plane ion, get the distance of the ion from the plane (convention is that 
the plane is at z=0 and that the z direction is perpendicular to the planar crystal)
"""
function get_d(filename::AbstractString, ion::AbstractString)
    return parse.(Float64, split(filter(line -> contains(line, "ion $ion"), readlines(filename))[1])[5])*bohr*40 
 end

 """
 Get the- possibly fixed- initial magnetization.  
 """
 function get_mag(filename::AbstractString)
    return parse(Float64, split(filter(line->contains(line, "elec-initial-magnetization"), readlines(filename))[1])[2])
end