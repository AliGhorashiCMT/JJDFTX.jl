"""
$(TYPEDSIGNATURES)
Returns the electron phonon coupling, h(ϵ) as defined in Ab initio phonon coupling and optical response of hot electrons in plasmonic metals
"""
function ephcoupling()
# Note: Our densities of states will not include 1/Volume factors 
# Therefore, the six dimensional integral is just a 1/(Nk*Nk') times an average over the Brillouin zone. 
# We use a histogramming method for the product of the two delta functions. 
# We also have an extra factor of 4 compared to the cited paper due to the spin degeneracy that we take into account here

#2*gfermi/gϵ^2 

end