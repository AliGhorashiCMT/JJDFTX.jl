module JJDFTX
#dependencies
using Base: Real
using PyPlot
using PyCall
using LinearAlgebra 
using Distances
using HCubature
using QuadGK
using DSP
using DelimitedFiles
using Setfield
using DocStringExtensions 
# ---------------------------------------------------------------------------------------- #
import Base: *, ∘, +

#= We define all constants required for future calculations.
hbar is given in ev*seconds. 
c, the speed of light, is given in angstroms/second
α is the fine structure constant (unitless)
e^2ϵ is actually e^2/ϵ and is given in units of ev*angstrom
the constant eV is the conversion from hartrees to eV.
Note that jdftx output is always in hartrees (atomic units)
=#

# ---------------------------------------------------------------------------------------- #
const ħ = 6.6e-16
const c = 3e18
const α = 1/137
const e²ϵ = 4π*ħ*c*α  
const bohrtoangstrom = 0.529177
const eV = 1/27.2114 
const kB = 8.617333262145e-5
const mₑ = 0.511e6/c^2
# ---------------------------------------------------------------------------------------- #
export ħ,c, α, e²ϵ, bohrtoangstrom, eV, kB, mₑ 
# ---------------------------------------------------------------------------------------- #

const np = PyNULL()
const spatial = PyNULL()
const interpol = PyNULL()
const pyintegrate = PyNULL()
const helper_scripts = PyNULL()
const ase = PyNULL()
export np, interpol, pyintegrate, spatial, helper_scripts
# ---------------------------------------------------------------------------------------- #
##Code that is loaded when jdftx_to_plot is loaded 
function __init__()
    copy!(np, pyimport_conda("numpy", "numpy"))
    copy!(interpol, pyimport_conda("scipy.interpolate", "scipy"))
    copy!(pyintegrate, pyimport_conda("scipy.integrate", "scipy"))
    copy!(spatial, pyimport_conda("scipy.spatial", "scipy"))
    try
        copy!(helper_scripts, pyimport("helper_scripts"))
        println("Congratulations- Your version of python is linked to helper_scripts.py")
    catch 
        println("Your version of python is not linked to helper_scripts.py")
    end  
    try 
        copy!(ase, pyimport_conda("ase", "ase", "conda-forge"))
        println("Congratulations- Your version of python is linked to the Atomic Simulation Environment")
    catch 
        println("Your version of python is not linked to the Atomic Simulation Environment")
    end
end

#= 
Properties of the unit cell, the real space lattice and the reciprocal lattice. 
Methods provided for normalizing the kvectors in the basis of reciprocal lattice vectors (necessary for jdftx).
=#
include("cell_properties.jl")
export unnormalize_kvector, normalize_kvector, brillouin_zone_area,
in_wigner_seitz, in_brillouin, reciprocal_vectors, ion_positions, plot_lattice, 
cell_vectors, unit_cell_area, unit_cell_volume, brillouin_zone_volume, brillouin_zone_volume_direct,
loadlattice, loadreciprocallattice, loadcellarea, loadcellvolume

include("phonon_properties.jl")
export phonon_dispersion, phonon_force_matrix

#=
A multitude of methods for calculating the density of states- either directly from DFT output or from wannier 
tight binding data. Functions also provided for the density of states of phonons. Lastly, a handy method has been provided
for calculation of the chemical potential at arbitrary filling. 
=#
include("density_of_states.jl")
export density_of_states, find_chemical_potential, phonon_dos, bands_overlayed_dos, 
bandstructkpoints2q, dos_properties

#=
Methods to find the susceptibility of materials and their logitudinal dielectric response (real and imaginary).
A variety of different methods are provided for ease of cross checking results. For instance, one may opt to find the real part of the susceptibility 
either through kramers kronig of the imaginary susceptibility or directly from 2d integration.
For Kramers-Kronig, we've also provided methods to find the imaginary susceptibility from the real susceptibility. 
This may be used to ensure results are reliable. 
=#
include("susceptibility_from_wannier.jl")
export direct_plasmon, return_2d_epsilons, return_2d_conductivity, confinement, ϵ, ImΠ

include("kramers_kronig.jl")
export im_polarization, kramers_kronig, kramers_kronig_scipy, kramers_kronig_quadgk, im_polarization_cubature, 
return_2d_epsilon, return_2d_epsilon_scipy, direct_epsilon,
direct_epsilon_cubature, return_2d_epsilon_quadgk,
kramers_kronig_reverse_scipy, kramers_kronig_reverse_quadgk, kramers_kronig_reverse, density_of_states_wannier_quad_check,
im_polarization_finite_temperature, im_polarization_mixedmesh

#=
Methods to plot band structures- either from direct DFT data or from wannier tight binding data 
=#
include("band_structures.jl")
export wannier_bands, wannier_vectors, plot_bands, plotmanybands, hwannier, plotwannierbands, plotbandsoverlayedwannier, 
label_plots

include("analytic_models.jl")
export levitov_kramers_kronig_epsilon, levitov_epsilon, levitov_im_polarization, levitov_integrand, levitov_energy,
graphene_bilayer_plasmon_modes, find_graphene_bilayer_plasmon_modes, find_graphene_plasmon, graphene_total_polarization, 
graphene_total_impolarization,
alevitov, Klevitov

include("export_wannier_hamiltonians.jl")
export export_momentum, export_hwannier, export_heph

#=
smoothing functions- useful for kramers kronig calculations for which a smooth imaginary susceptibility is preferable 
for reliable numerics. 
=#
include("smooth.jl")
export smooth

include("matrix_elements.jl")
export phmatrixelements, pwannier, momentum_matrix_elements, eph_matrix_elements, momentum_from_bloch, hephwannier

#=
Methods to write DFT input files. Note that these input files are specifically written with JDFTX in mind. 
=#
include("./write_input/write_scf.jl")
export write_scf, write_nscf, write_wannier, write_ionpos, write_lattice, write_phonon, write_kpoints, write_randkpoints

#=
Methods to calculate damping of plasmons up to second order in phonon interactions 
=#
include("./loss_calculations/plasmon_losses.jl")
export landau_damping, first_order_damping, second_order_damping

include("./read_prepared_data/read_prepared_data.jl")
export dft_graphene_dos_per_area, dft_graphene_phonon_dispersion, graphene_dos_check, graphene_wannier_impolarization

include("3d_methods.jl")
export im_epsilon_3d

include("non_wannier_methods.jl")
export nonwannier3dimepsilon, nonwannierimpol

include("band_projections.jl")
export bandprojections, plotbandprojections

include("density_plotting.jl")
export plot_density, plot_diffdensity, plot_wfns, wavefunctionoverlap, plot_dtot, plot_scfpotential

include("spectral_functions.jl")
export eliashberg, vFsquaredatmu, subsampling, eliashberg2, subsampling2, eliashbergresistivity, eliashbergresistivities

include("spectralfunctions2.jl")
export eliashberg3, returnfermikpoint, eliashberg4, eliashberg5

include("./loss_calculations/generalized_plasmon_losses.jl")
export landau_damping

include("./loss_calculations/marinko.jl")
export pri, ReS, ImS

include("heat_capacities.jl")
export lattice_heatcapacity, electron_heatcapacities, electron_heatcapacity

include("ephcoupling.jl")
export ephcoupling, ephcoupling2

include("loss_calculations/analyticlossmodels.jl")
export graphenetwoplasmonemission

include("loss_calculations/phononlosses.jl")
export forderphononloss

include("boltzmann.jl")
export drude_conductivity, interbandsigma, τ

include("transverse_modes.jl")
export te_plasmon

include("parse_output.jl")
export get_d, get_mag, list_energy, get_force, parse_wannier_band_ranges


include("wfns.jl")
export gvectors, return_cg
end # module

