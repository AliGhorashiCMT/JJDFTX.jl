# JJDFTX.jl
[![Build status][ci-status-img]][ci-status-url][![Coverage][codecov-img]][codecov-url]

JJDFTX.jl serves to bridge the gap between computational quantum chemistry (e.g. DFT) and photonics. The point of the package is to make easy the calculation of relevant physical quantities of materials that determine photonic properties such as plasmonic losses due to the electron-phonon interaction. Below, we show just a sample of the functionality of the package. 

Some basic plotting functionality demonstrated with hexagonal boron nitride: 

![hbn]

Spin densities may also be visualized and analyzed in the case of defect systems: 

![spindensity]

Phonon dispersions may be plotted using the phonon_force_matrix and the phonon_dispersion_path methods. 
An example below is of the monolayer graphene phonon dispersion along the Gamma-M-K-Gamma path.

![MLGPhonons]

Monolayer graphene plasmon dispersion: 

![mlgplas]

Two plasmon absorption in graphene as detailed in Jablan, Marinko, and Darrick E. Chang. "Multiplasmon absorption in graphene." Physical review letters 114.23 (2015): 236801.

![plasabs]

Phonon-assisted plasmon loss in graphene for a chemical potential of 0.135 eV calculated through a Fermi's golden rule. 
An analogous calculation is done in Jablan, Marinko, Hrvoje Buljan, and Marin Soljačić. "Plasmonics in graphene at infrared frequencies." Physical review B 80.24 (2009): 245435.

![phloss]

Twisted bilayer graphene plasmon dispersion:

![tbgplas]

Two plasmon modes in spatially separated graphene planes (50 nanometer separation at 0.4 eV doping. Units are always in terms of the chemical potential and the Fermi wavevector)

![twomodes]

Electron-Phonon Coupling in Aluminum (which may be cross checked with the reference Brown, Ana M., et al. "Ab initio phonon coupling and optical response of hot electrons in plasmonic metals." Physical Review B 94.7 (2016): 075120. )

![eph]

[MLGPhonon]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/MLGPhonon.png
[spindensity]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/SpinDensity.png
[hbn]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/hBNdensity.png
[eph]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/EphAl.png
[twomodes]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/TwoModes.png
[phloss]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/PhononPlasmon0135.png
[plasabs]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/2PAbs.png
[mlgplas]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/MLGPlasmon.png
[tbgplas]: https://github.com/AliGhorashiCMT/JJDFTX.jl/blob/main/imgs/tbg_graphene.png
[ci-status-img]:   https://github.com/AliGhorashiCMT/JJDFTX.jl/workflows/CI/badge.svg
[ci-status-url]:   https://github.com/AliGhorashiCMT/JJDFTX.jl/actions
[codecov-img]: https://codecov.io/gh/AliGhorashiCMT/JJDFTX.jl/branch/main/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/AliGhorashiCMT/JJDFTX.jl

