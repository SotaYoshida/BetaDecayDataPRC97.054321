# BetaDecayData_PRC97.054321

Complete data for [Phys. Rev. C 97, 054321 (2018)](https://link.aps.org/doi/10.1103/PhysRevC.97.054321).

- `halflifePn.csv` Halflives and beta-delayed one neutron emission probabilities
- `all_transit_FF.csv` First-Forbidden transitions
- `all_transit_GT.csv` Gamow-Teller transitions

For some nuclei, we consider multiple candidates for total J of the parent ground state.  
Please make sure which one you want to use.

**2023/12/03**

- added results of Mg isotopes
- `beta_Mg_script.jl` is a sample script for [NuclearToolkit.jl](https://github.com/SotaYoshida/NuclearToolkit.jl) (v 0.3.5) to calculate half-lives, Pn, etc. from KSHELL logfiles.
