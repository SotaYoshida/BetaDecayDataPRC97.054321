# BetaDecayData_PRC97.054321

Complete data for [Phys. Rev. C 97, 054321 (2018)](https://link.aps.org/doi/10.1103/PhysRevC.97.054321), Erratum [Phys. Rev. C 109, 029904](https://doi.org/10.1103/PhysRevC.109.029904).

- `halflifePn.csv` Halflives and beta-delayed one neutron emission probabilities
- `all_transit_FF.csv` First-Forbidden transitions (for nuclei considered in the paper)
- `all_transit_GT.csv` Gamow-Teller transitions (for nuclei considered in the paper)


For some nuclei, we consider multiple candidates for total J of the parent ground state.  
Please make sure which one you want to use.

**2024/11/19**

- added results of K -> Ca decays
- updated FF contributions on Mg -> Al decays

**2023/12/03**

- added results of Mg isotopes
- `beta_Mg_script.jl` is a sample script for [NuclearToolkit.jl](https://github.com/SotaYoshida/NuclearToolkit.jl) (v 0.3.5) to calculate half-lives, Pn, etc. from KSHELL logfiles.

