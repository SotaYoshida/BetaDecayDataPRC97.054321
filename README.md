# BetaDecayData_PRC97.054321

Complete data for [Phys. Rev. C 97, 054321 (2018)](https://link.aps.org/doi/10.1103/PhysRevC.97.054321), Erratum [Phys. Rev. C 109, 029904](https://doi.org/10.1103/PhysRevC.109.029904).

- `halflifePn.csv` Halflives and beta-delayed one neutron emission probabilities
- `all_transit_FF.csv` First-Forbidden transitions (for nuclei considered in the paper)
- `all_transit_GT.csv` Gamow-Teller transitions (for nuclei considered in the paper)


For some nuclei, we consider multiple candidates for total J of the parent ground state.  
Please make sure which one you want to use.

**2025/03/03**

The 45Cl(j3p)=>45Ar case is added (it was missing). It should also be noted that a couple of experimental works suggest that the ground state is to be $J=3/2^+$ instead of $J=1/2^+$ (g.s. predicted by the SDPF(SDG)-MU interaction):
- I. Cox et al., Phys. Rev. Lett. 132, 152503 (2024) [DOI](https://doi.org/10.1103/PhysRevLett.132.152503)
- V. Tripathi et al., Phys. Rev. C 109, 044320 (2024) [DOI](10.1103/PhysRevC.109.044320)

**2024/11/19**
- added results of K -> Ca decays
- updated FF contributions on Mg -> Al decays

**2023/12/03**

- added results of Mg isotopes
- `beta_Mg_script.jl` is a sample script for [NuclearToolkit.jl](https://github.com/SotaYoshida/NuclearToolkit.jl) (v 0.3.5) to calculate half-lives, Pn, etc. from KSHELL logfiles.

