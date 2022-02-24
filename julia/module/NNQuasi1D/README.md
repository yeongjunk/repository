## Description
Quasi 1D nearest neighbor hopping single-electron lattice model(bar and strip) with periodic boundary condition. The purpose of this module is to study the Anderson localization and to compare three different methods(ED, RGF, TMM) for computing localization lengths for this model.

This code is written in Julia-1.6.1.

Prerequisite custom modules: Lattice, PN, Binning  (written by Yeongjun Kim)

contact: yeongjun@ust.ac.kr

 - ED: exact diagonization
 - RGF: recursive green's function
 - TMM: transfer matrix method

## Reference
Mackinnon and Kramer, Z. Phys. B Condensed Matter 53, 1-13 (1983)
