<p align="center"><img src="assets/banner.svg" alt="The Sphere Packing Problem in Dimension 4 and the Twenty-Four-Cell Conjecture" width="100%"></p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.22766562"><img src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22766562-22d3ee?style=for-the-badge" alt="DOI"></a>
  <a href="papers/D4Paper.pdf"><img src="https://img.shields.io/badge/paper-PDF-f472b6?style=for-the-badge" alt="paper"></a>
  <img src="https://img.shields.io/badge/Lean_4-219_theorems-7c5cff?style=for-the-badge" alt="Lean 4">
  <img src="https://img.shields.io/badge/Julia-557_checks-34d399?style=for-the-badge" alt="Julia">
  <img src="https://img.shields.io/badge/license-MIT-f59e0b?style=for-the-badge" alt="MIT">
</p>

**Deep Bhattacharjee**, Ushashi Bhattacharya, Priyabrata Mandal (corresponding author), Shounak Bhattacharya.

| Proved | Open |
| :--- | :--- |
| Every packing of unit balls in R⁴ has density below **0.630996** | Optimality of D₄ (density π²/16) |
| Voronoi cell volume ≥ 8, equality only for the 24-cell, for all contact configurations | The twenty-four-cell conjecture in general |
| The same for at most 23 or at least 28 centres within √6, and for the proved radial profiles of 27 | (C): 25 and 26 centres, the remaining profiles of 27 |
| The same for 24 centres whose directions lie within 1/48 of the root system | (G): 24 centres whose directions are farther than 1/48 from the root system |

| Folder | Contents | Run |
| :--- | :--- | :--- |
| [`papers/`](papers) | the paper | |
| [`lean/`](lean) | Lean 4 checks of the finite arithmetic, no Mathlib | `sh lean/run_all.sh` |
| [`julia/`](julia) | second implementation of the certificate checks | `cd julia && julia --project=. run_all.jl --tier standard` |
