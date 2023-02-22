![simpledft.jl logo](https://gitlab.com/wangenau/simpledft.jl/-/raw/main/logo/simpledft_jl_logo.png)

# SimpleDFT.jl
[![language](https://img.shields.io/badge/language-Julia-green)](https://www.python.org)
[![license](https://img.shields.io/badge/license-APACHE2-lightgrey)](https://gitlab.com/wangenau/simpledft.jl/-/blob/main/LICENSE)

SimpleDFT.jl is a simple plane wave density functional theory (DFT) code.
It is a Julia implementation of the [DFT++](https://arxiv.org/abs/cond-mat/9909130) pragmas proposed by Thomas Arias et al.
Also, it serves as the minimalistic prototype for the [eminus](https://gitlab.com/wangenau/eminus) code,
which was introduced in the [master thesis](https://www.researchgate.net/publication/356537762_Domain-averaged_Fermi_holes_A_self-interaction_correction_perspective) of Wanja Timm Schulze to explain theory and software development compactly.
This version is a straightforward translation of the [SimpleDFT](https://gitlab.com/wangenau/simpledft) code from Python to Julia.
The resulting energy difference between all codes is well below the mEₕ range.

| SimpleDFT.jl | Description |
| --------- | ----------- |
| Language | Julia 1.8 |
| License | Apache 2.0 |
| Dependencies | Only FFTW |
| Basis set| Plane waves (PW) |
| DFT | Restricted Kohn-Sham (RKS) |

# Installation
All necessary dependencies can be installed in your Julia REPL using

```terminal
using Pkg
Pkg.add("FFTW")
```

# Examples
Example calculations, i.e., the H atom, He atom, and H2 molecule can be executed with

```terminal
julia examples.jl
```

# Simplifications
This code is about implementing plane wave DFT as simple as possible, while still being general.
To classify the shortcomings from a fully-fletched DFT code, the most important simplifications are listed below
* Restricted Kohn-Sham DFT only, no open-shell systems
* LDA only, no functionals using gradients
* All-electron Coulomb potential, no pseudopotentials to only treat valence electrons
* Gamma-point only, no band paths
* Steepest descent only, no sophisticated minimizations
* Random starting guess only, no calculated guess of the density
