![simpledft.jl logo](https://gitlab.com/esp42/sage/simpledft.jl/-/raw/main/logo/simpledft_jl_logo.png)

# SimpleDFT.jl
[![language](https://img.shields.io/badge/language-Julia-green)](https://www.python.org)
[![license](https://img.shields.io/badge/license-APACHE2-lightgrey)](https://gitlab.com/esp42/sage/simpledft.jl/-/blob/main/LICENSE)

SimpleDFT.jl is a simple plane wave density functional theory (DFT) code.
It is a Julia implementation of the [DFT++](https://arxiv.org/abs/cond-mat/9909130) pragmas proposed by Thomas Arias et al.
This version is a straightforward translation of the [SimpleDFT](https://gitlab.com/esp42/sage/simpledft) code from Python to Julia.

| SimpleDFT.jl | Description |
| --------- | ----------- |
| Language | Julia 1.8 |
| License | Apache 2.0 |
| Dependencies | Only FFTW and SpecialFunctions |
| Basis set| Plane waves (PW) |
| DFT | Restricted Kohn-Sham (RKS) |

# Installation
All necessary dependencies can be installed in your Julia terminal using

```terminal
using Pkg
Pkg.add("FFTW")
Pkg.add("SpecialFunctions")
```

# Examples
Example calculations, i.e., the H atom, He atom, and H2 molecule can be executed with

```terminal
julia examples.jl
```