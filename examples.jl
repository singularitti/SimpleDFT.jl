#=
Example calculations for the SimpleDFT.jl code, including H, He, and H2.
References for the DFT++ formulation: https://doi.org/10.1016/s0010-4655(00)00072-2
Code annotations reference this master thesis: http://dx.doi.org/10.13140/RG.2.2.27958.42568/2
More documentation can be found inside the eminus routines: https://gitlab.com/esp42/sage/eminus
=#
include("src/simpledft.jl")


function calculate(atoms)
    etot = runSCF(atoms)
    println("Etot($(atoms.atom)) = $(round(etot; digits=6)) Eh")
end


H_atom = Atoms(["H"], [0.0 0.0 0.0;], 16.0, 16.0, [1.0], [60, 60, 60], [1.0])
calculate(H_atom)
# Output:  Etot(["H"]) = -0.438419 Eh

He_atom = Atoms(["He"], [0.0 0.0 0.0;], 16.0, 16.0, [2.0], [60, 60, 60], [2.0])
calculate(He_atom)
# Output:  Etot(["He"]) = -2.632034 Eh

# Experimental geometry from CCCBDB: https://cccbdb.nist.gov/exp2x.asp?casno=1333740&charge=0
H2_atom = Atoms(["H", "H"], [0.0 0.0 0.0; 1.4 0.0 0.0], 16.0, 16.0, [1.0 1.0], [60, 60, 60], [2.0])
calculate(H2_atom)
# Output:  Etot(["H", "H"]) = -1.113967 Eh
