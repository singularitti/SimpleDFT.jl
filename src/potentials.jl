"""
All-electron Coulomb potential.
Thesis: Eq. 3.15 ff.
List. 3.17
"""
function coulomb(atoms::Atoms)
    Vcoul = -4.0 .* pi .* atoms.Z[1] ./ atoms.G2
    Vcoul[1] = 0.0
    return op_J(atoms, Vcoul .* atoms.Sf)
end
