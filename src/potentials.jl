function coulomb(atoms::Atoms)
    #= All-electron Coulomb potential.
    Thesis: Eq. 3.15 ff.
            List. 3.17
    =#
    Vcoul = -4 .* pi .* atoms.Z[1] ./ atoms.G2
    Vcoul[1] = 0
    return op_J(atoms, Vcoul .* atoms.Sf)
end
