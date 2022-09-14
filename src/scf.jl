"SCF struct that holds all calculation intermediate results."
mutable struct SCF
    atoms::Atoms
    pot::Array{ComplexF64}
    W::Matrix{ComplexF64}
    Y::Matrix{ComplexF64}
    n::Array{ComplexF64}
    phi::Array{ComplexF64}
    exc::Array{ComplexF64}
    vxc::Array{ComplexF64}
    Eewald::Float64

    "Initialize SCF struct."
    function SCF(atoms, pot, W)
        M = Matrix{ComplexF64}(undef, 0, 0)
        A = ComplexF64[]
        new(atoms, pot, W, M, A, A, A, A, 0.0)
    end
end


"SCF function to handle direct minimizations."
function runSCF(atoms::Atoms; Nit::Int64=1001, etol::Float64=1e-6)
    pot = coulomb(atoms)
    W = init_W(atoms)
    scf = SCF(atoms, pot, W)
    return run(scf; Nit=Nit, etol=etol)
end


"Run the self-consistent field (SCF) calculation."
function run(scf::SCF; Nit::Int64=1001, etol::Float64=1e-6)
    scf.Eewald = get_Eewald(scf.atoms)
    Etot = sd(scf, Nit; etol=etol)
    return Etot
end


"""
Generate random initial-guess coefficients as starting values.
Thesis: List. 3.18
"""
function init_W(atoms::Atoms; rand_seed::Int64=1234)
    W = pseudo_uniform((length(atoms.G2c), atoms.Nstate); seed=rand_seed)
    return orth(atoms, W)
end
