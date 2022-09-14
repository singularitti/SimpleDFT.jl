mutable struct SCF
    # SCF struct that holds all calculation intermediate results.
    atoms::Atoms
    pot::Array{ComplexF64}
    W::Matrix{ComplexF64}
    Y::Matrix{ComplexF64}
    n::Array{ComplexF64}
    phi::Array{ComplexF64}
    exc::Array{ComplexF64}
    vxc::Array{ComplexF64}
    Eewald::Float64

    function SCF(atoms, pot, W)
        #= Initialize SCF struct. =#
        M = Matrix{ComplexF64}(undef, 0, 0)
        A = ComplexF64[]
        new(atoms, pot, W, M, A, A, A, A, 0.0)
    end
end


function runSCF(atoms::Atoms)
    #= SCF function to handle direct minimizations. =#
    pot = coulomb(atoms)
    W = _init_W(atoms)
    scf = SCF(atoms, pot, W)
    return run(scf)
end


function run(scf::SCF; Nit::Int64=1001, etol::Float64=1e-6)
    #= Run the self-consistent field (SCF) calculation. =#
    scf.Eewald = get_Eewald(scf.atoms)
    Etot = sd(scf, Nit; etol=etol)
    return Etot
end


function _init_W(atoms::Atoms; rand_seed::Int64=1234)
    #= Generate random initial-guess coefficients as starting values.
    Thesis: List. 3.18
    =#
    W = pseudo_uniform((length(atoms.G2c), atoms.Nstate); seed=rand_seed)
    return orth(atoms, W)
end
