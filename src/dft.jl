"""
Solve the Poisson equation.
Thesis: Eq. 2.48
"""
function solve_poisson(atoms::Atoms, n::Array{ComplexF64})
    return -4.0 .* pi .* op_Linv(atoms, op_O(atoms, op_J(atoms, n)))
end


"""
Calculate the total electronic density.
Thesis: Eq. 2.36
List. 3.23
"""
function get_n_total(atoms::Atoms, Y::Matrix{ComplexF64})
    Yrs = op_I(atoms, Y)
    n = real(conj(Yrs) .* atoms.f .* Yrs)
    return sum(n, dims = 2)
end


"""
Orthogonalize coefficient matrix W.
Thesis: Eq. 2.34 ff.
"""
function orth(atoms::Atoms, W::Matrix{ComplexF64})
    U = sqrt(W' * op_O(atoms, W))
    return W * inv(U)
end


"""
Calculate the energy gradient with respect to W.
Thesis: Eq. 2.43
List. 3.24
"""
function get_grad(atoms::Atoms, W::Matrix{ComplexF64}, Y::Matrix{ComplexF64}, n::Array{ComplexF64}, phi::Array{ComplexF64}, vxc::Array{ComplexF64}, Vreciproc::Matrix{ComplexF64})
    F = Diagonal(atoms.f)
    HW = H(atoms, W, Y, n, phi, vxc, Vreciproc)
    WHW = W' * HW
    OW = op_O(atoms, W)
    U = W' * OW
    invU = inv(U)
    U12 = sqrt(invU)
    Ht = U12 * WHW * U12
    return (HW .- (OW * invU) * WHW) * (U12 * F * U12) .+ OW * (U12 * Q(Ht * F .- F * Ht, U))
end


"""
Left-hand side of the eigenvalue equation.
Thesis: Eq. 2.45 ff.
List. 3.26
"""
function H(atoms::Atoms, W::Matrix{ComplexF64}, Y::Matrix{ComplexF64}, n::Array{ComplexF64}, phi::Array{ComplexF64}, vxc::Array{ComplexF64}, Vreciproc::Matrix{ComplexF64})
    Veff = Vreciproc .+ op_Jdag(atoms, op_O(atoms, op_J(atoms, vxc) .+ phi))
    return -0.5 .* op_L(atoms, W) .+ op_Idag(atoms, Veff .* op_I(atoms, W))
end


"""
Operator needed to calculate gradients with non-constant occupations.
Thesis: Eq. 2.47
List. 3.25
"""
function Q(inp::Matrix{ComplexF64}, U::Matrix{ComplexF64})
    F = eigen(U)
    mu, V = F.values, F.vectors
    denom = sqrt.(mu) * ones(1, length(mu))
    denom2 = denom .+ denom'
    return V * ((V' * inp * V) ./ denom2) * V'
end
