"""
Overlap operator.
Thesis: List. 3.9
"""
function op_O(atoms::Atoms, W::Matrix{ComplexF64})
    return atoms.Omega .* W
end


"""
Laplacian operator.
Thesis: Eq. 3.10
List. 3.11
"""
function op_L(atoms::Atoms, W::Matrix{ComplexF64})
    if size(W)[1] == length(atoms.G2c)
        G2 = atoms.G2c
    else
        G2 = atoms.G2
    end
    return -atoms.Omega .* G2 .* W
end


"""
Inverse Laplacian operator.
Thesis: List. 3.12
"""
function op_Linv(atoms::Atoms, W::Array{ComplexF64})
    out = W ./ atoms.G2 ./ -atoms.Omega
    out[1] = 0.0
    return out
end


"""
Inverse Laplacian operator.
Thesis: List. 3.12
"""
function op_Linv(atoms::Atoms, W::Matrix{ComplexF64})
    out = W ./ atoms.G2 ./ -atoms.Omega
    out[1, :] .= 0.0
    return out
end


"""
Backwards transformation from reciprocal space to real-space.
Thesis: Eq. 3.11
List. 3.13
"""
function op_I(atoms::Atoms, W::Array{ComplexF64})
    n = prod(atoms.s)
    if size(W)[1] == length(atoms.G2)
        Wfft = W
    else
        Wfft = convert(Array{ComplexF64}, atoms.active)
        Wfft[Wfft.!=0] = W
    end
    Wfft = reshape(Wfft, atoms.s[1], atoms.s[2], atoms.s[3])
    Finv = vec(ifft(Wfft))
    return Finv .* n
end


"""
Backwards transformation from reciprocal space to real-space.
Thesis: Eq. 3.11
List. 3.13
"""
function op_I(atoms::Atoms, W::Matrix{ComplexF64})
    Finv = zeros(ComplexF64, length(atoms.G2), size(W)[2])
    for i = 1:size(W)[2]
        Finv[:, i] = op_I(atoms, W[:, i])
    end
    return Finv
end


"""
Forward transformation from real-space to reciprocal space.
Thesis: Eq. 3.12
List. 3.14
"""
function op_J(atoms::Atoms, W::Array{ComplexF64})
    n = prod(atoms.s)
    Wfft = reshape(W, atoms.s[1], atoms.s[2], atoms.s[3])
    Finv = vec(fft(Wfft))
    return Finv ./ n
end


"""
Backwards transformation from reciprocal space to real-space.
Thesis: Eq. 3.11
List. 3.13
"""
function op_J(atoms::Atoms, W::Matrix{ComplexF64})
    Finv = zeros(ComplexF64, length(atoms.G2), size(W)[2])
    for i = 1:size(W)[2]
        Finv[:, i] = op_J(atoms, W[:, i])
    end
    return Finv
end


"""
Conjugated backwards transformation from real-space to reciprocal space.
Thesis: Eq. 3.13
List. 3.15
"""
function op_Idag(atoms::Atoms, W::Array{ComplexF64})
    n = prod(atoms.s)
    F = op_J(atoms, W)
    F = F[atoms.active[:, 1]]
    return F .* n
end


"""
Conjugated backwards transformation from real-space to reciprocal space.
Thesis: Eq. 3.13
List. 3.15
"""
function op_Idag(atoms::Atoms, W::Matrix{ComplexF64})
    F = zeros(ComplexF64, length(atoms.G2c), size(W)[2])
    for i = 1:size(W)[2]
        F[:, i] = op_Idag(atoms, W[:, i])
    end
    return F
end


"""
Conjugated forward transformation from reciprocal space to real-space.
Thesis: Eq. 3.14
List. 3.16
"""
function op_Jdag(atoms::Atoms, W::Array{ComplexF64})
    n = prod(atoms.s)
    Finv = op_I(atoms, W)
    return Finv ./ n
end


"""
Conjugated forward transformation from reciprocal space to real-space.
Thesis: Eq. 3.14
List. 3.16
"""
function op_Jdag(atoms::Atoms, W::Matrix{ComplexF64})
    n = prod(atoms.s)
    Finv = op_I(atoms, W)
    return Finv ./ n
end
