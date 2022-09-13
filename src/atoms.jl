mutable struct Atoms
    # Atoms struct that holds all system and cell parameters.
    atom::Array{String}
    X::Matrix{Float64}
    a::Float64
    ecut::Float64
    Z::Array{Float64}
    s::Array{Int64}
    f::Array{Float64}
    Natoms::Int64
    Nstate::Int64
    R::Matrix{Float64}
    Omega::Float64
    r::Matrix{Float64}
    G::Matrix{Float64}
    G2::Array{Float64}
    active::Array{Bool}
    G2c::Array{Float64}
    Sf::Array{ComplexF64}

    function Atoms(atom, X, a, ecut, Z, s, f)
        #= Initialize and build all necessary parameters. =#
        M, N = _get_index_matrices(s)
        Natoms, Nstate, R, Omega, r = _set_cell(atom, a, s, f, M)
        G, G2, active, G2c, Sf = _set_G(X, ecut, R, N)
        new(atom, X, a, ecut, Z, s, f, Natoms, Nstate, R, Omega, r, G, G2, active, G2c, Sf)
    end
end


function _get_index_matrices(s::Array{Int64})
    #= Build index matrices M and N to build the real and reciprocal space samplings.
    Thesis: List 3.4
            List 3.5
    =#
    ms = 0:prod(s) - 1
    m1 = ms .% s[1]
    m2 = floor.(Int64, ms ./ s[1]) .% s[2]
    m3 = floor.(Int64, ms ./ (s[1] .* s[2])) .% s[3]
    M = [m1 m2 m3]

    n1 = m1 .- (m1 .> s[1] ./ 2) .* s[1]
    n2 = m2 .- (m2 .> s[2] ./ 2) .* s[2]
    n3 = m3 .- (m3 .> s[3] ./ 2) .* s[3]
    N = [n1 n2 n3]
    return M, N
end


function _set_cell(atom::Array{String}, a::Float64, s::Array{Int64}, f::Array{Float64}, M::Matrix{Int64})
    #= Build the unit cell and create the respective sampling.
    Thesis: Eq. 3.3
            List. 3.3
            Eq. 3.5
            List. 3.3
    =#
    Natoms = length(atom)
    Nstate = length(f)

    R = a .* Matrix(1.0I, 3, 3)
    Omega = abs(det(R))
    r = M * inv(Diagonal(s)) * R'
    return Natoms, Nstate, R, Omega, r
end


function _set_G(X::Matrix{Float64}, ecut::Float64, R::Matrix{Float64}, N::Matrix{Int64})
    #= Build G-vectors, build squared magnitudes G2, and generate the active space.
    Thesis: Eq. 3.8
            List. 3.5
            List. 3.6
            List. 3.7
            Eq. 3.9
            List. 3.8
    =#
    G = 2 .* pi .* N * inv(R)
    G2 = sum(G.^2; dims=2)
    active = G2 .<= 2 * ecut
    G2c = G2[active]
    Sf = sum(exp.(-1im .* G * X'), dims=2)
    return G, G2, active, G2c, Sf
end