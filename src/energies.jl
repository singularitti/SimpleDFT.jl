function get_E(scf::SCF)
    #= Calculate energy contributions.
    Thesis: Eq. 2.49
    =#
    Ekin = get_Ekin(scf.atoms, scf.Y)
    Ecoul = get_Ecoul(scf.atoms, scf.n, scf.phi)
    Exc = get_Exc(scf.atoms, scf.n, scf.exc)
    Een = get_Een(scf.n, scf.pot)
    return Ekin + Ecoul + Exc + Een + scf.Eewald
end


function get_Ekin(atoms::Atoms, W::Matrix{ComplexF64})
    #= Calculate the kinetic energy.
    Thesis: Eq. 2.37
    =#
    F = Diagonal(atoms.f)
    T = -0.5 * tr(F * (W' * op_L(atoms, W)))
    return real(T)
end


function get_Ecoul(atoms::Atoms, n::Array{ComplexF64}, phi::Array{ComplexF64})
    #= Calculate the Coulomb energy.
    Thesis: Eq. 2.40 + Eq. 2.41 (as in Eq. 2.49)
    =#
    Ecoul = 0.5 * (n' * op_Jdag(atoms, op_O(atoms, phi)))
    return real(Ecoul[1])
end


function get_Exc(atoms::Atoms, n::Array{ComplexF64}, exc::Array{ComplexF64})
    #= Calculate the exchange-correlation energy.
    Thesis: Eq. 2.39
    =#
    Exc = n' * op_Jdag(atoms, op_O(atoms, op_J(atoms, exc)))
    return real(Exc[1])
end


function get_Een(n::Array{ComplexF64}, Vreciproc::Matrix{ComplexF64})
    #= Calculate the electron-ion interaction.
    Thesis: Eq. 2.38
    =#
    Een = Vreciproc' * n
    return real(Een[1])
end


function get_Eewald(atoms::Atoms; gcut::Float64=2.0, gamma::Float64=1e-8)
    #= Calculate the Ewald energy.
    Thesis: Eq. A.12 ff.
    =#
    t1 = atoms.R[:, 1]
    t2 = atoms.R[:, 2]
    t3 = atoms.R[:, 3]
    t1m = sqrt(dot(t1, t1))
    t2m = sqrt(dot(t2, t2))
    t3m = sqrt(dot(t3, t3))

    RecVecs = 2 * pi * inv(atoms.R')
    g1 = RecVecs[:, 1]
    g2 = RecVecs[:, 2]
    g3 = RecVecs[:, 3]
    g1m = sqrt(dot(g1, g1))
    g2m = sqrt(dot(g2, g2))
    g3m = sqrt(dot(g3, g3))

    gexp = -log(gamma)
    nu = 0.5 * sqrt(gcut^2 / gexp)

    x = sum(atoms.Z.^2)
    totalcharge = sum(atoms.Z)

    Eewald = -nu * x / sqrt(pi)
    Eewald += -pi * (totalcharge^2) / (2 * atoms.Omega * nu^2)

    tmax = sqrt(0.5 * gexp) / nu
    mmm1 = round(Int64, tmax / t1m + 1.5)
    mmm2 = round(Int64, tmax / t2m + 1.5)
    mmm3 = round(Int64, tmax / t3m + 1.5)

    for ia = 1:atoms.Natoms
        for ja = 1:atoms.Natoms
            dX = atoms.X[ia, :] .- atoms.X[ja, :]
            ZiZj = atoms.Z[ia] * atoms.Z[ja]
            for i = -mmm1:mmm1
                for j = -mmm2:mmm2
                    for k = -mmm3:mmm3
                        if (ia != ja) || ((abs(i) + abs(j) + abs(k)) != 0)
                            T = i .* t1 .+ j .* t2 .+ k .* t3
                            rmag = sqrt(sum((dX .- T).^2))
                            Eewald += 0.5 * ZiZj * erfc(rmag * nu) / rmag
                        end
                    end
                end
            end
        end
    end

    mmm1 = round(Int64, gcut / g1m + 1.5)
    mmm2 = round(Int64, gcut / g2m + 1.5)
    mmm3 = round(Int64, gcut / g3m + 1.5)

    for ia = 1:atoms.Natoms
        for ja = 1:atoms.Natoms
            dX = atoms.X[ia, :] .- atoms.X[ja, :]
            ZiZj = atoms.Z[ia] * atoms.Z[ja]
            for i = -mmm1:mmm1
                for j = -mmm2:mmm2
                    for k = -mmm3:mmm3
                        if (abs(i) + abs(j) + abs(k)) != 0
                            G = i .* g1 .+ j .* g2 .+ k .* g3
                            GX = sum(G .* dX)
                            G2 = sum(G.^2)
                            x = 2 * pi / atoms.Omega * exp(-0.25 * G2 / nu^2) / G2
                            Eewald += x * ZiZj * cos(GX)
                        end
                    end
                end
            end
        end
    end
    return Eewald
end
