"Perform one SCF step for a DFT calculation."
function scf_step(scf::SCF)
    scf.Y = orth(scf.atoms, scf.W)
    scf.n = get_n_total(scf.atoms, scf.Y)
    scf.phi = solve_poisson(scf.atoms, scf.n)
    x, c = lda_slater_x(scf.n), lda_chachiyo_c(scf.n)
    scf.exc = x[1] + c[1]
    scf.vxc = x[2] + c[2]
    return get_E(scf)
end


"""
Steepest descent minimization algorithm.
Thesis: List. 3.21
Fig. 3.2
"""
function sd(scf::SCF, Nit::Int64; etol::Float64 = 1e-6, beta::Float64 = 1e-5)
    Elist = Float64[]

    for i = 1:Nit
        E = scf_step(scf)
        append!(Elist, E)
        print("Nit: $(i)  \tEtot: $(round(E; digits=6)) Eh\r")
        if i > 1 && abs(Elist[i-1] - Elist[i]) < etol
            println("\nSCF converged.")
            return E
        end
        g = get_grad(scf.atoms, scf.W, scf.Y, scf.n, scf.phi, scf.vxc, scf.pot)
        scf.W = scf.W .- beta .* g
    end
    println("\nSCF not converged!")
    return E
end
