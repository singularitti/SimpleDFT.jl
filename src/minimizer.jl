function scf_step(scf::SCF)
    #= Perform one SCF step for a DFT calculation. =#
    scf.Y = orth(scf.atoms, scf.W)
    scf.n = get_n_total(scf.atoms, scf.Y)
    scf.phi = solve_poisson(scf.atoms, scf.n)
    scf.exc = lda_slater_x(scf.n)[1] + lda_chachiyo_c(scf.n)[1]
    scf.vxc = lda_slater_x(scf.n)[2] + lda_chachiyo_c(scf.n)[2]
    return get_E(scf)
end


function sd(scf::SCF, Nit::Int64; etol::Float64=1e-6, beta::Float64=1e-5)
    #= Steepest descent minimization algorithm.
    Thesis: List. 3.21
            Fig. 3.2
    =#
    Elist = Float64[]

    for i in 1:Nit
        E = scf_step(scf)
        append!(Elist, E)
        print("Nit: $(i)  \tEtot: $(round(E; digits=6)) Eh\r")
        if i > 1 && abs(Elist[i - 1] - Elist[i]) < etol
            println("\nSCF converged.")
            return E
        end
        g = get_grad(scf.atoms, scf.W, scf.Y, scf.n, scf.phi, scf.vxc, scf.pot)
        scf.W = scf.W .- beta .* g
    end
    println("\nSCF not converged!")
    return E
end
