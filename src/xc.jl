"""
Slater exchange functional (spin-paired).
Thesis: Eq. 2.11
"""
function lda_x(n::Matrix{ComplexF64})
    f = -3.0 / 4.0 * (3.0 / (2.0 * pi))^(2.0 / 3.0)
    rs = (3.0 ./ (4.0 .* pi .* n)) .^ (1.0 / 3.0)

    ex = f ./ rs
    vx = 4.0 ./ 3.0 .* ex
    return ex, vx
end


"Chachiyo parametrization of the correlation functional (spin-paired)."
function lda_c_chachiyo(n::Matrix{ComplexF64})
    a = -0.01554535  # (log(2.0) - 1.0) / (2.0 * pi^2.0)
    b = 20.4562557

    rs = (3.0 ./ (4.0 .* pi .* n)) .^ (1.0 / 3.0)

    ec = a .* log.(1.0 .+ b ./ rs .+ b ./ rs .^ 2.0)
    vc = ec .+ a .* b .* (2.0 .+ rs) ./ (3.0 .* (b .+ b .* rs .+ rs .^ 2.0))
    return ec, vc
end


"""
Vosko-Wilk-Nusair parametrization of the correlation functional (spin-paired).
Not used, only for reference as it was used in the master thesis.
Thesis: Eq. 2.12 ff.
"""
function lda_c_vwn(n::Matrix{ComplexF64})
    A = 0.0310907
    b = 3.72744
    c = 12.9352
    x0 = -0.10498

    rs = (3.0 ./ (4.0 .* pi .* n)) .^ (1.0 / 3.0)

    x = sqrt.(rs)
    X = rs .+ b .* x .+ c
    Q = sqrt(4.0 * c - b^2.0)
    fx0 = b * x0 / (x0^2 + b * x0 + c)
    f3 = 2.0 * (2.0 * x0 + b) / Q
    tx = 2.0 .* x .+ b
    tanx = atan.(Q ./ tx)

    ec = A .* (log.(rs ./ X) .+ 2.0 .* b ./ Q .* tanx .- fx0 .* (log.((x .- x0) .^ 2.0 ./ X) .+ f3 .* tanx))

    tt = tx .^ 2.0 .+ Q .^ 2
    vc = ec .- x .* A ./ 6.0 .* (2.0 ./ x .- tx ./ X .- 4.0 .* b ./ tt .- fx0 .* (2.0 ./ (x .- x0) .- tx ./ X .- 4.0 .* (2.0 * x0 + b) ./ tt))
    return ec, vc
end
