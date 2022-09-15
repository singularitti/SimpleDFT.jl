"""
Slater exchange functional (spin-paired).
Thesis: Eq. 2.11
"""
function lda_slater_x(n::Matrix{ComplexF64}; alpha::Float64 = 2.0 / 3.0)
    rs = (3.0 ./ (4.0 .* pi .* n)) .^ (1.0 / 3.0)

    f = -9.0 / 8.0 * (3.0 / (2.0 * pi))^(2.0 / 3.0)

    ex = f .* alpha ./ rs
    vx = 4.0 ./ 3.0 .* ex
    return ex, vx
end


"Chachiyo parametrization of the correlation functional (spin-paired)."
function lda_chachiyo_c(n::Matrix{ComplexF64})
    rs = (3.0 ./ (4.0 .* pi .* n)) .^ (1.0 / 3.0)

    a = (log(2.0) - 1.0) / (2.0 * pi^2.0)
    b = 20.4562557

    ec = a .* log.(1.0 .+ b ./ rs .+ b ./ rs .^ 2.0)
    vc = ec .+ a .* b .* (2.0 .+ rs) ./ (3.0 .* (b .+ b .* rs .+ rs .^ 2.0))
    return ec, vc
end


"""
Vosko-Wilk-Nusair parametrization of the correlation functional (spin-paired).
Thesis: Eq. 2.12 ff.
"""
function lda_vwn_c(n::Matrix{ComplexF64})
    rs = (3.0 ./ (4.0 .* pi .* n)) .^ (1.0 / 3.0)

    a = 0.0310907
    b = 3.72744
    c = 12.9352
    x0 = -0.10498

    rs12 = sqrt.(rs)
    q = sqrt(4.0 * c - b * b)
    qx = atan.(q ./ (2.0 .* rs12 .+ b))
    f1 = 2.0 * b / q
    f2 = b * x0 / (x0 * x0 + b * x0 + c)
    f3 = 2.0 * (2.0 * x0 + b) / q
    fx = rs .+ b .* rs12 .+ c

    ec = a .* (log.(rs ./ fx) .+ f1 .* qx .- f2 .* (log.((rs12 .- x0) .^ 2.0 ./ fx) .+ f3 .* qx))
    tx = 2.0 .* rs12 .+ b
    tt = tx .* tx .+ q .* q
    vc = ec .- rs12 .* a ./ 6.0 .* (2.0 ./ rs12 .- tx ./ fx .- 4.0 .* b ./ tt .- f2 .* (2.0 ./ (rs12 .- x0) .- tx ./ fx .- 4.0 .* (2.0 .* x0 .+ b) ./ tt))
    return ec, vc
end
