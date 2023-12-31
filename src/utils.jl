"Lehmer random number generator, following MINSTD."
function pseudo_uniform(size::Tuple{Int64,Int64}; seed::Int64 = 1234)
    U = zeros(ComplexF64, size)
    mult = 48271
    mod = (2^31) - 1
    x = (seed * mult + 1) % mod
    for i = 1:size[1]
        for j = 1:size[2]
            x = (x * mult + 1) % mod
            U[i, j] = x / mod
        end
    end
    return U
end
