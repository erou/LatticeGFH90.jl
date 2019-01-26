#using LatticeGFH90
include("src/LatticeGFH90.jl")

function simple_bench(p, l)
    k, x = FiniteField(p, l, "x")
    @time A = tensor_algebra(k)
    @time h = solve_h90(A)
    return level(A)
end
