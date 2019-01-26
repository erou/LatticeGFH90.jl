#using LatticeGFH90
include("src/LatticeGFH90.jl")

function simple_bench(p, l)
    k, x = FiniteField(p, l, "x")
    @time A = tensor_algebra(k)
    @time h = solve_h90(A)
    return level(A)
end

function bench_solve_h90(p, m)

    path = string("../benchmarks/solve_h90-", p, ".txt")

    io = open(path, "w+")

    for j in 1:m
        if j%p == 0
            continue
        end
        a = level(j, p)
        if haskey(ZETAS, (p, a))
            println(j)
            k, x = FiniteField(p, j, "x")
            A, t = (@timed tensor_algebra(k))[1:2]
            t += @elapsed solve_h90(A)
            write(io, string(level(A), ",", degree(A), ",", t, "\n"))
        end
    end

    close(io)
end

function bench_embed(p, m)

    path = string("../benchmarks/embed-", p, ".txt")

    io = open(path, "w+")

    for j in 1:m
        if j%p == 0
            continue
        end
        a = level(2*j, p)
        if haskey(ZETAS, (p, a))
            println(j)
            k, x = FiniteField(p, j, "x")
            A = tensor_algebra(k)
            K, y = FiniteField(p, 2*j, "y")
            B = tensor_algebra(K)
            h = solve_h90(A)
            H90_ELEMENTS[(p, j)] = h
            g = solve_h90(B)
            H90_ELEMENTS[(p, 2*j)] = g
            t = @elapsed embed(k, K)
            write(io, string(level(A), ",", level(B), ",", degree(A), ",", t, "\n"))
        end
    end

    close(io)
end
