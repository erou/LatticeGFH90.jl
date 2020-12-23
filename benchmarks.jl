using Nemo: FiniteField
using LatticeGFH90
using LatticeGFH90: level, ZETAS, H90_ELEMENTS, tensor_algebra, solve_h90,
                    degree

path_nemo = dirname(dirname(pathof(LatticeGFH90)))
path_benchmark = joinpath(path_nemo, "benchmarks")

function simple_bench(p, l)
    make_zetas_conway(p)
    k, x = FiniteField(p, l, "x")
    @time A = tensor_algebra(k)
    @time h = solve_h90(A)
    return level(A)
end

function bench_solve_h90(p, N, s = "solve_h90")

    make_zetas_conway(p)
    path = joinpath(path_benchmark, string(s,"-", p, ".txt"))

    io = open(path, "w+")

    for j in 1:N
        if j%p == 0
            continue
        end
        a = level(j, p)
        if haskey(ZETAS, (p, a))
            println(j)
            k, x = FiniteField(p, j, "x")
            res = @timed tensor_algebra(k)
            A, t = res[1], res[2]
            t += @elapsed solve_h90(A)
            write(io, string(level(A), ",", degree(A), ",", t, "\n"))
        end
    end

    close(io)
end

function bench_embed(p, N, s = "embed")

    make_zetas_conway(p)
    path = joinpath(path_benchmark, string(s, "-", p, ".txt"))

    io = open(path, "w+")

    # Compute the H90 sol in A_l with l <= N
    for j in 1:N
        if j%p == 0
            continue
        end
        a = level(j, p)
        if haskey(ZETAS, (p, a))
            k, x = FiniteField(p, j, "x")
            A = tensor_algebra(k)
            h = solve_h90(A)
            H90_ELEMENTS[(p, j)] = h
        end
    end

    # Sort the A_l
    degrees = sort!(Int[y for (x,y) in keys(H90_ELEMENTS)])
    len = length(degrees)

    # Make all the possible embeddings
    for j in 1:len
        for i in (j+1):len
            l, m = degrees[j], degrees[i]
            if m % l == 0
                println(string(l, " ... ", m))
                k, x = FiniteField(p, l, "x")
                K, y = FiniteField(p, m, "y")
                a, b = level(l, p), level(m, p)

                t = @elapsed embed(k, K)
                write(io, string(a, ",", l, ",", b, ",", m, ",", t, "\n"))
            end
        end
    end
    close(io)
end

function bench_all(p, N)

    make_zetas_conway(p)
    path = joinpath(path_benchmark, string("solve-h90-", p, ".txt"))
    io = open(path, "w+")

    path2 = joinpath(path_benchmark, string("embed-", p, ".txt"))
    io2 = open(path2, "w+")

    # Compute the H90 sol in A_l with l <= N
    for j in 1:N
        if j%p == 0
            continue
        end
        a = level(j, p)
        if haskey(ZETAS, (p, a))

            println(j)
            k, x = FiniteField(p, j, "x")
            res1 = @timed tensor_algebra(k)
            A, t = res1[1], res1[2]

            res2 = @timed solve_h90(A)
            h, u = res2[1], res2[2] 
            write(io, string(level(A), ",", degree(A), ",", t+u, "\n"))
            H90_ELEMENTS[(p, j)] = h
        end
    end

    # Sort the A_l
    degrees = sort!(Int[y for (x,y) in keys(H90_ELEMENTS) if x == p])
    len = length(degrees)

    # Make all the possible embeddings
    for j in 1:len
        for i in (j+1):len
            l, m = degrees[j], degrees[i]
            if m % l == 0
                println(string(l, " ... ", m))
                k, x = FiniteField(p, l, "x")
                K, y = FiniteField(p, m, "y")
                a, b = level(l, p), level(m, p)

                v = @elapsed embed(k, K)
                write(io2, string(a, ",", l, ",", b, ",", m, ",", v, "\n"))
            end
        end
    end
    close(io)
    close(io2)
end
