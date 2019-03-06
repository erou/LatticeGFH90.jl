#######################################################################
#
# Basic functions
#
#######################################################################

export make_zetas_prim, make_zetas_conway

const H90_ELEMENTS = Dict{Tuple{Int, Int}, tensor_element}()
const EMBEDDINGS = Dict{Tuple{Int, Int, Int}, Any}() # Any for the "embedding" type
const ZETAS = Dict{Tuple{Int, Int}, nmod_poly}()

"""
    multiplicative_order(x::fq_nmod)

Compute the multiplicative order of `x`.
"""
function multiplicative_order(x::fq_nmod)
    N = Nemo.order(parent(x)) - 1
    order = ZZ(1)
    for (q, a) in Nemo.factor(N)
        tmp = x^divexact(N, q^a)
        while tmp != 1
            order *= q
            tmp = tmp^q
        end
    end
    return order
end

"Return `true` is `x` is primitive."
is_primitive(x::fq_nmod) = multiplicative_order(x) == order(parent(x)) - 1

"""
    divisors(n::Int)

Compute a list of the divisors of `n`.
"""
function divisors(n::Int)

    f = Primes.factor(n)

    nb_factor = 1
    for x in values(f)
        nb_factor *= x + 1
    end

    A = fill(1, nb_factor)
    stop = 1
    for (q, a) in Primes.factor(n)
        for j in 1:a
            for k in 1:stop
                A[j*stop + k] = q^j * A[k]
            end
        end
        stop *= a+1
    end
    return A
end

"""
    make_zetas_prim(p::Int, m::Int = 36)

Compute the minimal polynomial of the root ``ζ_{p^d-1}`` for 
all divisor ``d`` of `m`.

# Note
Take a **primitive element** at random.
"""
function make_zetas_prim(p::Int, m::Int = 36)
    k, x = FiniteField(p, m, "x")
    y = gen(k)
    while !is_primitive(y)
        y = rand(k)
    end

    for d in divisors(m)
        ZETAS[(p, d)] = minpoly(y^(divexact(ZZ(p)^m-1, ZZ(p)^d-1)))
    end
end

"""
    make_zetas_conway(p::Int, m::Int = 100)

Compute the minimal polynomial of the root ``ζ_{p^d-1}`` for 
all 1 <= d <= m for which we have a Conway polynomial.
"""
function make_zetas_conway(p::Int, m::Int = 100)
    for j in 1:m
        b, conw = import_conway(p, j)
        if b
            ZETAS[(p, j)] = conw
        end
    end
end

"""
    nth_root(x::fq_nmod, n::Int)

Compute the `n`-th root of `x`.
"""
function nth_root(x::fq_nmod, n::Int)
    R, T = PolynomialRing(parent(x), "T")
    P = T^n - x
    L = linfactor(P)
    return - coeff(L, 0)
end

"""
    complete_zeta(A::tensor_algebra)

Return the root ``ζ_{p^a-1})``, where ``a`` is the level of `A`.
"""
function complete_zeta(A::tensor_algebra)
    p::Int = characteristic(A.L) 
    l, a = degree(A), level(A)
    K, Z = FiniteField(ZETAS[(p, a)], "Z")
    z = change_basis_inverse(A.R, Z, Z^(divexact(ZZ(p)^a-1, l))) # /!\
    return z
end

"""
    coeff!(x::fq_nmod, j::Int, c::UInt)

Set the `j`-th coeff of `x` to `c`.
"""
function coeff!(x::fq_nmod, j::Int, c::UInt)
    ccall((:nmod_poly_set_coeff_ui, :libflint), Nothing,
          (Ref{fq_nmod}, Int, UInt), x, j, c)
end

"""
    scalar(x::tensor_element)

View `x` as a scalar element.
"""
function scalar(x::tensor_element)
    A = parent(x)
    R = A.R
    s = R()
    for j in 0:level(A)-1
        coeff!(s, j, coeff(coeff(x.elem, j), 0))
    end
    return s
end

"""
    solve_h90(A::tensor_algebra)

Solve H90* in `A` for the root generating the cyclotomic part of `A`.

# Remark
The solution ``α ∈  A`` is such that
``α^l=(1 ⊗ ζ_{p^a-1})^a``, where ``a`` is the level of ``A``.
"""
function solve_h90(A::tensor_algebra)
    g = _solve_h90(A)
    l, a = degree(A), level(A)
    c = scalar(g^l)
    z = complete_zeta(A)^a
    r = nth_root(c^-1 * z, l)
    return r * g
end

"""
    compute_emb(a::fq_nmod, b::fq_nmod)

Compute the embedding ``a\\mapsto b``.
"""
function compute_emb(a::fq_nmod, b::fq_nmod)
    P = minpoly(a)
    L = parent(a)
    K, x = FiniteField(P, "x")
    
    g = change_basis_inverse(K, gen(L), a)
    G = change_basis_direct(g, b)

    f(z::fq_nmod) = change_basis_direct(z, G)
    return f
end

"""
    derive_emb(hl::tensor_element, HL::tensor_element)

Compute the embedding sending the first coordinate of `hl` in the 
base ``(1 ⊗ ζ_{l}^j)_j`` to the first coordinate of `HL` in the 
base ``(1 ⊗ ζ_{m}^(jm/l)_j``.
"""
function derive_emb(hl::tensor_element, HL::tensor_element)
    Al, Am = parent(hl), parent(HL)
    l, m = degree(Al), degree(Am)
    a = coeff(hl, 0)
    b = left(Am)()
    zm = gen(right(Am))
    g = zm^(divexact(m, l))
    for j in 0:level(Am)-1
        if coeff(HL, j) != 0
            tmp = change_basis_inverse(right(Al), zm^j, g)
            b += coeff(HL, j)*coeff(tmp, 0)
        end
    end
    return compute_emb(a, b)
end

"""
    standard_polynomial(p, l)

Compute the `l`-th standard polynomial in characteristic `p`.
"""
function standard_polynomial(p, l)
    k, x = FiniteField(p, l, "x")
    A = tensor_algebra(k)
    h = solve_h90(A)
    return minpoly(coeff(h, 0))
end

"""
Some test function.
"""
function test_fn(HL::tensor_element, Al::tensor_algebra)
    Am = parent(HL)
    l, m = degree(Al), degree(Am)
    b = left(Am)()
    zm = gen(right(Am))
    g = zm^(divexact(m, l))

    for j in 0:level(Am)-1
        if coeff(HL, j) != 0
            tmp = change_basis_inverse(right(Al), zm^j, g)
            b += coeff(HL, j)*coeff(tmp, 0)
        end
    end
    return b
end
