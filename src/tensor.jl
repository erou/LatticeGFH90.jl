################################################################
#
# Tensor type
#
################################################################

export tensor_algebra

################################################################
#
# Level
#
################################################################

"Compute the image of `n` by the Euler Totient function ``φ``"
phi_euler(n::Int) = ccall((:n_euler_phi, :libflint), Int, (Int,), n)

"""
level(n::Int, p::Int)

Compute the level associated with the degree `n` in characteristic `p`.

# Note
This is *the order* of `p` in the multiplicative group of ``Z/nZ``.
"""
function level(n::Int, p::Int)
    phi = phi_euler(n)
    order = 1
    for (q, a) in Primes.factor(phi)
        tmp = powmod(p, divexact(phi,q^a), n)
        while tmp != 1
            order *= q
            tmp = powmod(tmp, q, n)
        end
    end
    return order
end

################################################################
#
# Struct definitions
#
################################################################

mutable struct tensor_algebra
    L::FqNmodFiniteField
    R::FqNmodFiniteField

    function tensor_algebra(p::Int, l::Int)
        A = new()
        a = level(l, p)
        !haskey(ZETAS, (p, a)) && throw("NotImplementedError")
        k, x = FiniteField(ZETAS[(p, a)], "x")
        L = FiniteField(p, l, string("x", l))[1]
        R = FiniteField(minpoly(x^(divexact(ZZ(p)^a-1, l))), string("ζ", l))[1]
        A.L = L
        A.R = R
        return A
    end

    function tensor_algebra(L::FqNmodFiniteField)
        A = new()
        p::Int = characteristic(L)
        l::Int = degree(L)
        a = level(l, p)
        !haskey(ZETAS, (p, a)) && throw("NotImplementedError")
        k, x = FiniteField(ZETAS[(p, a)], "x")
        R = FiniteField(minpoly(x^(divexact(ZZ(p)^a-1, l))), string("ζ", l))[1]
        A.L = L
        A.R = R
        return A
    end
end

mutable struct tensor_element
    elem::fq_nmod_poly
    parent::tensor_algebra
end

################################################################
#
# Show
#
################################################################

show(io::IO, x::tensor_element) = show(io, x.elem)

function show(io::IO, A::tensor_algebra)
    show(io, A.L)
    print(io, " ⊗ ")
    show(io, A.R)
end

################################################################
#
# Basic operations
#
################################################################

degree(A::tensor_algebra) = degree(A.L)
level(A::tensor_algebra) = degree(A.R)
parent(x::tensor_element) = x.parent

################################################################
#
# Callable tensor_algebra object
#
################################################################

function (A::tensor_algebra)()
    res = PolynomialRing(A.L, string("ζ", degree(A)))[1]()
    return tensor_element(res, A)
end

(A::tensor_algebra)(x::fq_nmod_poly) = tensor_element(x, A)

function (A::tensor_algebra)(x::Array{fq_nmod, 1})
    res = PolynomialRing(A.L, string("ζ", degree(A)))[1](x)
    return tensor_element(res, A)
end

function (A::tensor_algebra)(x::Array{fmpz, 1})
    res = PolynomialRing(A.L, string("ζ", degree(A)))[1](x)
    return tensor_element(res, A)
end

function (A::tensor_algebra)(x::Array{T, 1}) where {T <: Integer}
    res = PolynomialRing(A.L, string("ζ", degree(A)))[1](x)
    return tensor_element(res, A)
end

################################################################
#
# Binary operations
#
################################################################

function *(x::tensor_element, y::tensor_element)
    A = parent(x)
    z = A()
    ccall((:tensor_mul, :libembed), Nothing, (Ref{fq_nmod_poly}, Ref{fq_nmod_poly},
          Ref{fq_nmod_poly}, Ref{FqNmodFiniteField}, Ref{FqNmodFiniteField}),
          z.elem, x.elem, y.elem, A.L, A.R)
    return z
end

function *(x::fq_nmod, y::tensor_element)
    A = parent(y)
    z = A()
    ccall((:tensor_scalar_mul_r, :libembed), Nothing, (Ref{fq_nmod_poly},
          Ref{fq_nmod}, Ref{fq_nmod_poly}, Ref{FqNmodFiniteField},
          Ref{FqNmodFiniteField}), z.elem, x, y.elem, A.L, A.R)
    return z
end

*(x::tensor_element, y::fq_nmod) = y*x

################################################################
#
# Powering
#
################################################################

function ^(x::tensor_element, n::Int)
    n < 0 && throw(DomainError())
    A = parent(x)
    y = x

    bit = ~((~UInt(0)) >> 1)
    while (UInt(bit) & n) == 0
        bit >>= 1
    end
    bit >>= 1
    while bit != 0
        y = y*y
        if (UInt(bit) & n) != 0
            y *= x
        end
        bit >>= 1
    end
    return y
end

################################################################
#
# Comparison
#
################################################################

==(x::tensor_element, y::tensor_element) = x.elem == y.elem
