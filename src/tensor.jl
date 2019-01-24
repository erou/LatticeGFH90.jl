################################################################
#
# Tensor type
#
################################################################

import Base: show, *, ^, ==

export tensor_algebra, tensor_element

################################################################
#
# Struct definitions
#
################################################################

mutable struct tensor_algebra
    L::FqNmodFiniteField
    R::FqNmodFiniteField
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
# one line functions
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
    res = PolynomialRing(A.L, "T")[1]()
    return tensor_element(res, A)
end

(A::tensor_algebra)(x::fq_nmod_poly) = tensor_element(x, A)

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

################################################################
#
# Comparison
#
################################################################

==(x::tensor_element, y::tensor_element) = x.elem == y.elem
