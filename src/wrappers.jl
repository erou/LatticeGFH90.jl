##################################
#
# Wrappers for C code
#
##################################

export solve_h90

function solve_h90(A::tensor_algebra)
    res = A()
    ccall((:solve_h90, :libembed), Nothing, (Ref{fq_nmod_poly},
          Ref{FqNmodFiniteField}, Ref{FqNmodFiniteField}), res.elem, A.L, A.R)
    return res
end

function lift_h90(x::fq_nmod, A::tensor_algebra)
    res = A()
    ccall((:lift_h90, :libembed), Nothing, (Ref{fq_nmod_poly},
          Ref{fq_nmod}, Ref{FqNmodFiniteField}, Ref{FqNmodFiniteField}),
          res.elem, x, A.L, A.R)
    return res
end

function minpoly(x::fq_nmod)
    result = PolynomialRing(ResidueRing(ZZ, characteristic(parent(x))), "T")()
    ccall((:minpoly, :libembed), Nothing, (Ref{fq_nmod_poly}, Ref{fq_nmod}, 
          Ref{FqNmodFiniteField}), result, x, parent(x))
    return result
end

function phi_euler(n::Int) = ccall((:n_euler_phi, :libflint), Int, (Int,), n)
