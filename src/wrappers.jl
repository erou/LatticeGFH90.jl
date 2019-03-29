##########################################################
#
# Wrappers for C code
#
##########################################################

function _solve_h90(A::tensor_algebra)
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
    result = PolynomialRing(ResidueRing(ZZ, Int(characteristic(parent(x)))), "T")[1]()
    ccall((:minpoly, :libembed), Nothing, (Ref{nmod_poly}, Ref{fq_nmod}, 
          Ref{FqNmodFiniteField}), result, x, parent(x))
    return result
end

function nth_root2(x::fq_nmod, n::Int)
    F = parent(x)
    y = F()
    ccall((:nth_root, :libembed), Nothing, (Ref{fq_nmod}, Ref{fq_nmod}, UInt,
          Ref{FqNmodFiniteField}), y, x, UInt(n), F)
    return y
end

function linfactor(x::fq_nmod_poly)
    R = parent(x)
    y = R()
    ccall((:linfactor, :libembed), Nothing, (Ref{fq_nmod_poly},
          Ref{fq_nmod_poly}, Ref{FqNmodFiniteField}), y, x, base_ring(R))
    return y
end

function change_basis_inverse(toK::FqNmodFiniteField, a::fq_nmod, g::fq_nmod)
    fromL = parent(a)
    res = toK()
    ccall((:change_basis_inverse, :libembed), Nothing, (Ref{fq_nmod},
          Ref{fq_nmod}, Ref{fq_nmod}, Ref{FqNmodFiniteField},
          Ref{FqNmodFiniteField}), res, a, g, fromL, toK)
    return res
end

function change_basis_direct(a::fq_nmod, g::fq_nmod)
    fromK = parent(a)
    toL = parent(g)
    res = toL()
    ccall((:change_basis_direct, :libembed), Nothing, (Ref{fq_nmod},
          Ref{fq_nmod}, Ref{FqNmodFiniteField}, Ref{fq_nmod},
          Ref{FqNmodFiniteField}), res, a, fromK, g, toL)
    return res
end

function import_conway(p::Int, d::Int)
    poly = PolynomialRing(ResidueRing(ZZ, p), "T")[1]()
    b = ccall((:import_conway, :libembed), Bool, (Ref{nmod_poly}, Ref{fmpz},
              UInt), poly, ZZ(p), d)
    return b, poly
end

function change_basis_prepare(toK::FqNmodFiniteField, g::fq_nmod)
    fromL = parent(g)
    deriv_inv = toK()
    trace_one = fromL()
    ccall((:change_basis_prepare, :libembed), Nothing, (Ref{fq_nmod},
          Ref{fq_nmod}, Ref{fq_nmod}, Ref{FqNmodFiniteField},
          Ref{FqNmodFiniteField}), deriv_inv, trace_one, g, fromL, toK)
    return deriv_inv, trace_one
end

function change_basis_inverse_precomp(a::fq_nmod, g::fq_nmod,
                                      deriv_inv::fq_nmod, trace_one::fq_nmod)
    fromL = parent(a)
    toK = parent(deriv_inv)
    res = toK()
    ccall((:change_basis_inverse_precomp, :libembed), Nothing, (Ref{fq_nmod},
          Ref{fq_nmod}, Ref{fq_nmod}, Ref{FqNmodFiniteField},
          Ref{FqNmodFiniteField}, Ref{fq_nmod}, Ref{fq_nmod}), res, a, g, fromL,
          toK, deriv_inv, trace_one)
    return res
end

function _transpose(x::fq_nmod_poly, to::FqNmodFiniteField)
    from = base_ring(parent(x))
    y = PolynomialRing(to, "T")[1]()
    ccall((:_transpose, :libembed), Nothing, (Ref{fq_nmod_poly},
          Ref{fq_nmod_poly}, Ref{FqNmodFiniteField}, Ref{FqNmodFiniteField}), y,
         x, from, to)
    return y
end

_transpose(x::tensor_element) = _transpose(x.elem, parent(x).R)

function change_basis_inverse_and_project_precomp(h::tensor_element, b::Int, g::fq_nmod,
                                               deriv_inv::fq_nmod, trace_one::fq_nmod)
    toK = parent(deriv_inv)
    res = toK()
    ccall((:change_basis_inverse_and_project_precomp_jl, :libembed), Nothing,
          (Ptr{Nothing}, Ptr{Nothing}, Int, Ref{fq_nmod},
           Ref{FqNmodFiniteField}, Ref{FqNmodFiniteField}, Ref{fq_nmod},
           Ref{fq_nmod}), res.coeffs, h.elem.coeffs, b, g, parent(g), toK, deriv_inv,
           trace_one)
    res.length = degree(toK)
    return res
end

function change_basis_inverse_and_project(HL::tensor_element, Al::tensor_algebra)

    Am = parent(HL)
    l, m = degree(Al), degree(Am)
    res = left(Am)()
    g = gen(right(Am))^(divexact(m, l))
    ccall((:change_basis_inverse_and_project_jl, :libembed), Nothing,
          (Ref{fq_nmod}, Ref{fq_nmod_poly}, Int, Ref{fq_nmod},
           Ref{FqNmodFiniteField}, Ref{FqNmodFiniteField}), res,
          _transpose(HL), m, g, right(Am), right(Al))
    return res
end
