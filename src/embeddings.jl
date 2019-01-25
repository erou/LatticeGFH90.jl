########################################################
#
# Embeddings
#
########################################################

"""
    embed(K::FqNmodFiniteField, L::FqNmodFiniteField)

Compute an embedding from `K` to `L`.
"""
function embed(K::FqNmodFiniteField, L::FqNmodFiniteField)
    l, m = degree(K), degree(L)
    p::Int = characteristic(K)

    if haskey(EMBEDDINGS, (p, l, m))
        return EMBEDDINGS[(p, l, m)]
    end

    if haskey(H90_ELEMENTS, (p, m))
        hm = H90_ELEMENTS[(p, m)]
    else
        Am = tensor_algebra(L)
        hm = solve_h90(Am)
        H90_ELEMENTS[(p, m)] = hm
    end

    if haskey(H90_ELEMENTS, (p, l))
        hl = H90_ELEMENTS[(p, l)]
    else
        Al = tensor_algebra(K)
        hl = solve_h90(Al)
        H90_ELEMENTS[(p, l)] = hl
    end

    a, b = level(Al), level(Am)

    pow = divexact((b-a)*ZZ(p)^(a+b) - b*ZZ(p)^b + a*ZZ(p)^a, (ZZ(p)^a-1)*l)
    z = (complete_zeta(Am)^pow)^-1
    HL = z * hm^(divexact(m, l))
    f = derive_emb(hl, HL)
    EMBEDDINGS[(p, l, m)] = f
    return f
end
