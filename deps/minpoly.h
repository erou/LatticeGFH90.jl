#ifndef MINPOLY_H
#define MINPOLY_H

#include <flint/fq_nmod.h>
#include <flint/nmod_poly_mat.h>

/**
 * Computes the power projection:
 * \[\langle v, 1 \rangle, \langle v, h \rangle, \dots, \langle v, h^{l - 1} \rangle\]
 * where $\langle ., . \rangle$ is the inner product, considering $h^i$ as the vector
 * of its coefficients.
 * 
 * @param modulus_inv_rev	1 / rev(m + 1, modulus) mod x^{m - 1} where m = deg(modulus)
 */
void project_powers(mp_limb_t *result, const mp_limb_t *a, slong l, const nmod_poly_t h,
		    const nmod_poly_t modulus, const nmod_poly_t modulus_inv_rev);

/**
 * Computes the transposed modular product of {@code a} and {@code b} modulo {@code mod}: b°a.
 * 
 * @param a			a vector of length deg(mod)
 * @param b			a polynomial of degree < deg(mod)
 * @param mod			the modulus
 * @param mod_inv_rev		1/rev(mod, n) mod x^{n - 1} where n = deg(mod)
 */
void transposed_mulmod(mp_limb_t *result, const mp_limb_t *a, const nmod_poly_t b, const nmod_poly_t mod,
		       const nmod_poly_t mod_rev_inv);

/**
 * Computes the minimal polynomial of {@code f} modulo {@code ctx}.
 * 
 * @param result	the minimal polynomial of {@code f}
 */
void minpoly(nmod_poly_t result,
	     const fq_nmod_t f, const fq_nmod_ctx_t ctx);

#endif
