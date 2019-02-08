#ifndef AE_H
#define AE_H

#include <flint/nmod_poly.h>

/*
 * Compute A(σ)(g) mod f, naively.
 *
 * A, g, f are polynomials, σ is the Frobenius morphism of the base
 * field, f_inv is the precomputed inverse of the revers of f.
 */
void compose_naive(nmod_poly_t res,
		   const nmod_poly_t A, const nmod_poly_t g,
		   const nmod_poly_t f, const nmod_poly_t f_inv);

/*
 * Compute A(σ)(g) mod f, using the Kaltofen-Shoup algorithm.
 *
 * A, g, f are polynomials, σ is the Frobenius morphism of the base
 * field, f_inv is the precomputed inverse of the revers of f.
 */
void compose(nmod_poly_t res,
	     const nmod_poly_t A, const nmod_poly_t g,
	     const nmod_poly_t f, const nmod_poly_t f_inv);

#endif
