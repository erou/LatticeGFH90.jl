#ifndef TENSOR_H
#define TENSOR_H

#include <flint/fq_nmod_poly.h>

static inline slong tensor_degree(const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  return fq_nmod_ctx_degree(L);
}
static inline slong tensor_level(const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  return fq_nmod_ctx_degree(R);
}
static inline mp_limb_t tensor_p(const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  return L->mod.n;
}
static inline mp_limb_t tensor_pinv(const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  return L->mod.ninv;
}

/*
 * Rewrites a, seen as an element of (k[x]/from)[y]/to to an element
 * of (k[y]/to)[x]/from.
 *
 * Requires res to be a zero polynomial.
 */
void _transpose(fq_nmod_poly_t res, const fq_nmod_poly_t a,
		const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to);

/*
 * Computes a*b
 */
void tensor_mul(fq_nmod_poly_t res,
		const fq_nmod_poly_t a, const fq_nmod_poly_t b,
		const fq_nmod_ctx_t L, const fq_nmod_ctx_t R);

/*
 * Computes a*b, where a ∈ R and b ∈ L⊗R.
 */
void tensor_scalar_mul_r(fq_nmod_poly_t res,
			 const fq_nmod_t a, const fq_nmod_poly_t b,
			 const fq_nmod_ctx_t L, const fq_nmod_ctx_t R);

#endif
