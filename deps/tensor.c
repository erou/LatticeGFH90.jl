#include "tensor.h"

/*
 * Rewrites a, seen as an element of (k[x]/from)[y]/to to an element
 * of (k[y]/to)[x]/from.
 *
 * Requires res to be a zero polynomial.
 */
void _transpose(fq_nmod_poly_t res, const fq_nmod_poly_t a,
	       const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to) {
  fq_nmod_poly_fit_length(res, ctx_from->modulus->length - 1, ctx_to);
  _fq_nmod_poly_set_length(res, ctx_from->modulus->length - 1, ctx_to);
  for (slong i = 0; i < a->length; i++) {
    for (slong j = 0; j < a->coeffs[i].length; j++) {
      nmod_poly_set_coeff_ui(res->coeffs + j, i, a->coeffs[i].coeffs[j]);
    }
  }
  _fq_nmod_poly_normalise(res, ctx_to);
}

/*
 * Computes a*b
 */
void tensor_mul(fq_nmod_poly_t res,
		const fq_nmod_poly_t a, const fq_nmod_poly_t b,
		const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  fq_nmod_poly_t temp;
  
  fq_nmod_poly_mul(res, a, b, L);
    
  fq_nmod_poly_init(temp, R);
  _transpose(temp, res, L, R);
  for (slong i = 0; i < temp->length; i++) {
    fq_nmod_reduce(temp->coeffs + i, R);
  }
  _fq_nmod_poly_set_length(res, 0, L);
  _transpose(res, temp, R, L);

  fq_nmod_poly_clear(temp, R);
}

/*
 * Computes a*b, where a ∈ R and b ∈ L⊗R.
 */
void tensor_scalar_mul_r(fq_nmod_poly_t res,
			 const fq_nmod_t a, const fq_nmod_poly_t b,
			 const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  fq_nmod_poly_t temp;
  
  fq_nmod_poly_init(temp, R);
  _transpose(temp, b, L, R);
  for (slong i = 0; i < temp->length; i++) {
    fq_nmod_mul(temp->coeffs + i, a, temp->coeffs + i, R);
  }
  _fq_nmod_poly_set_length(res, 0, L);
  _transpose(res, temp, R, L);

  fq_nmod_poly_clear(temp, R);
}
