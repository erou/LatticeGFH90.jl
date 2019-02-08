#include "h90.h"
#include "AE.h"

/*
 * Compute the polynomial Σ a_i X^i of degree deg(R) with
 *
 *   a_{k-1} = a
 *   a_i = a_{i+1}^p + R[i+1] a
 *
 * Can use to reconstruct a H90 solution from its (k-1)-th coefficient.
 */
void lift_h90(fq_nmod_poly_t res,
	      const fq_nmod_t a, const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  slong k = tensor_level(L, R);

  if (fq_nmod_is_zero(a, L)) {
    fq_nmod_poly_zero(res, L);
    return;
  }
  
  fq_nmod_t temp;
  fq_nmod_init(temp, L);
  
  fq_nmod_poly_set_coeff(res, k-1, a, L);
  for (slong i = k-2; i >= 0; i--) {
    fq_nmod_frobenius(temp, res->coeffs + i + 1, 1, L);
    fq_nmod_mul_ui(res->coeffs + i, a, nmod_poly_get_coeff_ui(R->modulus, i+1), L);
    fq_nmod_add(res->coeffs + i, res->coeffs + i, temp, L);
  }

  fq_nmod_clear(temp, L);
}


/*
 * Evaluate
 * 
 *   Σ σ^i(a) ⊗ ζ^(-i)
 *
 * where a ∈ L, σ is the Frobenius of L, and ζ is the canonical root
 * of R.
 *
 * phi_div_R is the polynomial (X^n - 1) / R, where n is the order of
 * ζ and the degree of L.
 */
void _eval_cycle(fq_nmod_poly_t res, const fq_nmod_t a,
		 const fq_nmod_ctx_t L, const fq_nmod_ctx_t R, const nmod_poly_t phi_div_R) {
  nmod_poly_t b;
  
  nmod_poly_init(b, tensor_p(L, R));
  // eval phi_div_R(σ)(a) mod L
  compose(b, phi_div_R, a, L->modulus, L->inv);
  lift_h90(res, b, L, R);

  nmod_poly_clear(b);
}

/*
 * Find a solution res to the hilbert 90 equation
 * 
 *   (σ⊗1)(res) = ζ · res
 *
 * in the algebra L⊗R, where ζ is the canonical root of R
 */
void solve_h90(fq_nmod_poly_t res,
	       const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  fq_nmod_t a;
  nmod_poly_t g;
  flint_rand_t state;
  slong n = tensor_degree(L, R);

  // Set g to (X^n - 1) / R
  // Can't this be improved?
  nmod_poly_init2_preinv(g, tensor_p(L, R), tensor_pinv(L, R), n+1);
  nmod_poly_set_coeff_ui(g, n, 1);
  nmod_poly_set_coeff_ui(g, 0, -1);
  nmod_poly_div(g, g, R->modulus);

  fq_nmod_init(a, L);
  flint_randinit(state);
  do {
    fq_nmod_randtest(a, state, L);
    _eval_cycle(res, a, L, R, g);
  } while(fq_nmod_poly_is_zero(res, L));

  flint_randclear(state);
  fq_nmod_clear(a, L);
  nmod_poly_clear(g);
}

/*
 * Tests if x∈L⊗R is a solution of Hilbert 90 of order deg(L).
 */
int is_h90(const fq_nmod_poly_t x, const fq_nmod_ctx_t L, const fq_nmod_ctx_t R) {
  int ret = 1;
  slong k = tensor_level(L, R);
  fq_nmod_t temp1, temp2;
  fq_nmod_init(temp1, L);
  fq_nmod_init(temp2, L);
  
  for (slong i = 0; i < k; i++) {
    fq_nmod_mul_ui(temp1, x->coeffs + k - 1,
		   nmod_poly_get_coeff_ui(R->modulus, i), L);
    if (i > 0) {
      fq_nmod_sub(temp1, temp1, x->coeffs + i - 1, L);
    }
    fq_nmod_frobenius(temp2, x->coeffs + i, 1, L);
    fq_nmod_add(temp1, temp1, temp2, L);
    if (!fq_nmod_is_zero(temp1, L)) {
      ret = 0;
      break;
    }
  }
  
  fq_nmod_clear(temp2, L);
  fq_nmod_clear(temp1, L);
  return ret;
}
