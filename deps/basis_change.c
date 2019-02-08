#include <flint/fmpz.h>
#include <flint/nmod_vec.h>
#include "basis_change.h"
#include "minpoly.h"

/**
 * Computes the dual representation of the given polynomial {@code a} modulo {@code ctx}
 */
void monomial_to_dual(mp_limb_t *dual, const nmod_poly_t a, const fq_nmod_ctx_t ctx) {
  nmod_poly_t temp;
  nmod_poly_init(temp, ctx->modulus->mod.n);
  
  slong m = nmod_poly_degree(ctx->modulus);
  
  nmod_poly_derivative(temp, ctx->modulus);
  fq_nmod_mul(temp, temp, a, ctx);
  nmod_poly_reverse(temp, temp, m);
  nmod_poly_mullow(temp, temp, ctx->inv, m);
  
  for (slong i = 0; i < m; i++)
    dual[i] = nmod_poly_get_coeff_ui(temp, i);
  
  nmod_poly_clear(temp);
}

void dual_to_monomial_precomp_deriv_inv(nmod_poly_t result, const mp_limb_t *dual,
					const fq_nmod_ctx_t ctx, const nmod_poly_t deriv_inv) {
  nmod_poly_t temp1;
  nmod_poly_t temp2;
  nmod_poly_init(temp1, ctx->modulus->mod.n);
  nmod_poly_init(temp2, ctx->modulus->mod.n);

  slong m = nmod_poly_degree(ctx->modulus);
  for (slong i = 0; i < m; i++)
    nmod_poly_set_coeff_ui(temp1, i, dual[i]);
  
  nmod_poly_reverse(temp2, ctx->modulus, m + 1);
  nmod_poly_mullow(temp1, temp1, temp2, m);
  
  nmod_poly_reverse(temp2, temp1, m);
  
  fq_nmod_mul(result, deriv_inv, temp2, ctx);
  
  nmod_poly_clear(temp1);
  nmod_poly_clear(temp2);
}

/**
 * Computes the representation of the given dual {@code dual}
 * in monomial basis modulo {@code ctx} 
 */
void dual_to_monomial(nmod_poly_t result, const mp_limb_t *dual, const fq_nmod_ctx_t ctx) {
  nmod_poly_t deriv_inv;

  nmod_poly_init(deriv_inv, ctx->modulus->mod.n);
  // compute 1 / modulus'
  nmod_poly_derivative(deriv_inv, ctx->modulus);
  fq_nmod_inv(deriv_inv, deriv_inv, ctx);

  dual_to_monomial_precomp_deriv_inv(result, dual, ctx, deriv_inv);

  nmod_poly_clear(deriv_inv);
}

/*
 * See doc for change_basis_inverse.
 *
 * Takes two additional precomputed arguments:
 *
 * - to_deriv_inv is 1/P' mod P, where P is the modulus of ctx_to,
 *
 * - an optional trace_one element, that is multiplied by a prior to
 *   conversion.
 */
void change_basis_inverse_precomp(fq_nmod_t res,
				  const fq_nmod_t a, const fq_nmod_t g,
				  const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to,
				  const fq_nmod_t to_deriv_inv, const fq_nmod_t trace_one) {
  slong degree = fq_nmod_ctx_degree(ctx_from);
  
  mp_limb_t *dual = _nmod_vec_init(degree); //new mp_limb_t[degree];
  for (slong i = 0; i < degree; i++)
    dual[i] = 0;

  // compute the dual basis image of a
  if (trace_one) {
    fq_nmod_t temp;
    fq_nmod_init(temp, ctx_from);
    fq_nmod_mul(temp, a, trace_one, ctx_from);
    monomial_to_dual(dual, temp, ctx_from);
    fq_nmod_clear(temp, ctx_from);
  } else {
    monomial_to_dual(dual, a, ctx_from);
  }
  
  // compute the power projection <dual, ctx_to>
  project_powers(dual, dual, degree, g, ctx_from->modulus, ctx_from->inv);

  // compute result such that result(f) = g
  dual_to_monomial_precomp_deriv_inv(res, dual, ctx_to, to_deriv_inv);

  _nmod_vec_clear(dual);
}

/*
 * Given the image g of the canonical generator of ctx_to in ctx_from,
 * compute the inverse modulo ctx_to of the derivative of ctx_to, and
 * an element with relative ctx_from/ctx_to trace 1.
 *
 * These elements are needed as precomputations for
 * change_basis_inverse_precomp, as well as
 * change_basis_inverse_and_project.
 */
void change_basis_prepare(fq_nmod_t deriv_inv, fq_nmod_t trace_one,
			  const fq_nmod_t g,
			  const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to) {
  slong degree = fq_nmod_ctx_degree(ctx_from);
  slong deg_low = fq_nmod_ctx_degree(ctx_to);
  fmpz_t deg_ratio;
  
  // put 1 / modulus' in deriv_inv
  nmod_poly_derivative(deriv_inv, ctx_to->modulus);
  fq_nmod_inv(deriv_inv, deriv_inv, ctx_to);

  // Compute an element of relative trace 1
  fmpz_init_set_ui(deg_ratio, degree / deg_low);
  if (fmpz_invmod(deg_ratio, deg_ratio, &(ctx_to->p))) {
    fq_nmod_set_fmpz(trace_one, deg_ratio, ctx_from);
  } else {
    fq_nmod_t x_trace;
    fq_nmod_init(x_trace, ctx_to);
    fq_nmod_one(trace_one, ctx_from);
    do {
      nmod_poly_shift_left(trace_one, trace_one, 1);
      change_basis_inverse_precomp(x_trace, trace_one, g, ctx_from, ctx_to, deriv_inv, NULL);
    } while(fq_nmod_is_zero(x_trace, ctx_to));
    change_basis_direct(x_trace, x_trace, ctx_to, g, ctx_from);
    fq_nmod_div(trace_one, trace_one, x_trace, ctx_from);
    fq_nmod_clear(x_trace, ctx_from);
  }
}

/*
 * Given a ∈ ctx_from, and given the image g of the canonical
 * generator of ctx_to in ctx_from, compute the image of a in ctx_to.
 *
 * ctx_to can be smaller than ctx_from, but not bigger.
 */
void change_basis_inverse(fq_nmod_t res,
			  const fq_nmod_t a, const fq_nmod_t g,
			  const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to) {
  nmod_poly_t deriv_inv;
  fq_nmod_t trace_one;

  fq_nmod_init(trace_one, ctx_from);
  nmod_poly_init(deriv_inv, ctx_to->modulus->mod.n);

  change_basis_prepare(deriv_inv, trace_one, g, ctx_from, ctx_to);
  change_basis_inverse_precomp(res, a, g, ctx_from, ctx_to, deriv_inv, trace_one);
  
  fq_nmod_clear(trace_one, ctx_from);
  nmod_poly_clear(deriv_inv);
}

/*
 * Given a ∈ ctx_from, and given the image g of the canonical
 * generator of ctx_from in ctx_to, compute the image of a in ctx_to.
 *
 * I.e., compute a(g) mod ctx_to.
 *
 * ctx_from is unused and a NULL pointer can be passed in its place.
 */
void change_basis_direct(fq_nmod_t res,
			 const fq_nmod_t a, const fq_nmod_ctx_t ctx_from_unused,
			 const fq_nmod_t g, const fq_nmod_ctx_t ctx_to) {
  nmod_poly_compose_mod_brent_kung_preinv(res, a, g, ctx_to->modulus, ctx_to->inv);
}


/*
 * Given an element g ∈ ctx with minpoly min, compute the vector of
 * the linear form ⎣.⎦_g ∘ Tr, where Tr is the trace from ctx to k(g),
 * and ⎣.⎦_g is the projection on the first component in the basis g.
 */
void project_tr(mp_limb_t * res, const fq_nmod_t g,
		const fq_nmod_ctx_t minpoly, const nmod_poly_t minpoly_deriv_inv,
		const fq_nmod_ctx_t ctx) {
  fq_nmod_t hprime, gshift;
  slong d = fq_nmod_ctx_degree(ctx);

  fq_nmod_init(hprime, ctx);
  fq_nmod_init(gshift, minpoly);
  
  // rev(ctx' · ((minpoly ÷ x) / minpoly')(g)) / rev(ctx)  mod x^d
  
  nmod_poly_derivative(hprime, ctx->modulus);
  nmod_poly_shift_right(gshift, minpoly->modulus, 1);
  fq_nmod_mul(gshift, gshift, minpoly_deriv_inv, minpoly);
  nmod_poly_compose_mod_brent_kung_preinv(gshift, gshift, g, ctx->modulus, ctx->inv);
  fq_nmod_mul(gshift, hprime, gshift, ctx);
  nmod_poly_reverse(gshift, gshift, d);
  _nmod_poly_mullow(res, gshift->coeffs, gshift->length,
		    ctx->inv->coeffs, ctx->inv->length, d, ctx->mod);

  fq_nmod_clear(hprime, ctx);
  fq_nmod_clear(gshift, minpoly);
}

/*
 * Inputs:
 *  - an array `polys` of lenght `n` of element of `ctx_to`
 *     expressed in `ctx_from`,
 *  - the canonical generator `g` of `ctx_to` expressed in `ctx_from`,
 *  - the values `to_deriv_inv` and `trace_one` precomputed 
 *    by `change_basis_prepare`.
 * 
 * Output the constant coefficients of `polys` expressed in the basis
 * g.
 *
 * Note, the output `res` must be pre-allocated to length n.
 */
void change_basis_inverse_and_project_precomp(mp_limb_t * res,
					      const fq_nmod_struct * polys, slong n,
					      const fq_nmod_t g, const fq_nmod_ctx_t ctx_from,
					      const fq_nmod_ctx_t ctx_to,
					      const fq_nmod_t to_deriv_inv,
					      const fq_nmod_t trace_one) {
  slong d = fq_nmod_ctx_degree(ctx_from);
  mp_limb_t form[d];
  
  project_tr(form, g, ctx_to, to_deriv_inv, ctx_from);
  transposed_mulmod(form, form, trace_one, ctx_from->modulus, ctx_from->inv);

  slong nlimbs = _nmod_vec_dot_bound_limbs(d, ctx_from->mod);
  for (slong i = 0; i < n; i++) {
    res[i] = _nmod_vec_dot(form, polys[i].coeffs, polys[i].length, ctx_from->mod, nlimbs);
  }
}

/*
 * Inputs:
 *  - an array `polys` of lenght `n` of element of `ctx_to`
 *     expressed in `ctx_from`,
 *  - the canonical generator `g` of `ctx_to` expressed in `ctx_from`.
 * 
 * Output the constant coefficients of `polys` expressed in the basis
 * g.
 *
 * Note, the output `res` must be pre-allocated to length n.
 */
void change_basis_inverse_and_project(mp_limb_t * res,
				      const fq_nmod_struct * polys, slong n,
				      const fq_nmod_t g, const fq_nmod_ctx_t ctx_from,
				      const fq_nmod_ctx_t ctx_to) {
  nmod_poly_t deriv_inv;
  fq_nmod_t trace_one;

  fq_nmod_init(trace_one, ctx_from);
  nmod_poly_init(deriv_inv, ctx_to->modulus->mod.n);

  change_basis_prepare(deriv_inv, trace_one, g, ctx_from, ctx_to);
  change_basis_inverse_and_project_precomp(res, polys, n, g, ctx_from, ctx_to, deriv_inv, trace_one);
  
  fq_nmod_clear(trace_one, ctx_from);
  nmod_poly_clear(deriv_inv);
}
