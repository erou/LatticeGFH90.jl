#ifndef BASIS_CHANGE_H
#define BASIS_CHANGE_H

#include <flint/fq_nmod.h>

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
			  const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to);

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
				  const fq_nmod_t to_deriv_inv, const fq_nmod_t trace_one);


/*
 * Given a ∈ ctx_from, and given the image g of the canonical
 * generator of ctx_to in ctx_from, compute the image of a in ctx_to.
 *
 * ctx_to can be smaller than ctx_from, but not bigger. In this case,
 * deg(ctx_from)/deg(ctx_to) must not be divisible by the
 * characteristic, othwerise 0 is returned.
 */
void change_basis_inverse(fq_nmod_t res,
			  const fq_nmod_t a, const fq_nmod_t g,
			  const fq_nmod_ctx_t ctx_from, const fq_nmod_ctx_t ctx_to);

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
			 const fq_nmod_t g, const fq_nmod_ctx_t ctx_to);


/*
 * Inputs:
 *  - an array `polys` of length `n` of element of `ctx_to`
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
					      const fq_nmod_t trace_one);

/*
 * Inputs:
 *  - an array `polys` of length `n` of element of `ctx_to`
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
				      const fq_nmod_ctx_t ctx_to);

#endif
