#ifndef H90_H
#define H90_H

#include <flint/fq_nmod_poly.h>

#include "tensor.h"

/*
 * Compute the polynomial Σ a_i X^i of degree deg(A.R) with
 *
 *   a_{k-1} = a
 *   a_i = a_{i+1}^p + A.R[i+1] a
 *
 * Can use to reconstruct a H90 solution from its (k-1)-th coefficient.
 */
void lift_h90(fq_nmod_poly_t res,
	      const fq_nmod_t a, const fq_nmod_ctx_t L, const fq_nmod_ctx_t R);

/*
 * Find a solution res to the hilbert 90 equation
 * 
 *   (σ⊗1)(res) = ζ · res
 *
 * in the algebra L⊗R, where ζ is the canonical root of R
 */
void solve_h90(fq_nmod_poly_t res,
	       const fq_nmod_ctx_t L, const fq_nmod_ctx_t R);
/*
 * Tests if x∈A is a soultion of Hilbert 90 of order ord(A).
 */
int is_h90(const fq_nmod_poly_t x, const fq_nmod_ctx_t L, const fq_nmod_ctx_t R);

#endif
