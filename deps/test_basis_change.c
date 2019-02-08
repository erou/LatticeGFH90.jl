#include "basis_change.h"
#include "minpoly.h"

int main() {
  flint_rand_t state;
  fq_nmod_ctx_t ctx, ctx_g;
  fq_nmod_t a, g, aa, aaa, temp;
  nmod_poly_t modulus;
  fmpz_t p;
  int primes[5] = {2, 7, 11, 31, 59};
  int degrees[4] = {2, 7, 11, 15};
  slong deg, deg1, deg2;
  
  flint_printf("Testing basis change\n");
  flint_randinit(state);
  for (int i = 1 ; i < 40 ; i++) {
    fmpz_init(p);
    fmpz_set_si(p, primes[i % 5]);
    deg1 = degrees[i % 4];
    deg2 = degrees[i % 3];
    deg = deg1*deg2;
    fq_nmod_ctx_init(ctx, p, deg, "x");

    deg = fq_nmod_ctx_degree(ctx);
    
    fq_nmod_init(a, ctx);
    fq_nmod_init(g, ctx);
    fq_nmod_randtest_not_zero(a, state, ctx);

    nmod_poly_init(modulus, ctx->mod.n);
    do {
      fq_nmod_randtest_not_zero(g, state, ctx);
      minpoly(modulus, g, ctx);
    } while (nmod_poly_degree(modulus) != deg);
    fq_nmod_ctx_init_modulus(ctx_g, modulus, "z");

    fq_nmod_init(aa, ctx_g);
    fq_nmod_init(aaa, ctx);
    change_basis_inverse(aa, a, g, ctx, ctx_g);
    change_basis_direct(aaa, aa, ctx_g, g, ctx);
    
    if (!fq_nmod_equal(a, aaa, ctx)) {
      flint_printf("Elements differ (%d): ", deg);
      fq_nmod_print_pretty(a, ctx);
      flint_printf("   !=   ");
      fq_nmod_print_pretty(aaa, ctx);
      flint_printf("\n");      
    } else {
      flint_printf(".%d.", deg);
    }

    fq_nmod_init(temp, ctx);
    do {
      fq_nmod_randtest_not_zero(g, state, ctx);
      fq_nmod_set(temp, g, ctx);
      for (int j = 1; j < deg1; j++) {
	fq_nmod_frobenius(temp, temp, deg2, ctx);
	fq_nmod_mul(g, g, temp, ctx);
      }
      minpoly(modulus, g, ctx);
    } while (nmod_poly_degree(modulus) != deg2);
    fq_nmod_ctx_init_modulus(ctx_g, modulus, "z");
    
    fq_nmod_zero(a, ctx);
    for (int j = 0; j < deg2; j++) {
      fq_nmod_mul(a, a, g, ctx);
      nmod_poly_set_coeff_ui(a, 0, n_randlimb(state));
    }

    change_basis_inverse(aa, a, g, ctx, ctx_g);
    change_basis_direct(aaa, aa, ctx_g, g, ctx);
    
    if (!fq_nmod_equal(a, aaa, ctx)) {
      flint_printf("Elements differ (%d,%d): ", deg, deg2);
      fq_nmod_print_pretty(a, ctx);
      flint_printf("   !=   ");
      fq_nmod_print_pretty(aaa, ctx);
      flint_printf("\n");      
    } else {
      flint_printf(".%d,%d,%d.", deg, deg2, deg1 % primes[i%5]);
    }

    mp_limb_t coeff;
    change_basis_inverse_and_project(&coeff, a, 1, g, ctx, ctx_g);

    if (coeff != nmod_poly_get_coeff_ui(aa, 0)) {
      flint_printf("Constant coefficients differ %d != %d\n",
		   coeff, nmod_poly_get_coeff_ui(aa, 0));
    } else {
      flint_printf("!");
    }
    
    nmod_poly_clear(modulus);
    fq_nmod_clear(aaa, ctx);
    fq_nmod_clear(aa, ctx_g);
    fq_nmod_clear(a, ctx);
    fq_nmod_clear(g, ctx);
    fq_nmod_ctx_clear(ctx_g);
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(p);
  }

  flint_randclear(state);
  flint_printf("done\n");
}
