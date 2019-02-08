#include "minpoly.h"

int main() {
  flint_rand_t state;
  fq_nmod_ctx_t ctx;
  fq_nmod_t a, eval;
  nmod_poly_t min;
  
  flint_printf("Testing minpoly\n");
  flint_randinit(state);
  for (int i = 1 ; i < 40 ; i++) {
    fq_nmod_ctx_randtest(ctx, state);
    
    fq_nmod_init(a, ctx);
    fq_nmod_randtest(a, state, ctx);
    nmod_poly_init(min, ctx->mod.n);
    minpoly(min, a, ctx);

    fq_nmod_init(eval, ctx);
    fq_nmod_zero(eval, ctx);
    for (int i = min->length - 1; i >= 0; i--) {
      fq_nmod_mul(eval, eval, a, ctx);
      eval->coeffs[0] = nmod_add(eval->coeffs[0], min->coeffs[i], eval->mod);
    }
    
    if (!fq_nmod_is_zero(eval, ctx)) {
      flint_printf("Not the minimal polynomial: (");
      nmod_poly_print_pretty(min, "X");
      flint_printf(")(");
      fq_nmod_print_pretty(a, ctx);
      flint_printf(") = ");
      fq_nmod_print_pretty(eval, ctx);
      flint_printf("\n");      
    } else {
      flint_printf(".");
    }

    nmod_poly_clear(min);
    fq_nmod_clear(eval, ctx);
    fq_nmod_clear(a, ctx);
    fq_nmod_ctx_clear(ctx);
  }

  flint_randclear(state);
  flint_printf("done\n");
}
