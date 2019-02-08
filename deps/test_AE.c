#include <flint/fq_nmod.h>

#include "AE.h"

int main() {
  fq_nmod_ctx_t ctx;
  fq_nmod_t A, g;
  nmod_poly_t res1, res2;
  flint_rand_t state;

  flint_printf("Testing automorphism evaluation (AE)\n");
  flint_randinit(state);
  for (int i = 1 ; i < 40 ; i++) {
    fq_nmod_ctx_randtest(ctx, state);
    fq_nmod_init(A, ctx);
    fq_nmod_init(g, ctx);
    nmod_poly_init(res1, ctx->mod.n);
    nmod_poly_init(res2, ctx->mod.n);

    fq_nmod_randtest_not_zero(A, state, ctx);
    fq_nmod_randtest_not_zero(g, state, ctx);
      
    compose_naive(res1, A, g, ctx->modulus, ctx->inv);
    compose(res2, A, g, ctx->modulus, ctx->inv);

    if (!nmod_poly_equal(res1, res2)) {
      flint_printf("Not equal:");
      nmod_poly_print(res1);
      flint_printf("\n");
      nmod_poly_print(res2);
      flint_printf("\n");
    } else {
      flint_printf(".");
    }
      
    nmod_poly_clear(res1);
    nmod_poly_clear(res2);      
    fq_nmod_clear(g, ctx);
    fq_nmod_clear(A, ctx);
    fq_nmod_ctx_clear(ctx);
  }

  flint_randclear(state);
  flint_printf("done\n");
}
