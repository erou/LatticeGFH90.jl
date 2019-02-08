#include <flint/fq_nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/fmpz_poly.h>

#include "h90.h"

int main() {
  fq_nmod_ctx_t L, R;
  flint_rand_t state;
  fmpz_poly_t cyclo;
  nmod_poly_t fact;
  nmod_poly_factor_t f;
  fq_nmod_poly_t h90;

  flint_printf("Testing H90\n");
  flint_randinit(state);
  for (int i = 1 ; i < 40 ; i++) {
    fq_nmod_ctx_randtest(L, state);
    slong n = tensor_degree(L, R);
    if (n == 1) {
      fq_nmod_ctx_clear(L);
      continue;
    }
    
    fmpz_poly_init(cyclo);
    fmpz_poly_cyclotomic(cyclo, n);
    nmod_poly_init(fact, tensor_p(L, R));
    for (int i = fmpz_poly_degree(cyclo); i >= 0; i--) {
      nmod_poly_set_coeff_ui(fact, i, fmpz_poly_get_coeff_ui(cyclo, i));
    }

    nmod_poly_factor_init(f);
    nmod_poly_factor(f, fact);
    nmod_poly_set(fact, f->p);
    fq_nmod_ctx_init_modulus(R, fact, "Z");

    fq_nmod_poly_init(h90, L);
    solve_h90(h90, L, R);
    
    if (!is_h90(h90, L, R)) {
      flint_printf("Not an H90 solution of order %d: ", n);
      fq_nmod_poly_print(h90, L);
      flint_printf("\n");
    } else {
      flint_printf(".%d.", n);
    }

    nmod_poly_clear(fact);
    fmpz_poly_clear(cyclo);
    fq_nmod_poly_clear(h90, L);
    nmod_poly_factor_clear(f);
    fq_nmod_ctx_clear(L);
    fq_nmod_ctx_clear(R);
  }

  flint_randclear(state);
  flint_printf("done\n");
}
