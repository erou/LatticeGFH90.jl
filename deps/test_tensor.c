#include <time.h>
#include "tensor.h"

int main() {
  fq_nmod_ctx_t L, R;
  flint_rand_t state;
  fq_nmod_poly_t a, b, c, d;
  nmod_poly_t r;
  fq_nmod_poly_t mod, inv;
  fmpz_t coeff;
  int primes[3] = {3, 11, 31};
  time_t diff;
  
  flint_printf("Testing tensor mul\n");
  flint_randinit(state);
  for (int i = 100 ; i < 150 ; i++) {
    fmpz_init(coeff);
    fmpz_set_si(coeff, primes[i % 3]);
    fq_nmod_ctx_init(L, coeff, i, "x");
    nmod_poly_init(r, tensor_p(L, R));
    nmod_poly_randtest_monic(r, state, tensor_degree(L, R) + 1);
    fq_nmod_ctx_init_modulus(R, r, "z");

    fq_nmod_poly_init(a, L);
    fq_nmod_poly_init(b, L);
    fq_nmod_poly_init(c, L);
    fq_nmod_poly_init(d, L);
    
    fq_nmod_poly_randtest(a, state, tensor_level(L, R), L);
    fq_nmod_poly_randtest(b, state, tensor_level(L, R), L);
    diff = clock();
    tensor_mul(c, a, b, L, R);
    diff += -clock();

    fq_nmod_poly_init(mod, L);
    for (slong i = 0; i < R->modulus->length; i++) {
      fmpz_set_si(coeff, R->modulus->coeffs[i]);
      fq_nmod_poly_set_coeff_fmpz(mod, i, coeff, L);
    }
    fq_nmod_poly_init(inv, L);
    for (slong i = 0; i < R->inv->length; i++) {
      fmpz_set_si(coeff, R->inv->coeffs[i]);
      fq_nmod_poly_set_coeff_fmpz(inv, i, coeff, L);
    }
    diff -= clock();
    fq_nmod_poly_mulmod_preinv(d, a, b, mod, inv, L);
    diff += clock();
    
    if (!fq_nmod_poly_equal(c, d, L)) {
      flint_printf("Polynomials differ: ");
      fq_nmod_poly_print_pretty(c, "Z", L);
      flint_printf("   !=   ");
      fq_nmod_poly_print_pretty(d, "Z", L);
      flint_printf("\n");
    } else {
      flint_printf(".");
    }
    if (diff < -1000) {
      flint_printf("Warning: naive mul is faster by %g clocks", (double)diff);
    }

    nmod_poly_clear(r);
    fmpz_clear(coeff);
    fq_nmod_poly_clear(mod, L);
    fq_nmod_poly_clear(inv, L);
    fq_nmod_poly_clear(a, L);
    fq_nmod_poly_clear(b, L);
    fq_nmod_poly_clear(c, L);
    fq_nmod_poly_clear(d, L);
    fq_nmod_ctx_clear(L);
    fq_nmod_ctx_clear(R);
  }

  flint_randclear(state);
  flint_printf("done\n");
}
