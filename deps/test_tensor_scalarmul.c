#include "tensor.h"

int main() {
  fq_nmod_ctx_t L, R;
  flint_rand_t state;
  fq_nmod_t a;
  fq_nmod_poly_t b, c, d, one;
  nmod_poly_t r;
  fmpz_t coeff;
  int primes[3] = {3, 11, 31};
  
  flint_printf("Testing tensor scalar mul\n");
  flint_randinit(state);
  for (int i = 2 ; i < 50 ; i++) {
    fmpz_init(coeff);
    fmpz_set_si(coeff, primes[i % 3]);
    fq_nmod_ctx_init(L, coeff, i, "x");
    nmod_poly_init(r, tensor_p(L, R));
    nmod_poly_randtest_monic(r, state, tensor_degree(L, R) + 1);
    fq_nmod_ctx_init_modulus(R, r, "z");

    fq_nmod_init(a, R);
    fq_nmod_poly_init(b, L);
    fq_nmod_poly_init(c, L);
    fq_nmod_poly_init(d, L);
    fq_nmod_poly_init(one, L);
    
    fq_nmod_one(a, R);
    fq_nmod_poly_set_coeff(one, 0, a, L);

    fq_nmod_randtest(a, state, R);
    fq_nmod_poly_randtest(b, state, tensor_level(L, R), L);

    tensor_scalar_mul_r(c, a, b, L, R);

    tensor_scalar_mul_r(d, a, one, L, R);
    tensor_mul(d, d, b, L, R);
    
    if (!fq_nmod_poly_equal(c, d, L)) {
      flint_printf("Polynomials differ: ");
      fq_nmod_poly_print_pretty(c, "Z", L);
      flint_printf("   !=   ");
      fq_nmod_poly_print_pretty(d, "Z", L);
      flint_printf("\n");
    } else {
      flint_printf(".");
    }

    fmpz_clear(coeff);
    nmod_poly_clear(r);
    fq_nmod_clear(a, R);
    fq_nmod_poly_clear(b, L);
    fq_nmod_poly_clear(c, L);
    fq_nmod_poly_clear(d, L);
    fq_nmod_poly_clear(one, L);
    fq_nmod_ctx_clear(L);
    fq_nmod_ctx_clear(R);
  }

  flint_randclear(state);
  flint_printf("done\n");
}
