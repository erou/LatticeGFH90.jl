#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>

#include "AE.h"

typedef struct {
  long n, k, m, len_f_inv;
  mp_ptr f, f_inv, giant; 
  nmod_mat_t A;  // baby steps
  nmod_t mod;
} brent_kung_struct;

/*------------------------------------------------------*/
/* helper function for modular composition              */
/* computes res[i] = polys[i](arg) mod f, 0 <= i < len2 */
/*------------------------------------------------------*/
void nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(nmod_poly_struct * res,
							 const brent_kung_struct * data,
							 const nmod_poly_struct * polys,
							 slong len2) {
  nmod_mat_t B, C;
  mp_ptr t, h;
  slong len1;

  h = _nmod_vec_init(2*data->n);
  t = _nmod_vec_init(2*data->n);

  for (long i = 0; i < len2; i++){
    nmod_poly_init2_preinv(res + i, data->mod.n, data->mod.ninv, data->n);
    _nmod_poly_set_length(res + i, data->n);
  }

  nmod_mat_init(B, data->k * len2, data->m, data->mod.n);
  nmod_mat_init(C, data->k * len2, data->n, data->mod.n);

  /* Set rows of B to the segments of polys */
  for (long j = 0; j < len2; j++) {
    len1 = (polys + j)->length;
    for (long i = 0; i < len1 / data->m; i++)
      _nmod_vec_set(B->rows[i + j * data->k], (polys + j)->coeffs + i * data->m, data->m);
    long i = len1 / data->m;
    _nmod_vec_set(B->rows[i + j * data->k], (polys + j)->coeffs + i * data->m, len1 % data->m);
  }

  nmod_mat_mul(C, B, data->A);
      
  /* Evaluate block composition using the Horner scheme */
  for (long j = 0; j < len2; j++) {
    _nmod_vec_set(t, C->rows[j * data->k], data->n);
    for (long i = data->n; i < 2*data->n; i++)
      t[i] = 0;
    for (long i = 1; i < data->k; i++){
      _nmod_poly_mul(h, C->rows[j * data->k + i], data->n, data->giant + (i-1)*data->n,
		     data->n, data->mod);
      _nmod_poly_add(t, h, 2*data->n-1, t, 2*data->n-1, data->mod);
    }
    _nmod_poly_divrem_newton_n_preinv(h, (res + j)->coeffs, t, 2*data->n-1, data->f,
				      data->n+1, data->f_inv, data->len_f_inv, data->mod);
  }
      


  for (long i = 0; i < len2; i++)
    _nmod_poly_normalise(res + i);
    
  _nmod_vec_clear(h);
  _nmod_vec_clear(t);

  nmod_mat_clear(B);
  nmod_mat_clear(C);
}


/*------------------------------------------------------*/
/* setup for modular composition by arg mod poly        */
/* precomputes baby steps and giant steps               */      
/*------------------------------------------------------*/
void _nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(brent_kung_struct * data,
							  mp_srcptr arg, slong len_arg,
							  mp_srcptr poly, slong len_poly,
							  mp_srcptr polyinv, slong len_poly_inv,
							  nmod_t mod_in,
							  long sz){
  data->n = len_poly - 1;
  data->f = _nmod_vec_init(len_poly);
  flint_mpn_copyi(data->f, poly, len_poly);

  data->len_f_inv = len_poly_inv;
  data->f_inv = _nmod_vec_init(data->len_f_inv);
  flint_mpn_copyi(data->f_inv, polyinv, data->len_f_inv);

  data->m = sz;
  data->k = data->n / data->m + 1;
  data->mod = mod_in;
  nmod_mat_init(data->A, data->m, data->n, data->mod.n);

  /* baby steps: set rows of A to powers of arg */
  data->A->rows[0][0] = UWORD(1);
  _nmod_vec_set(data->A->rows[1], arg, len_arg);
  flint_mpn_zero(data->A->rows[1] + len_arg, data->n - len_arg);
  for (long i = 2; i < data->m; i++)
    _nmod_poly_mulmod_preinv(data->A->rows[i], data->A->rows[i - 1], data->n,
			     data->A->rows[1], data->n, data->f, data->n+1,
			     data->f_inv, data->len_f_inv, data->mod);

  /* giant steps */
  data->giant = _nmod_vec_init(data->k*data->n);
  _nmod_poly_mulmod_preinv(data->giant, data->A->rows[data->m - 1], data->n,
			   data->A->rows[1], data->n, data->f, data->n+1,
			   data->f_inv, data->len_f_inv, data->mod);
  for (long i = 1; i < data->k-1; i++)
    _nmod_poly_mulmod_preinv(data->giant + i*data->n, data->giant + (i-1)*data->n,
			     data->n, data->giant, data->n, data->f, data->n+1,
			     data->f_inv, data->len_f_inv, data->mod);
}

/*------------------------------------------------------*/
/* setup for modular composition by arg mod poly        */
/* precomputes baby steps and giant steps               */      
/*------------------------------------------------------*/
void nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(brent_kung_struct * data,
							 const nmod_poly_t arg,
							 const nmod_poly_t poly,
							 const nmod_poly_t polyinv,
							 long sz){

  _nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(data,
						       arg->coeffs, arg->length,
						       poly->coeffs, poly->length,
						       polyinv->coeffs, polyinv->length,
						       arg->mod, sz);
}


void _automorphism_evaluation_compose(mp_ptr res, 
				      mp_srcptr a, slong len_a,
				      mp_srcptr g, 
				      mp_srcptr f, slong len_f,
				      mp_srcptr f_inv, slong len_f_inv, nmod_t mod){


  nmod_mat_t A, B, C;
  mp_ptr t, h;
  slong i, n, m;

  n = len_f - 1;
  m = 0.5*n_sqrt(len_a) + 1;
  //    m = len_a;
  //    m = len_a/2;
  long k = (len_a-1) / m + 1;

  nmod_mat_init(A, m, n, mod.n);
  nmod_mat_init(B, k, m, mod.n);
  nmod_mat_init(C, k, n, mod.n);

  h = _nmod_vec_init(n);
  t = _nmod_vec_init(n);

  /* Set rows of B to the segments of a */
  for (i = 0; i < k-1; i++)
    _nmod_vec_set(B->rows[i], a + i*m, m);
  _nmod_vec_set(B->rows[i], a + i*m, (len_a - (k-1)*m));

  /* Set rows of A to frobeniuses of g */
  flint_mpn_copyi(A->rows[0], g, n);
  for (i = 1; i < m; i++)
    _nmod_poly_powmod_ui_binexp_preinv (A->rows[i], A->rows[i-1], mod.n, 
					f, len_f, f_inv, len_f_inv, mod);

  nmod_mat_mul(C, B, A);

  // get x^{p^m} mod f
  _nmod_poly_powmod_x_ui_preinv(h, mod.n, f, len_f, f_inv, len_f_inv, mod);
  for (i = 1; i < m; i++){
    _nmod_poly_powmod_ui_binexp_preinv (t, h, mod.n, f, len_f, f_inv, len_f_inv, mod);
    _nmod_vec_set(h, t, n);
  }

  /* Evaluate block composition using the Horner scheme */
    
  brent_kung_struct compose;
  _nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(&compose, h, n, f, len_f, f_inv, len_f_inv, mod, 2*n_sqrt(n)+1);

  nmod_poly_t input;
  nmod_poly_init2_preinv(input, mod.n, mod.ninv, n);
  nmod_poly_t output;
  _nmod_vec_set(input->coeffs, C->rows[k - 1], n);
  nmod_poly_init2_preinv(output, mod.n, mod.ninv, n);
  input->length = n;
  _nmod_poly_normalise(input);
    
  for (i = k- 2; i >= 0; i--){
    nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(output, &compose, input, 1);
    _nmod_poly_add(input->coeffs, output->coeffs, n, C->rows[i], n, mod);
    input->length = n;
    _nmod_poly_normalise(input);
  }

  _nmod_vec_set(res, input->coeffs, n);
  nmod_poly_clear(output);
  nmod_poly_clear(input);

  // _nmod_vec_set(res, C->rows[k - 1], n);
  // for (i = k- 2; i >= 0; i--){
  //   _nmod_poly_compose_mod_brent_kung_preinv(t, res, n, h, f, len_f, f_inv, len_f_inv, mod);
  // }


  _nmod_vec_clear(h);
  _nmod_vec_clear(t);

  nmod_mat_clear(A);
  nmod_mat_clear(B);
  nmod_mat_clear(C);
}


/*
 * Compute A(σ)(g) mod f, using the Kaltofen-Shoup algorithm.
 *
 * A, g, f are polynomials, σ is the Frobenius morphism of the base
 * field, f_inv is the precomputed inverse of the revers of f.
 */
void compose(nmod_poly_t res,
	     const nmod_poly_t A, const nmod_poly_t g, const nmod_poly_t f, const nmod_poly_t f_inv){

  slong len_A = A->length;
  slong len_g = g->length;
  slong len_f = f->length;
  slong len = len_f - 1;

  mp_ptr ptr2;

  if (len_f == 0) {
    flint_printf("Exception (Nmod_poly_automorphism_evaluation::compose). Division by zero.\n");
    abort();
  }

  if (len_A == 0 || len_f == 1) {
    nmod_poly_zero(res);
    return;
  }

  if (len_A == 1){
    nmod_poly_scalar_mul_nmod(res, g, A->coeffs[0]);
    return;
  }

  if (res == A || res == g || res == f || res == f_inv){
    flint_printf("Exception (Nmod_poly_automorphism_evaluation::compose). Aliasing not supported.\n");
    abort();
  }

  ptr2 = _nmod_vec_init(len);

  if (len_g <= len){
    flint_mpn_copyi(ptr2, g->coeffs, len_g);
    flint_mpn_zero(ptr2 + len_g, len - len_g);
  }
  else {
    _nmod_poly_rem(ptr2, g->coeffs, len_g, f->coeffs, len_f, res->mod);
  }

  nmod_poly_fit_length(res, len);
  _automorphism_evaluation_compose(res->coeffs,
				   A->coeffs, len_A, 
				   ptr2, 
				   f->coeffs, len_f,
				   f_inv->coeffs, f_inv->length, res->mod);
  res->length = len;
  _nmod_poly_normalise(res);
  _nmod_vec_clear(ptr2);
}


/*
 * Compute A(σ)(g) mod f, naively.
 *
 * A, g, f are polynomials, σ is the Frobenius morphism of the base
 * field, f_inv is the precomputed inverse of the revers of f.
 */
void compose_naive(nmod_poly_t res,
		   const nmod_poly_t A, const nmod_poly_t g, const nmod_poly_t f, const nmod_poly_t f_inv){

  nmod_poly_t g_loc, tmp;
  nmod_poly_init(g_loc, f->mod.n);
  nmod_poly_set(g_loc, g);
  nmod_poly_init(tmp, f->mod.n);

  {
    nmod_poly_scalar_mul_nmod(tmp, g_loc, nmod_poly_get_coeff_ui(A, 0));
    nmod_poly_set(res, tmp);
  }
  for (slong i = 1; i <= nmod_poly_degree(A); i++) {
    nmod_poly_powmod_ui_binexp_preinv(g_loc, g_loc, f->mod.n, f, f_inv);
    nmod_poly_scalar_mul_nmod(tmp, g_loc, nmod_poly_get_coeff_ui(A, i));
    nmod_poly_add(res, res, tmp);
  }

  nmod_poly_clear(tmp);
  nmod_poly_clear(g_loc);
}
