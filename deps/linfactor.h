#ifndef LINFACTOR_H
#define LINFACTOR_H

#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>

void linfactor(fq_nmod_poly_t linfactor, const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx); 

#endif
