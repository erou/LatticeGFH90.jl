#include "linfactor.h"

void linfactor(fq_nmod_poly_t linfactor, const fq_nmod_poly_t input,
               const fq_nmod_ctx_t ctx){
    if (input->length == 2)
    {
        fq_nmod_poly_set(linfactor, input, ctx);
    }
    else
    {
        flint_rand_t state;
        ulong deflation;
        fq_nmod_poly_t pol;

        flint_randinit(state);
        fq_nmod_poly_init(pol, ctx);
        deflation = fq_nmod_poly_deflation(input, ctx);

        if (deflation == 1 || deflation == fq_nmod_poly_degree(input, ctx))
        {
            fq_nmod_poly_set(pol, input, ctx);

            while (fq_nmod_poly_degree(linfactor, ctx) != 1) {
                while (!fq_nmod_poly_factor_equal_deg_prob(linfactor, state, pol, 1, ctx))
                {
                };
                fq_nmod_poly_set(pol, linfactor, ctx);
            }   

        }
        else
        {
            fq_nmod_poly_deflate(pol, input, deflation, ctx);

            while (fq_nmod_poly_degree(pol, ctx) != 1) {
                while (!fq_nmod_poly_factor_equal_deg_prob(linfactor, state, pol, 1, ctx))
                {
                };
                fq_nmod_poly_set(pol, linfactor, ctx);
            }   
            
            fq_nmod_poly_inflate(pol, linfactor, deflation, ctx);

            while (fq_nmod_poly_degree(pol, ctx) != 1) {
                while (!fq_nmod_poly_factor_equal_deg_prob(linfactor, state, pol, 1, ctx))
                {
                };
                fq_nmod_poly_set(pol, linfactor, ctx);
            }   

        }

        flint_randclear(state);
        fq_nmod_poly_clear(pol, ctx);

    }
}
