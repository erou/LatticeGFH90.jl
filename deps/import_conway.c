#include "import_conway.h" 

int flint_conway_polynomials [] = {
#include "CPimport.h"
  0
};

int import_conway(nmod_poly_t poly, const fmpz_t p, slong d)
{
    unsigned int position;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        return 0;
    }

    for (position = 0; flint_conway_polynomials[position] != 0; position += 3+flint_conway_polynomials[position+1])
    {
        /* Different prime? */
        if (fmpz_cmp_ui(p, flint_conway_polynomials[position]))
            continue;

        /* Same degree? */
        if (d == flint_conway_polynomials[position+1])
        {
            nmod_poly_t mod;
            slong i;

            nmod_poly_init(mod, fmpz_get_ui(p));
            
            /* Copy the polynomial */
            for (i = 0; i < d; i++)
            {
                int coeff = flint_conway_polynomials[position+2+i];                
                nmod_poly_set_coeff_ui(mod, i, coeff);
            }

            nmod_poly_set_coeff_ui(mod, d, 1);

            nmod_poly_set(poly, mod);
            nmod_poly_clear(mod);
            return 1;
        }
    }

    return 0;
}
