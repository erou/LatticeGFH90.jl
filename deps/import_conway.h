#ifndef IMPORT_CONWAY_H
#define IMPORT_CONWAY_H

#include <flint/nmod_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <string.h>

int import_conway(nmod_poly_t poly, const fmpz_t, slong d);

#endif
