# Benchmarks

## Why are there new files and old files?

There are two types of benchmarks, the only difference is that the ones
containing `new` were made using the more optimized function 
`change_basis_inverse_and_project`. 

## What is benchmarked?

We benchmark two routines:
  - in the files `solve_h90`, we measure the time needed to compute the Kummer
    algebra A_l and to find a solution of
    Hilbert 90 in A_l.
    - the first column is the level of the algebra
    - the second column is the degree of the algebra
    - the third column is the time needed 
  - in the files `embed`, we measure the time needed to compute an embedding
    from GF(p^l) in GF(p^m), when the algebras A_l, A_m and the solutions of
    Hilbert 90 have already been computed
    - the first column is the level of the algebra A_l
    - the second column is the degree of the algebra A_l
    - the third column is the level of the algebra A_m
    - the fourth column is the degree of the algebra A_m
    - the fifth column is the time needed 
