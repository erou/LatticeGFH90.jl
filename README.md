# LatticeGFH90.jl

A [Julia](https://julialang.org/)/[Nemo](http://nemocas.org/) package to make lattices of compatibly embedded finite fields
based on Hilbert's theorem 90 and Allombert's algorithm.

## How to: install

### The Julia packages

- Install Julia (preferably the last version, to avoid suprises...)
- Type `]` to go in package mode (`Backspace` to quit)

```
(v1.1) pkg> add Nemo

(v1.1) pkg> add "https://github.com/erou/LatticeGFH90.jl"

julia> using LatticeGFH90
```

- You're *almost* done: the code relies on a library called `libembed`

### The C library `libembed`

- The sources are in [(yet) another repo](https://github.com/erou/gf-h90-lattice/tree/master/implementation)
- Once you have the sources:
```
$ gcc -fPIC -c AE.c h90.c minpoly.c tensor.c nth-root.c linfactor.c basis_change.c
$ gcc -shared AE.o h90.o minpoly.o tensor.o nth-root.o basis_change.o -o libembed.so
```
- Copy the headers and the shared library where you usually keep these things
  (typically `/usr/include` and `/usr/lib`.

If this does not work, it *might* be that you have to link `libembed` against the
Flint version used by Nemo.
