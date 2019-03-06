# LatticeGFH90.jl
[![Build Status](https://travis-ci.org/erou/LatticeGFH90.jl.svg?branch=master)](https://travis-ci.org/erou/LatticeGFH90.jl)[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/erou/gfh90-binder/master?filepath=GFH90-demo.ipynb)

A [Julia](https://julialang.org/)/[Nemo](http://nemocas.org/) package to make lattices of compatibly embedded finite fields
based on Hilbert's theorem 90 and Allombert's algorithm.

## How to: install

### The Julia packages

- Install Julia (preferably the last version, to avoid suprises...)
- Type `]` to go in package mode (`Backspace` to quit)

```
(v1.1) pkg> add Nemo

(v1.1) pkg> add "https://github.com/erou/LatticeGFH90.jl"

julia> using Nemo, LatticeGFH90
```

- You're done!

**Note:** this has not been tested on a wide variety of systems, if you
encounter a problem during install, please tell us!

## Testing the package

A longer demo, with explanations, is available on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/erou/gfh90-binder/master?filepath=GFH90-demo.ipynb)

```
julia> make_zetas_conway(5)

julia> k3, x3 = FiniteField(5, 3, "x3")
(Finite field of degree 3 over F_5, x3)

julia> k6, x6 = FiniteField(5, 6, "x6")
(Finite field of degree 6 over F_5, x6)

julia> k12, x12 = FiniteField(5, 12, "x12")
(Finite field of degree 12 over F_5, x12)

julia> f = embed(k3, k6);

julia> g = embed(k6, k12);

julia> h = embed(k3, k12);

julia> x = rand(k3); g(f(x)) == h(x)
true
```
