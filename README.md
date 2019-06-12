[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtoolsDeforLinear.jl/latest/

# FinEtoolsDeforLinear: Linear stress analysis application

`FinEtools` is a package for basic operations on finite element meshes.
`FinEtoolsDeforLinear` is a package using `FinEtools` to solve linear stress analysis problems.
Included is statics and dynamics (modal analysis, steady-state vibration).

## News

- 06/11/2019: Applications are now separated  out from the `FinEtools` package.

[Past news](oldnews.md)

## How to run

The [FinEtools](https://github.com/PetrKryslUCSD/FinEtools.jl) package is
needed. The entire setup of `FinEtoolsDeforLinear` can be performed with
```julia
] activate .; instantiate
```

The package `FinEtoolsDeforLinear` can be tested as
```julia
] activate .; instantiate; test
```

There are a number of examples covering statics and dynamics. The examples may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
