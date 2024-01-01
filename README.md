[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl/actions)
[![Code Coverage](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsDeforLinear.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PetrKryslUCSD/FinEtoolsDeforLinear.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsDeforLinear.jl/latest)
[![Codebase Graph](https://img.shields.io/badge/Codebase-graph-green.svg)](https://octo-repo-visualization.vercel.app/?repo=PetrKryslUCSD/FinEtoolsDeforLinear.jl)

# FinEtoolsDeforLinear: Linear stress analysis application


[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl.git) is a package
for basic operations on finite element meshes. `FinEtoolsDeforLinear` is a
package using `FinEtools` to solve linear stress analysis problems. Included is
statics and dynamics (modal analysis, steady-state vibration).

## News

- 12/31/2023: Update for Julia 1.10.
- 12/22/2023: Merge the tutorials into the package tree.
- 10/23/2023: Remove dependency on FinEtools predefined types (except for the data dictionary in the algorithm module).
- 06/21/2023: Update for FinEtools 7.0.
- 08/15/2022: Updated all examples.
- 08/09/2022: Updated all 3D dynamic examples.
- 04/03/2022: Examples now have their own project environment.

[Past news](oldnews.md)

## Tutorials

There are a number of tutorials explaining the use of this package.
Check out the [index](https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl/blob/main/tutorials/index.md). The  tutorials themselves can be executed as
follows:

- Download the package or clone it.
```
git clone https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git
```
- Change into the `tutorials` folder: `cd .\FinEtoolsDeforLinear.jl\tutorials`.
- Start Julia: `julia`.
- Activate the environment:
```
using Pkg; Pkg.activate("."); Pkg.instantiate();
```
- Execute the desired tutorial. Here `name.jl` is the name of the tutorial file:
```
include("name.jl")
```

## Examples

Many examples of solving for static and dynamic stress response with continuum FE models are available.
Begin with changing your working directory to the `examples` folder. Activate
and instantiate the examples environment.
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
There are a number of examples covering statics and dynamics. The examples may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
