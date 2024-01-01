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

## How to test the package

Here is a record of a session to install this package and test it. You should
see something similar. The git bash running on Windows 10 was used in this
example.

Clone the repo:
```
$ git clone https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git
Cloning into 'FinEtoolsDeforLinear.jl'...
...
Resolving deltas: 100% (87/87), done.
```
Change your working directory, and run Julia:
```
PetrKrysl@Spectre MINGW64 /tmp/exp
$ cd FinEtoolsDeforLinear.jl/

PetrKrysl@Spectre MINGW64 /tmp/exp/FinEtoolsDeforLinear.jl (master)
$ ~/AppData/Local/Julia-1.2.0-rc1/bin/julia.exe
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.2.0-rc1.0 (2019-05-30)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |
```
Activate and instantiate the environment:
```
(v1.2) pkg> activate .; instantiate
[ Info: activating environment ...
```
Test the package:
```
(FinEtoolsDeforLinear) pkg> test
   Testing FinEtoolsDeforLinear
 Resolving package versions...
Test Summary:        | Pass  Total
...
   Testing FinEtoolsDeforLinear tests passed
```

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
