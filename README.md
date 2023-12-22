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

- 12/22/2023: Merge the tutorials into the package tree.
- 10/23/2023: Remove dependency on FinEtools predefined types (except for the data dictionary in the algorithm module).
- 06/21/2023: Update for FinEtools 7.0.
- 08/15/2022: Updated all examples.
- 08/09/2022: Updated all 3D dynamic examples.
- 04/03/2022: Examples now have their own project environment.
- 02/08/2021: Updated dependencies for Julia 1.6 and FinEtools 5.0.
- 08/23/2020: Added a separate tutorial package, [FinEtoolsDeforLinearTutorials.jl](https://petrkryslucsd.github.io/FinEtoolsDeforLinearTutorials.jl)).
- 08/17/2020: Added tutorials to the documentation.
- 04/02/2020: The examples still need to be updated, some don't work, sorry.
- 01/23/2020: Dependencies have been updated to work with Julia 1.3.1.


[Past news](oldnews.md)

## How to test the package

Here is a record of a session to install this package and test it. You should
see something similar. The git bash running on Windows 10 was used in this
example.

Clone the repo:
```
$ git clone https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git
Cloning into 'FinEtoolsDeforLinear.jl'...
remote: Enumerating objects: 207, done.
remote: Counting objects: 100% (207/207), done.
remote: Compressing objects: 100% (121/121), done.
remote: Total 207 (delta 87), reused 191 (delta 74), pack-reused 0
Receiving objects: 100% (207/207), 3.28 MiB | 6.71 MiB/s, done.
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
[ Info: activating environment at `C:\Users\PETRKR~1\AppData\Local\Temp\exp\FinEtoolsDeforLinear.jl\Project.toml`.
   Cloning default registries into `C:\Users\PetrKrysl\.julia`
   Cloning registry from "https://github.com/JuliaRegistries/General.git"
     Added registry `General` to `C:\Users\PetrKrysl\.julia\registries\General`
   Cloning git-repo `https://github.com/PetrKryslUCSD/FinEtools.jl.git`
  Updating git-repo `https://github.com/PetrKryslUCSD/FinEtools.jl.git`
 Installed Missings ──────────── v0.4.1
 Installed Requires ──────────── v0.5.2
 Installed Compat ────────────── v2.1.0
 Installed Crayons ───────────── v4.0.0
 Installed OrderedCollections ── v1.1.0
 Installed DefaultApplication ── v0.1.3
 Installed Arpack ────────────── v0.3.1
 Installed BinaryProvider ────── v0.5.4
 Installed SortingAlgorithms ─── v0.3.1
 Installed ArgCheck ──────────── v1.0.1
 Installed DocStringExtensions ─ v0.7.0
 Installed Tokenize ──────────── v0.5.4
 Installed DataStructures ────── v0.15.0
 Installed MacroTools ────────── v0.5.0
 Installed StatsBase ─────────── v0.30.0
 Installed CSTParser ─────────── v0.6.0
 Installed Parameters ────────── v0.10.3
 Installed PGFPlotsX ─────────── v0.3.8
  Building Arpack ───→ `C:\Users\PetrKrysl\.julia\packages\Arpack\cu5By\deps\build.log`
  Building PGFPlotsX → `C:\Users\PetrKrysl\.julia\packages\PGFPlotsX\PZlVQ\deps\build.log`
```
Test the package:
```
(FinEtoolsDeforLinear) pkg> test
   Testing FinEtoolsDeforLinear
 Resolving package versions...
Test Summary:        | Pass  Total
Linear deformation 1 |   12     12
 44.937931 seconds (79.93 M allocations: 6.428 GiB, 6.69% gc time)
Test Summary:        | Pass  Total
Linear deformation 2 |   19     19
 15.568998 seconds (32.58 M allocations: 2.429 GiB, 6.54% gc time)
Test Summary:        | Pass  Total
Linear deformation 3 |   30     30
  5.180744 seconds (10.47 M allocations: 600.374 MiB, 3.84% gc time)
Test Summary:        | Pass  Total
Linear deformation 4 |   31     31
 25.287681 seconds (32.96 M allocations: 4.108 GiB, 6.55% gc time)
Test Summary:        | Pass  Total
Linear deformation 5 |   22     22
  7.057942 seconds (15.94 M allocations: 984.675 MiB, 4.71% gc time)
Test Summary:        | Pass  Total
Linear deformation 6 |   29     29
 23.001916 seconds (57.78 M allocations: 6.001 GiB, 11.74% gc time)
Test Summary:        | Pass  Total
Linear deformation 7 |   36     36
  3.932681 seconds (13.96 M allocations: 799.952 MiB, 6.10% gc time)
   Testing FinEtoolsDeforLinear tests passed
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
