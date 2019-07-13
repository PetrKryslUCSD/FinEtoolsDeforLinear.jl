# Guide

The
[`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest/index.html)
package is used here to solve linear stress analysis (deformation)
problems.

## Modules

The package `FinEtoolsDeforLinear` has the following structure:

- `FinEtoolsDeforLinear` is the  top-level module.
- Linear deformation:  `AlgoDeforLinearModule` (algorithms),
  `DeforModelRedModule` (model-reduction definitions, 3D, plane strain
  and stress, and so on), `FEMMDeforLinearBaseModule`,
  `FEMMDeforLinearModule`, `FEMMDeforLinearMSModule`,
  `FEMMDeforWinklerModule` (FEM machines to evaluate the matrix and
  vector quantities),  `MatDeforModule`, `MatDeforElastIsoModule`,
  `MatDeforElastOrthoModule` (elastic material models).


## Linear deformation FEM  machines

For  the base machine for linear deformation, `FEMMDeforLinearBase`,
assumes standard isoparametric  finite elements. It evaluates  the
interior integrals:

- The stiffness matrix, the mass matrix.

- The load vector corresponding to thermal strains.

Additionally:

- Function to inspect  integration points.

The FEM machine `FEMMDeforLinear` simply stores the data required by the
base `FEMMDeforLinearBase`.

The machine `FEMMDeforWinkler` is specialized for the boundary integrals
for bodies  supported  on continuously distributed springs:

- Compute the stiffness matrix corresponding to the springs.

The  mean-strain FEM machine `FEMMDeforLinearMS` implements advanced
hexahedral and tetrahedral elements based on multi-field theory and
energy-sampling  stabilization. It provides functions to compute:

- The stiffness matrix, the mass matrix.

- The load vector corresponding to thermal strains.

Additionally it defines:

- Function to inspect  integration points.

## Materials for linear deformation analysis

The module `MatDeforModule` provides functions to convert between vector
and matrix (tensor) representations of stress and strain. Further,
functions to rotate stress and strain between different coordinate
systems (based upon the model-reduction type, 3-D, 2-D, or 1-D) are
provided.

Currently  there are material types for isotropic and orthotropic linear
elastic materials. The user may add  additional material types by
deriving from `AbstractMatDefor` and equipping them with three methods:
(1) compute the tangent moduli, (2) update the material state, (3)
compute the thermal strain.

For full generality, material types  should implement these methods for
fully three-dimensional, plane strain and plane stress, 2D axially
symmetric, and one-dimensional deformation models.

## Linear deformation algorithms

There are algorithms for

- Linear static analysis;
- Export  of the deformed shape for visualization;
- Export  of the nodal and elementwise stress fields for visualization;
- Modal (free-vibration) analysis;
- Export  of modal shapes for visualization;
- Subspace-iteration method implementation.

### Model data

Model data is a dictionary, with string keys, and arbitrary values.
The documentation string for each method of an algorithm lists the required input.
For instance, for the method `linearstatics` of the `AlgoDeforLinearModule`, the
`modeldata` dictionary needs to provide key-value pairs for the finite element node set, and
the regions, the boundary conditions, and so on.

The `modeldata` may be also supplemented with additional key-value pairs inside an algorithm
and returned for further processing by other algorithms.
