"""
FinEtoolsDeforLinear (C) 2017-2019, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics. Package for linear stress analysis problems.
"""
module FinEtoolsDeforLinear

__precompile__(true)

include("allmodules.jl")

# Exports follow:

###########################################################################
# Linear deformation functionality
###########################################################################
using .DeforModelRedModule: AbstractDeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,    DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D, nstressstrain, nthermstrain, stresscomponentmap, Blmat!, divmat, vgradmat
# Exported: types  for model reduction in stress analysis
export AbstractDeforModelRed, DeforModelRed1D, DeforModelRed2DStrain,    DeforModelRed2DStress, DeforModelRed2DAxisymm, DeforModelRed3D
# Exported: num stresses/strains,  number of thermal strains, and map of  the numbering of stress components
export nstressstrain, nthermstrain, stresscomponentmap
# Exported: strain-displacement matrix and divergence matrix for all model-reduction types
export Blmat!, divmat, vgradmat

using .MatDeforModule: AbstractMatDefor
# using .MatDeforModule: AbstractMatDefor, strain2x2tto3v!, strain3vto2x2t!, strain3x3tto6v!, strain6vto3x3t!, strain9vto3x3t!, strain3x3tto9v!, strain9vto6v!, strain6vto9v!, stress2x2to3v!,  stress3vto2x2t!, stress3vto3x3t!, stress4vto3x3t!, stress6vto3x3t!, stress3x3tto6v!, stress9vto6v!,  stress6vto9v!
using .MatDeforModule: rotstressvec!
using .MatDeforModule: dett, strainvdet, strainvtr, strainttov!, strainvtot!, stressttov!, stressvtot!
using .MatDeforModule: tens4checksymmetry, tens4symmto6x6t!, tens4symm6x6tot!, tens4dot2!, tens4ijkl!, tens4ikjl!, tens4iljk!, tens4identity!, tens4transposor!, tens4tracor!, tens4symmetrizor!, tens4skewor!, tens4deviator!
# Exported: abstract type for  models of deformation,  conversion methods  for strain and stress, transformations  of strain and stress
export AbstractMatDefor
export rotstressvec!
export dett, strainvdet, strainvtr, strainttov!, strainvtot!, stressttov!, stressvtot!
export tens4checksymmetry, tens4symmto6x6t!, tens4symm6x6tot!, tens4dot2!, tens4ijkl!, tens4ikjl!, tens4iljk!, tens4identity!, tens4transposor!, tens4tracor!, tens4symmetrizor!, tens4skewor!, tens4deviator!

using .MatDeforLinearElasticModule: AbstractMatDeforLinearElastic
# Exported: type of  isotropic elastic material
export AbstractMatDeforLinearElastic

using .MatDeforElastIsoModule: MatDeforElastIso
# Exported: type of  isotropic elastic material
export MatDeforElastIso

using .MatDeforElastOrthoModule: MatDeforElastOrtho
# Exported: type of orthotropic elastic material
export MatDeforElastOrtho

using .FEMMDeforLinearBaseModule: AbstractFEMMDeforLinear, stiffness, nzebcloadsstiffness, thermalstrainloads, mass, inspectintegpoints
# Exported: abstract type for linear information, discretization methods for the abstract type
export AbstractFEMMDeforLinear, stiffness, nzebcloadsstiffness, thermalstrainloads, mass, inspectintegpoints

using .FEMMDeforLinearModule: FEMMDeforLinear
# Exported: type for linear deformation
export FEMMDeforLinear

using .FEMMDeforWinklerModule: FEMMDeforWinkler, surfacenormalspringstiffness
# Exported: type for distributed-spring support, discretization method
export FEMMDeforWinkler, surfacenormalspringstiffness

using .FEMMDeforLinearMSModule: FEMMDeforLinearMSH8, FEMMDeforLinearMST10, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
# Exported: type for mean-strain solid elements, discretization methods
export FEMMDeforLinearMSH8, FEMMDeforLinearMST10, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

using .FEMMDeforSurfaceDampingModule: FEMMDeforSurfaceDamping, dampingABC
#Exported: type for surface damping (absorbing boundary conditions)
export FEMMDeforSurfaceDamping, dampingABC

using .FEMMDeforLinearNICEModule: FEMMDeforLinearNICEH8, FEMMDeforLinearNICET4, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
# Exported: type for NICE (Nodally-integrated continuum elements) solid elements, discretization methods
export FEMMDeforLinearNICEH8, FEMMDeforLinearNICET4, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

using .FEMMDeforLinearESNICEModule: FEMMDeforLinearESNICET4, FEMMDeforLinearESNICEH8, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints
# Exported: type for ESICE (Energy-sampling stabilized nodally-integrated continuum elements) solid elements, discretization methods
export FEMMDeforLinearESNICET4, FEMMDeforLinearESNICEH8, stiffness, nzebcloadsstiffness, thermalstrainloads, inspectintegpoints

end # module
