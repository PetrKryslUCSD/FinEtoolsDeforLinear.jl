# # Moving sphere in an infinite fluid

# ## Description

# A rigid sphere in an infinite volume of fluid accelerates alternately in the
# positive and negative x-direction, generating positive pressure ahead of it,
# negative pressure behind. Time-dependent simulation. Described in [1].

# ## References

# [1] Krysl P, Hawkins AD, Schilt C, Cranford TW (2012): Angular Oscillation of
# Solid Scatterers in Response to Progressive Planar Acoustic Waves: Do Fish
# Otoliths Rock?. PLOS ONE 7(8): e42591. https://doi.org/10.1371/journal.pone.0042591

# ![](sphere_dipole.png)

# ## Goals

# - Show how to generate hexahedral mesh, mirroring and merging together parts.
# - Execute transient simulation by the trapezoidal-rule time stepping of [1].

##
# ## Definitions

# # Vibration example  solved with FinEtools and Abaqus

# In this example we solve for the free-vibration modes of unit cube  of almost incompressible material.

# The solution with the `FinEtools` package is compared with a commercial software  solution, and hence we also export the model to Abaqus.

# ## Reference:
# Puso MA, Solberg J (2006) A stabilized nodally integrated tetrahedral. International Journal for Numerical Methods in Engineering 67: 841-867.

# This is the finite element toolkit itself.
using FinEtools

# The linear stress analysis application is implemented in this package.
using FinEtoolsDeforLinear

using FinEtoolsDeforLinear.AlgoDeforLinearModule

# Convenience import.
using FinEtools.MeshExportModule

# The eigenvalue problem is solved with the Lanczos algorithm from this package.
using Arpack

# The material properties and dimensions are defined with physical units.
E = 1*phun("PA");
nu = 0.499;
rho = 1*phun("KG/M^3");
a = 1*phun("M"); # length of the side of the cube

# We generate a mesh of  5 x 5 x 5 serendipity 20-node hexahedral elements in a regular grid.
fens,fes  = H20block(a, a, a, 5, 5, 5);

# The problem is solved in three dimensions and hence we create the  displacement field as three-dimensional with three displacement components per node. The degrees of freedom are then numbered  (note that no essential boundary conditions are applied since the cube is free-floating).
geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
numberdofs!(u);

# The model is fully three-dimensional, and hence the material model  and the FEMM created below need to refer to an appropriate model-reduction scheme.
MR = DeforModelRed3D
material = MatDeforElastIso(MR, rho, E, nu, 0.0);

# Note that we compute the stiffness  and the mass matrix using different FEMMs. The difference  is only the quadrature rule chosen: in order to make the mass matrix  non-singular, the accurate  Gauss rule  needs to be used, whereas for the stiffness matrix we want to avoid the excessive stiffness  and therefore  the reduced Gauss rule is used.
femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)), material);
K = stiffness(femm, geom, u)
femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
M = mass(femm, geom, u);

# The free vibration problem  can now be solved.   In order for the eigenvalue solver  to work well, we apply mass-shifting (otherwise the first matrix given to the solver would be singular). We specify the number of eigenvalues to solve for, and we  guess the frequency  with which to shift as 0.01 Hz.
neigvs = 20 # how many eigenvalues
OmegaShift = (0.01*2*pi)^2; # The frequency with which to shift

# The `eigs` routine can now be invoked to solve for a given number of frequencies from the smallest-magnitude end of the spectrum. Note that the mass shifting  needs to be undone when the solution is obtained.
d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
d = d .- OmegaShift;
fs = real(sqrt.(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")
# The first nonzero frequency, frequency 7, should be around .263 Hz.

# The computed mode can be visualized in Paraview. Use the  "Animation view" to produce moving pictures for the mode.
mode = 7
scattersysvec!(u, v[:,mode])
File =  "unit_cube_modes.vtk"
vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
@async run(`"paraview.exe" $File`);

# Finally  we export the model to Abaqus.  Note that we specify the mass density (necessary for dynamics).
AE = AbaqusExporter("unit_cube_modes_h20");
HEADING(AE, "Vibration modes of unit cube  of almost incompressible material.");
COMMENT(AE, "The  first six frequencies are rigid body modes.");
COMMENT(AE, "The  first nonzero frequency (7) should be around 0.26 Hz");
PART(AE, "part1");
END_PART(AE);
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
COMMENT(AE, "The hybrid form of the serendipity hexahedron is chosen because");
COMMENT(AE, "the material is  nearly incompressible.");
ELEMENT(AE, "C3D20RH", "AllElements", 1, connasarray(fes))
ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
END_INSTANCE(AE);
END_ASSEMBLY(AE);
MATERIAL(AE, "elasticity")
ELASTIC(AE, E, nu)
DENSITY(AE, rho)
STEP_FREQUENCY(AE, neigvs)
END_STEP(AE)
close(AE)

# It remains is to load the model into Abaqus and execute it as a job. Alternatively Abaqus can be called on the input file to carry out the analysis at the command line as
# ```
# abaqus job=unit_cube_modes_h20.inp
# ```
# The output database `unit_cube_modes_h20.odb` can then be loaded for postprocessing, for instance from the command line as
# ```
# abaqus viewer database=unit_cube_modes_h20.odb
# ```
# Don't forget to compare the computed frequencies and the mode shapes.  For instance, the first six frequencies should be nearly 0, and the seventh frequency should be approximately  0.262 Hz. There may be  very minor differences due to the fact that  the
# FinEtools formulation is purely displacement-based, whereas the Abaqus model is hybrid (displacement plus pressure).
