module unit_cube_tet_examples
using FinEtools
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using LinearAlgebra
using Arpack

function unit_cube_esnice()
    println("""
    Vibration modes of unit cube  of almost incompressible material.

    Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    tetrahedral. International Journal for Numerical Methods in
    Engineering 67: 841-867.
    """)
    t0 = time()

    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 16;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (0.01*2*pi)^2;

    MR = DeforModelRed3D
    fens,fes  = T4block(a,b,h, na,nb,nh)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")

    mode = 17
    scattersysvec!(u, v[:,mode])
    File =  "unit_cube_esnice.vtk"
    vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
    @async run(`"paraview.exe" $File`)
    true

end # unit_cube_esnice

function allrun()
    println("#####################################################")
    println("# unit_cube_esnice ")
    unit_cube_esnice()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing