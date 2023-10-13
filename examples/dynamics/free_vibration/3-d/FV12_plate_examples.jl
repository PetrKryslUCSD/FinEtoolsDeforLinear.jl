module FV12_plate_examples
using Statistics
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using SubSIt
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
# using PGFPlotsX
using Test

E = 200e3*phun("MPa");
nu = 0.3;
rho = 8000*phun("KG/M^3");
L = 10.00*phun("M"); t = 0.05*phun("M");
nL = 8; nt = 4;
neigvs = 14                   # how many eigenvalues
OmegaShift = (1.0*2*pi)^2;
# Fundamental frequency
f_analytical =[0 0 0  0 0 0 1.622 2.360 2.922 4.190 4.190 7.356 7.356 7.668];

function FV12_plate_esnice()
    global a,b,h,na,nb,nh
    MR = DeforModelRed3D
    fens, fes = T4block(L,L,t,nL,nL,nt);

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)

    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    println("f/f_analytical % = $(fs[7:10]  ./ f_analytical[7:10] .* 100) %")

    vectors = []
    for i = 7:length(fs)
        scattersysvec!(u, v[:, i])
        push!(vectors, ("Mode_$i", deepcopy(u.values)))
    end
    File  =   "rectangular_plate_esnice.vtk"
    vtkexportmesh(File, connasarray(fes), fens.xyz, FinEtools.MeshExportModule.VTK.T4; vectors = vectors)
    @async run(`"paraview.exe" $File`)

    d,v,nconv = SubSIt.ssit(K+OmegaShift*M, M; nev=neigvs)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")

    true
end # function

function allrun()
    println("#####################################################")
    println("# FV12_plate_esnice ")
    FV12_plate_esnice()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing
    
