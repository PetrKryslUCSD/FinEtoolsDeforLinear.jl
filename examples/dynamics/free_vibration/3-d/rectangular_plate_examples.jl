module rectangular_plate_examples
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
using PGFPlotsX
using Test

E = 210e3*phun("MPa");
nu = 0.3;
rho = 7850*phun("KG/M^3");
a=4.00*phun("M"); b=1.00*phun("M"); h= 0.1*phun("M");
na= 10; nb=  5; nh =2;
na= 2*10; nb=  2*5; nh =8;
neigvs = 20                   # how many eigenvalues
OmegaShift = (10.0*2*pi)^2;
# Fundamental frequency
f_analytical =[0 0 0  0 0 0 33.78 82.28  92.99 170.06];

function rectangular_plate_esnice()
    global a,b,h,na,nb,nh
    MR = DeforModelRed3D
    fens, fes = T4block(a,b,h,na,nb,nh);

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    @show minimum(vec(femm.nphis)), maximum(vec(femm.nphis))
    @pgf a = Axis({
            xlabel = "Entity",
            ylabel = "Stabilization factor",
            grid="major",
            legend_pos  = "north east"
        },
        Plot({"only marks", mark="+"}, Table([:x => vec(1:count(fes)), :y => vec(femm.ephis)]))
        )
    display(a)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    println("f/f_analytical = $(fs[7:10]  ./ f_analytical[7:10] .* 100) %")

    vectors = []
    for i = 7:length(fs)
        scattersysvec!(u, v[:, i])
        push!(vectors, ("Mode_$i", deepcopy(u.values)))
    end
    File  =   "rectangular_plate_esnice.vtk"
    vtkexportmesh(File, connasarray(fes), fens.xyz, FinEtools.MeshExportModule.VTK.T4; vectors = vectors)
    @async run(`"paraview.exe" $File`)


    true
end # function

function allrun()
    println("#####################################################")
    println("# rectangular_plate_esnice ")
    rectangular_plate_esnice()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing
