module uxo_mode_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule

using LinearAlgebra
using Arpack
using CSV

function uxo_mode_esnice_t4()
    E = 70000 * phun("MPa")
    nu = 0.33
    rho = 2700 * phun("KG/M^3")
    radius = 0.5 * phun("ft")
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0 * 2 * pi)^2

    MR = DeforModelRed3D
    output = import_NASTRAN(joinpath(@__DIR__, "UXO.nas"))
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    @show count(fens)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d, v, nev, nconv = eigs(Symmetric(K + OmegaShift * M),
        Symmetric(M);
        nev = neigvs,
        which = :SM,
        explicittransform = :none)
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    println("Eigenvalues: $fs [Hz]")

    true
end # uxo_modes

# function uxo_mode_swept_h8()

#     E = 70000*phun("MPa");
#     nu = 0.33;
#     rho = 2700*phun("KG/M^3");
#     radius = 0.5*phun("ft");
#     neigvs = 20                   # how many eigenvalues
#     OmegaShift = (1.0*2*pi)^2;

#     MR = DeforModelRed3D
#     xyzrows = CSV.File(joinpath(@__DIR__, "UXO-swept-mesh-xyz.csv"), header=0)
#     xyz = fill(0.0, length(xyzrows), 3)
#     for i in 1:size(xyz, 1)
#         xyz[i, :] .= xyzrows[i]
#     end
#     connrows = CSV.File(joinpath(@__DIR__, "UXO-swept-mesh-conn.csv"), header=0)
#     conn = fill(0, length(connrows), 8)
#     for i in 1:size(conn, 1)
#         conn[i, :] .= connrows[i]
#     end
#     fens = FENodeSet(xyz)
#     fes = FESetH8(conn)

#     # fens.xyz .*= phun("mm") # The input is provided in SI(mm) units

#     geom = NodalField(fens.xyz)
#     u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

#     numberdofs!(u)

#     material = MatDeforElastIso(MR, rho, E, nu, 0.0)

#     femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
#     # femm = FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)
#     associategeometry!(femm,  geom)
#     K  = stiffness(femm, geom, u)
#     M = mass(femm, geom, u)
#     d,v,nev = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
#     d = d .- OmegaShift;
#     fs = real(sqrt.(complex(d)))/(2*pi)
#     println("Eigenvalues: $fs [Hz]")

#     File =  "uxo_mode_swept_h8.vtk"
#     vtkexportmesh(File, fens, fes)

#     for mode = 1:7
#         scattersysvec!(u, v[:,mode])
#         File =  "uxo_mode-$(mode).vtk"
#         vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
#     end

#     true

# end # uxo_modes

function allrun()
    println("#####################################################")
    println("# uxo_mode_nice_t4 omitted due to significant size")
    # uxo_mode_esnice_t4()
    # println("#####################################################")
    # println("# uxo_mode_esnice_t4 ")
    # uxo_mode_swept_h8()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module uxo_mode_examples
nothing
