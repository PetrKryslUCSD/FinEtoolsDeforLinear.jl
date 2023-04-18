
module malum_cyl_mode_esnice_t4
using FinEtools
using FinEtoolsDeforLinear
using Test
using Arpack
using LinearAlgebra
using DataDrop
using InteractiveUtils
function test()
    # Aluminum cylinder free vibration, mesh imported from Abaqus
    # Mesh converted from quadratic tetrahedra to linear tetrahedra
    # NICE tetrahedral elements used
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    Eigenvalues =   [0.0, 0.0, 0.0, 1.8846e-5, 7.35917e-5, 0.000119445, 2498.15, 2498.88, 2513.31, 4082.65, 4585.99, 4586.42, 4987.01, 6648.02, 6648.48, 6679.04, 6682.16, 6777.89, 6780.59, 6799.36]

    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    @show boundingbox(fens.xyz)
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    @edit stiffness(femm, SysmatAssemblerSparse(), geom, u)
    K  = stiffness(femm, SysmatAssemblerSparse(), geom, u)
    M = mass(femm, geom, u)

    Kref = DataDrop.retrieve_matrix("rfK")
    @show norm(Kref - K)

    Mref = DataDrop.retrieve_matrix("rfM")
    @show norm(Mref - M)

    DataDrop.empty_hdf5_file("K")
    DataDrop.empty_hdf5_file("M")
    DataDrop.store_matrix("K", K)
    DataDrop.store_matrix("M", M)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    println("Eigenvalues: $fs [Hz]")
    @show norm(vec(fs) .- vec(Eigenvalues)) < 1.0e-3*maximum(vec(Eigenvalues))

    nothing
end
end
using .malum_cyl_mode_esnice_t4
malum_cyl_mode_esnice_t4.test()

