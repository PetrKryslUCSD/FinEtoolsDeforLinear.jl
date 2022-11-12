module test_alum_cyl_mode_examples

using FinEtools
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule: ssit
using LinearAlgebra
using Arpack
using Test

# Mesh alum_cyl.inp
# Abaqus with the Standard isoparametric C3D4 tetrahedron:
C3D4 = [1 0
2 0
3 6.63586E-005
4 0.000171053
5 0.000211299
6 0.000244378
7 2564.63
8 2568.09
9 2597.26
10 4094.38
11 4714.36
12 4717.19
13 5181.98
14 6865.13
15 6868.17
16 6962.86
17 6965.67
18 7024.97
19 7029.44
20 7108.54]
# Abaqus with the standard quadratic tetrahedron:
C3D10 = [1 0
2 0
3 0
4 0.000139365
5 0.000221551
6 0.000291805
7 2546.81
8 2546.81
9 2560.69
10 4100
11 4693.55
12 4693.56
13 5121.57
14 6841.21
15 6841.24
16 6914.22
17 6914.23
18 6950.64
19 6950.66
20 7000.64]

neigvs = 24

function alum_cyl_mode_esnice_t4()

    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    OmegaShift = (10.0*2*pi)^2;
    
    MR = DeforModelRed3D
    output = import_ABAQUS(joinpath(@__DIR__, "alum_cyl.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM, explicittransform=:none)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    # @show     v' * M * v
    @test norm(fs - [0.00000e+00, 0.00000e+00, 0.00000e+00, 5.54160e-06, 8.64750e-05, 1.18749e-04, 2.49815e+03, 2.49888e+03, 2.51331e+03, 4.08265e+03, 4.58599e+03, 4.58642e+03, 4.98701e+03, 6.64802e+03, 6.64848e+03, 6.67904e+03, 6.68216e+03, 6.77789e+03, 6.78059e+03, 6.79936e+03, 6.80400e+03, 7.38167e+03, 7.45600e+03, 7.47771e+03]) < 0.01
    

    true
    
end # alum_cyl_modes

function alum_cyl_mode_esnice_t4_ssit()

    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    OmegaShift = (10.0*2*pi)^2;
    
    MR = DeforModelRed3D
    output = import_ABAQUS(joinpath(@__DIR__, "alum_cyl.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    p = 2 * neigvs
    v0 = [i==j ? one(FFlt) : zero(FFlt) for i=1:size(K,1), j=1:p]
    d, v, nconv, niter, lamberr = ssit(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, 
        v0 = v0, tol = 1.0e-7, maxiter = 301, verbose=false)
    d = d[1:neigvs]
    v = v[:, 1:neigvs]
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    # @show     v' * M * v1
    @test norm(fs - [0.00000e+00, 0.00000e+00, 0.00000e+00, 5.54160e-06, 8.64750e-05, 1.18749e-04, 2.49815e+03, 2.49888e+03, 2.51331e+03, 4.08265e+03, 4.58599e+03, 4.58642e+03, 4.98701e+03, 6.64802e+03, 6.64848e+03, 6.67904e+03, 6.68216e+03, 6.77789e+03, 6.78059e+03, 6.79936e+03, 6.80400e+03, 7.38167e+03, 7.45600e+03, 7.47771e+03]) < 0.01
    true
    
end # alum_cyl_modes


function alum_cyl_mode_esnice_t4_ssit2()

    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    OmegaShift = (10.0*2*pi)^2;
    
    MR = DeforModelRed3D
    output = import_ABAQUS(joinpath(@__DIR__, "alum_cyl.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    
    numberdofs!(u)
    
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    
    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d, v, nconv, niter, lamberr = ssit(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, 
        tol = 1.0e-7, maxiter = 301, verbose=false)
    d = d[1:neigvs]
    v = v[:, 1:neigvs]
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    # @show     v' * M * v1
    @test norm(fs - [0.00000e+00, 0.00000e+00, 0.00000e+00, 5.54160e-06, 8.64750e-05, 1.18749e-04, 2.49815e+03, 2.49888e+03, 2.51331e+03, 4.08265e+03, 4.58599e+03, 4.58642e+03, 4.98701e+03, 6.64802e+03, 6.64848e+03, 6.67904e+03, 6.68216e+03, 6.77789e+03, 6.78059e+03, 6.79936e+03, 6.80400e+03, 7.38167e+03, 7.45600e+03, 7.47771e+03]) < 0.01
    
    
    true
    
end # alum_cyl_modes


alum_cyl_mode_esnice_t4()
alum_cyl_mode_esnice_t4_ssit()
alum_cyl_mode_esnice_t4_ssit2()


end # module 
nothing