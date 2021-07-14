module m4thsymmcheck1
using FinEtools
using FinEtoolsDeforLinear: tens4checksymmetry, tens4symmt6x6tot!, tens4symmtto6x6t!
using LinearAlgebra
using Test
function test()
	C = rand(6, 6)
	C = C + C'
	Co = deepcopy(C)
	t = fill(0.0, 3, 3, 3, 3)
	tens4symmt6x6tot!(t, C)
	@test tens4checksymmetry(t)
	tens4symmtto6x6t!(C, t)
	@assert norm(C - C') <= eps(1.0)
	@assert norm(C - Co) <= eps(1.0)
end
end
using .m4thsymmcheck1
m4thsymmcheck1.test()




module mlumpedmass1
using FinEtools
using FinEtoolsDeforLinear
using Test
import Arpack: eigs
function test()
    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 10;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (0.01*2*pi)^2;

    MR = DeforModelRed3D
    fens,fes  = H20block(a,b,h, na,nb,nh)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material=MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)), material)
    K = stiffness(femm, geom, u)
    
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    @test abs(fs[7]-0.26259869196259) < 1.0e-5

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
    M = mass(femm, SysmatAssemblerSparseHRZLumpingSymm(), geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    @test abs(fs[7]- 0.2598164590380608) < 1.0e-5
    # println("Eigenvalues: $fs [Hz]")

    # mode = 17
    # scattersysvec!(u, v[:,mode])
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])

    
    true
end
end
using .mlumpedmass1
mlumpedmass1.test()

