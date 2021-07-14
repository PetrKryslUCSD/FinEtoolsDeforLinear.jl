module mbaruniax1
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 30e6*phun("psi") # Young's modulus
    alpha = 5.5e-6*phun("in")/phun("in")/phun("F")
    for nu in [0.0 0.3 0.45]

        MR = DeforModelRed1DStress
        material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
        D = fill(0.0, 1, 1)
        t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
        tangentmoduli!(material,  D,  t, dt, loc, label)
        @test D[1, 1] ≈ E
    end
    
    true
end
end
using .mbaruniax1
mbaruniax1.test()

module mbaruniax2
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear.MatDeforLinearElasticModule: tangentmoduli!, update!, thermalstrain!
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 30e6*phun("psi") # Young's modulus
    alpha = 5.5e-6*phun("in")/phun("in")/phun("F")
    for nu in [0.0 0.3 0.45]

        MR = DeforModelRed1DStrain
        material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
        D = fill(0.0, 1, 1)
        t::FFlt, dt::FFlt, loc::FFltMat, label::FInt = 0.0, 0.0, [0.0 0.0 0.0], 0
        tangentmoduli!(material,  D,  t, dt, loc, label)
        lambda = E * nu / (1 + nu) / (1 - 2*(nu));
        mu = E / 2. / (1+nu);
        @test D[1, 1] ≈ lambda + 2*mu
    end
    
    true
end
end
using .mbaruniax2
mbaruniax2.test()
