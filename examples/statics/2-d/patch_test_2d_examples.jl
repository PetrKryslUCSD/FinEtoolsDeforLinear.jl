module patch_test_2d_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
import LinearAlgebra: cholesky

function q4_stress()
    println("Q4. Plane stress.")

    E = 1.0
    nu = 1.0 / 3
    alpha, beta, gamma, delta, eta, phi = 1.0 / 30,
    1.0 / 34,
    -1.0 / 21,
    -1.0 / 51,
    -1.0 / 26,
    -1.0 / 35
    ux(x, y) = alpha + beta * x + gamma * y
    uy(x, y) = delta + eta * x + phi * y
    MR = DeforModelRed2DStress

    fens = FENodeSet([1.0 -0.3;
        2.3 -0.3;
        2.3 0.95;
        1.0 0.95;
        1.4 0.05;
        1.9 -0.03;
        1.7 0.5;
        1.3 0.6])
    fes = FESetQ4([1 2 6 5; 6 2 3 7; 7 3 4 8; 8 4 1 5; 5 6 7 8])

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply prescribed displacements to exterior nodes
    for i in 1:4
        setebc!(u, [i], 1, ux(fens.xyz[i, :]...))
        setebc!(u, [i], 2, uy(fens.xyz[i, :]...))
    end

    applyebc!(u)
    numberdofs!(u)

    material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material)

    F = nzebcloadsstiffness(femm, geom, u)
    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K \ (F)
    scattersysvec!(u, U[:])

    for i in 5:8
        uexact = [ux(fens.xyz[i, :]...), uy(fens.xyz[i, :]...)]
        println("u.values[$i, :] = $(u.values[i, :]), uexact = [$(uexact)]")
    end

    File = "a.vtk"
    vtkexportmesh(File, fes.conn, geom.values,
        FinEtools.MeshExportModule.Q4; vectors = [("u", u.values)])

    true
end # cookstress

function q4_stress_export()
    println("Q4. Plane stress.")

    E = 1.0
    nu = 1.0 / 3
    alpha, beta, gamma, delta, eta, phi = 1.0 / 30,
    1.0 / 34,
    -1.0 / 21,
    -1.0 / 51,
    -1.0 / 26,
    -1.0 / 35
    ux(x, y) = alpha + beta * x + gamma * y
    uy(x, y) = delta + eta * x + phi * y
    MR = DeforModelRed2DStress

    fens = FENodeSet([1.0 -0.3;
        2.3 -0.3;
        2.3 0.95;
        1.0 0.95;
        1.4 0.05;
        1.9 -0.03;
        1.7 0.5;
        1.3 0.6])
    fes = FESetQ4([1 2 6 5; 6 2 3 7; 7 3 4 8; 8 4 1 5; 5 6 7 8])

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply prescribed displacements to exterior nodes
    for i in 1:4
        setebc!(u, [i], 1, ux(fens.xyz[i, :]...))
        setebc!(u, [i], 2, uy(fens.xyz[i, :]...))
    end

    applyebc!(u)
    numberdofs!(u)

    material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material)

    F = nzebcloadsstiffness(femm, geom, u)
    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K \ (F)
    scattersysvec!(u, U[:])

    for i in 5:8
        uexact = [ux(fens.xyz[i, :]...), uy(fens.xyz[i, :]...)]
        println("u.values[$i, :] = $(u.values[i, :]), uexact = [$(uexact)]")
    end

    AE = AbaqusExporter("q4_stress_export")
    HEADING(AE, "q4_stress_export")
    COMMENT(AE, "")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "CPS4", "AllElements", connasarray(fes))
    NSET_NSET(AE, "clamped", 1:4)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", 1.0)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE,
        "ASSEM1.INSTNC1",
        1:4,
        fill(true, 4, 2),
        [[ux(fens.xyz[i, :]...) for i in 1:4] [uy(fens.xyz[i, :]...) for i in 1:4]])
    END_STEP(AE)
    close(AE)

    true
end # cookstress

function allrun()
    println("#####################################################")
    println("# q4_stress ")
    q4_stress()
    println("#####################################################")
    println("# q4_stress_export ")
    q4_stress_export()
    return true
end # function allrun

end # module patch_test_2d_examples
