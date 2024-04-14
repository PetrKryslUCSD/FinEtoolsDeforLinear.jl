module NAFEMS_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using LinearAlgebra
using Arpack

function NAFEMS_FV32_algo()
    println("""
    FV32: Cantilevered tapered membrane
    This is a test recommended by the National Agency for Finite Element Methods and
    Standards (U.K.): Test FV32 from NAFEMS publication TNSB, Rev. 3, “The Standard
    NAFEMS Benchmarks,” October 1990.

    Reference solution: 44.623	130.03	162.70	246.05	379.90	391.44 for the first
    six modes.
    """)

    t0 = time()

    E = 200 * phun("GPA")
    nu = 0.3
    rho = 8000 * phun("KG/M^3")
    L = 10 * phun("M")
    W0 = 5 * phun("M")
    WL = 1 * phun("M")
    nL, nW = 28, 14 # How many element edges per side?
    neigvs = 10                   # how many eigenvalues
    Reffs = [44.623 130.03 162.70 246.05 379.90 391.44]

    fens, fes = Q8block(1.0, 2.0, nL, nW)
    for i = 1:count(fens)
        xi, eta = fens.xyz[i, :]
        eta = eta - 1.0
        fens.xyz[i, :] .= (xi * L, eta * (1.0 - 0.8 * xi) * W0 / 2)
    end
    # File =  "mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    # Make the region
    MR = DeforModelRed2DStress
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict(
        "femm" => FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material),
        "femm_mass" => FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3)), material),
    )

    nl1 = selectnode(fens; plane = [1.0 0.0 0.0], thickness = L / 1.0e4)
    ebc1 = FDataDict("node_list" => nl1, "component" => 1, "displacement" => 0.0)
    ebc2 = FDataDict("node_list" => nl1, "component" => 2, "displacement" => 0.0)

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "regions" => [region1],
        "essential_bcs" => [ebc1 ebc2],
        "neigvs" => neigvs,
    )

    # Solve
    modeldata = AlgoDeforLinearModule.modal(modeldata)

    fs = modeldata["omega"] / (2 * pi)
    println("Frequencies: $(fs[1:6]) [Hz]")
    println("Percentage frequency errors: $((vec(fs[1:6]) - vec(Reffs))./vec(Reffs)*100)")

    modeldata["postprocessing"] = FDataDict("file" => "FV32-modes", "mode" => 1:10)
    modeldata = AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)

    true
end # NAFEMS_FV32_algo

function NAFEMS_FV32_algo_interior()
    println("""
    FV32: Cantilevered tapered membrane
    This is a test recommended by the National Agency for Finite Element Methods and
    Standards (U.K.): Test FV32 from NAFEMS publication TNSB, Rev. 3, “The Standard
    NAFEMS Benchmarks,” October 1990.

    Reference solution: 44.623	130.03	162.70	246.05	379.90	391.44 for the first
    six modes.
    """)

    t0 = time()

    E = 200 * phun("GPA")
    nu = 0.3
    rho = 8000 * phun("KG/M^3")
    L = 10 * phun("M")
    W0 = 5 * phun("M")
    WL = 1 * phun("M")
    nL, nW = 28, 14 # How many element edges per side?
    neigvs = 10                   # how many eigenvalues
    Reffs = [44.623 130.03 162.70 246.05 379.90 391.44]

    fens, fes = Q8block(1.0, 2.0, nL, nW)
    for i = 1:count(fens)
        xi, eta = fens.xyz[i, :]
        eta = eta - 1.0
        fens.xyz[i, :] .= (xi * L, eta * (1.0 - 0.8 * xi) * W0 / 2)
    end
    # File =  "mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    # Make the region
    MR = DeforModelRed2DStress
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict(
        "femm" => FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material),
        "femm_mass" => FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3)), material),
    )

    nl1 = selectnode(fens; plane = [1.0 0.0 0.0], thickness = L / 1.0e4)
    ebc1 = FDataDict("node_list" => nl1, "component" => 1, "displacement" => 0.0)
    ebc2 = FDataDict("node_list" => nl1, "component" => 2, "displacement" => 0.0)
    nl2 = connectednodes(meshboundary(fes))
    ebc3 = FDataDict("node_list" => nl2, "component" => 1, "displacement" => 0.0)
    ebc4 = FDataDict("node_list" => nl2, "component" => 2, "displacement" => 0.0)

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "regions" => [region1],
        "essential_bcs" => [ebc1 ebc2 ebc3 ebc4],
        "neigvs" => neigvs,
    )

    # Solve
    modeldata = AlgoDeforLinearModule.modal(modeldata)

    fs = modeldata["omega"] / (2 * pi)
    println("Frequencies: $(fs[1:6]) [Hz]")
    println("Percentage frequency errors: $((vec(fs[1:6]) - vec(Reffs))./vec(Reffs)*100)")

    modeldata["postprocessing"] = FDataDict("file" => "FV32-modes", "mode" => 1:10)
    modeldata = AlgoDeforLinearModule.exportmode(modeldata)
    @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)

    true
end # NAFEMS_FV32_algo

function allrun()
    println("#####################################################")
    println("# NAFEMS_FV32_algo ")
    NAFEMS_FV32_algo()
    println("#####################################################")
    println("# NAFEMS_FV32_algo_interior ")
    NAFEMS_FV32_algo_interior()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
