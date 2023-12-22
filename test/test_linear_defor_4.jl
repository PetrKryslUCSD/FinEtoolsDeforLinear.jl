

module mmmCookmstressmorthommm
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


    E = 1.0
    nu = 1.0 / 3
    width = 48.0
    height = 44.0
    thickness = 1.0
    free_height = 16.0
    Mid_edge = [48.0, 52.0]# Location of tracked  deflection
    magn = 1.0 / (free_height * thickness)# Density of applied load
    convutip = 23.97
    n = 30# number of elements per side
    tolerance = minimum([width, height]) / n / 1000.0#Geometrical tolerance

    fens, fes = T3block(width, height, n, n)

    # Reshape into a trapezoidal panel
    for i = 1:count(fens)
        fens.xyz[i, 2] =
            fens.xyz[i, 2] +
            (fens.xyz[i, 1] / width) *
            (height - fens.xyz[i, 2] / height * (height - free_height))
    end

    # Clamped edge of the membrane
    l1 = selectnode(fens; box = [0.0, 0.0, -Inf, Inf], inflate = tolerance)
    ess1 = FDataDict("displacement" => 0.0, "component" => 1, "node_list" => l1)
    ess2 = FDataDict("displacement" => 0.0, "component" => 2, "node_list" => l1)

    # Traction on the opposite edge
    boundaryfes = meshboundary(fes)
    Toplist =
        selectelem(fens, boundaryfes, box = [width, width, -Inf, Inf], inflate = tolerance)
    el1femm =
        FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
    flux1 = FDataDict("traction_vector" => [0.0, +magn], "femm" => el1femm)

    # Make the region
    MR = DeforModelRed2DStress
    material = MatDeforElastOrtho(MR, 0.0, E, nu, 0.0)
    region1 = FDataDict(
        "femm" =>
            FEMMDeforLinear(MR, IntegDomain(fes, TriRule(1), thickness), material),
    )

    modeldata = FDataDict(
        "fens" => fens,
        "regions" => [region1],
        "essential_bcs" => [ess1, ess2],
        "traction_bcs" => [flux1],
    )

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

    u = modeldata["u"]
    geom = modeldata["geom"]

    # Extract the solution
    nl = selectnode(
        fens,
        box = [Mid_edge[1], Mid_edge[1], Mid_edge[2], Mid_edge[2]],
        inflate = tolerance,
    )
    theutip = u.values[nl, :]
    # println("displacement =$(theutip[2]) as compared to converged $convutip")
    @test abs(theutip[2] - 23.79623002934245) < 1.0e-3

    modeldata["postprocessing"] =
        FDataDict("file" => "cookstress-ew", "quantity" => :Cauchy, "component" => :xy)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.07275125002229098, 0.1309644027374448],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] =
        FDataDict("file" => "cookstress-ew-vm", "quantity" => :vm, "component" => 1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [0.014136291979824203, 0.4181381117075297],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] = FDataDict(
        "file" => "cookstress-ew-pressure",
        "quantity" => :pressure,
        "component" => 1,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
    # display(norm([minimum(fld.values), maximum(fld.values)] - [-0.12996180159464202, 0.2436183072544255]))
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.08664120106309468, 0.16241220483628366],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] = FDataDict(
        "file" => "cookstress-ew-princ1",
        "quantity" => :princCauchy,
        "component" => 1,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.10280794467415574, 0.24831137383158813],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] = FDataDict(
        "file" => "cookstress-ew-princ3",
        "quantity" => :princCauchy,
        "component" => 2,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.4398236425818378, 0.022575693063719465],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    AE = AbaqusExporter("Cookstress_algo_stress")
    HEADING(AE, "Cook trapezoidal panel, plane stress")
    COMMENT(AE, "Converged free mid-edge displacement = 23.97")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    COMMENT(AE, "We are assuming three node triangles in plane-stress")
    COMMENT(AE, "CPS3 are pretty poor-accuracy elements, but here we don't care about it.")
    @test nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
    ELEMENT(
        AE,
        "CPS3",
        "AllElements",
        connasarray(modeldata["regions"][1]["femm"].integdomain.fes),
    )
    NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
    bfes = modeldata["traction_bcs"][1]["femm"].integdomain.fes
    COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary")
    COMMENT(AE, "have two nodes each and also that they are the same length.")
    COMMENT(AE, "Then the concentrated loads below will be correctly lumped.")
    nl = connectednodes(bfes)
    F = zeros(count(modeldata["fens"]))
    for ix = 1:count(bfes)
        for jx = 1:2
            F[bfes.conn[ix][jx]] += 1.0 / n / 2 / thickness
        end
    end
    for ixxxx = 1:length(F)
        if F[ixxxx] != 0.0
            CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
        end
    end
    END_STEP(AE)
    close(AE)

    nlines = 0
    open(AE.filename) do f
        s = readlines(f)
        nlines = length(s)
    end
    @test nlines == 2886
    rm(AE.filename)
end
end
using .mmmCookmstressmorthommm
mmmCookmstressmorthommm.test()

module mmmCookmstressisommm
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


    E = 1.0
    nu = 1.0 / 3
    width = 48.0
    height = 44.0
    thickness = 1.0
    free_height = 16.0
    Mid_edge = [48.0, 52.0]# Location of tracked  deflection
    magn = 1.0 / (free_height * thickness)# Density of applied load
    convutip = 23.97
    n = 30# number of elements per side
    tolerance = minimum([width, height]) / n / 1000.0#Geometrical tolerance

    fens, fes = T3block(width, height, n, n)

    # Reshape into a trapezoidal panel
    for i = 1:count(fens)
        fens.xyz[i, 2] =
            fens.xyz[i, 2] +
            (fens.xyz[i, 1] / width) *
            (height - fens.xyz[i, 2] / height * (height - free_height))
    end

    # Clamped edge of the membrane
    l1 = selectnode(fens; box = [0.0, 0.0, -Inf, Inf], inflate = tolerance)
    ess1 = FDataDict("displacement" => 0.0, "component" => 1, "node_list" => l1)
    ess2 = FDataDict("displacement" => 0.0, "component" => 2, "node_list" => l1)

    # Traction on the opposite edge
    boundaryfes = meshboundary(fes)
    Toplist =
        selectelem(fens, boundaryfes, box = [width, width, -Inf, Inf], inflate = tolerance)
    el1femm =
        FEMMBase(IntegDomain(subset(boundaryfes, Toplist), GaussRule(1, 2), thickness))
    flux1 = FDataDict("traction_vector" => [0.0, +magn], "femm" => el1femm)

    # Make the region
    MR = DeforModelRed2DStress
    material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)
    region1 = FDataDict(
        "femm" =>
            FEMMDeforLinear(MR, IntegDomain(fes, TriRule(1), thickness), material),
    )

    modeldata = FDataDict(
        "fens" => fens,
        "regions" => [region1],
        "essential_bcs" => [ess1, ess2],
        "traction_bcs" => [flux1],
    )

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

    u = modeldata["u"]
    geom = modeldata["geom"]

    # Extract the solution
    nl = selectnode(
        fens,
        box = [Mid_edge[1], Mid_edge[1], Mid_edge[2], Mid_edge[2]],
        inflate = tolerance,
    )
    theutip = u.values[nl, :]
    # println("displacement =$(theutip[2]) as compared to converged $convutip")
    @test abs(theutip[2] - 23.79623002934245) < 1.0e-3

    modeldata["postprocessing"] =
        FDataDict("file" => "cookstress-ew", "quantity" => :Cauchy, "component" => :xy)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.07275125002229098, 0.1309644027374448],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] =
        FDataDict("file" => "cookstress-ew-vm", "quantity" => :vm, "component" => 1)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of vm = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [0.014136291979824203, 0.4181381117075297],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] = FDataDict(
        "file" => "cookstress-ew-pressure",
        "quantity" => :pressure,
        "component" => 1,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of pressure = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.08664120106309468, 0.16241220483628366],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] = FDataDict(
        "file" => "cookstress-ew-princ1",
        "quantity" => :princCauchy,
        "component" => 1,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of princCauchy Max = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.10280794467415574, 0.24831137383158813],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    modeldata["postprocessing"] = FDataDict(
        "file" => "cookstress-ew-princ3",
        "quantity" => :princCauchy,
        "component" => 2,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.4398236425818378, 0.022575693063719465],
    ) < 1.0e-5
    # File = modeldata["postprocessing"]["exported"][1]["file"]
    # @async run(`"paraview.exe" $File`)
    try
        rm(modeldata["postprocessing"]["exported"][1]["file"])
    catch
    end

    AE = AbaqusExporter("Cookstress_algo_stress")
    HEADING(AE, "Cook trapezoidal panel, plane stress")
    COMMENT(AE, "Converged free mid-edge displacement = 23.97")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    COMMENT(AE, "We are assuming three node triangles in plane-stress")
    COMMENT(AE, "CPS3 are pretty poor-accuracy elements, but here we don't care about it.")
    @test nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
    ELEMENT(
        AE,
        "CPS3",
        "AllElements",
        connasarray(modeldata["regions"][1]["femm"].integdomain.fes),
    )
    NSET_NSET(AE, "clamped", modeldata["essential_bcs"][1]["node_list"])
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", thickness)
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.clamped", 2)
    bfes = modeldata["traction_bcs"][1]["femm"].integdomain.fes
    COMMENT(AE, "Concentrated loads: we are assuming that the elements on the boundary")
    COMMENT(AE, "have two nodes each and also that they are the same length.")
    COMMENT(AE, "Then the concentrated loads below will be correctly lumped.")
    nl = connectednodes(bfes)
    F = zeros(count(modeldata["fens"]))
    for ix = 1:count(bfes)
        for jx = 1:2
            F[bfes.conn[ix][jx]] += 1.0 / n / 2 / thickness
        end
    end
    for ixxxx = 1:length(F)
        if F[ixxxx] != 0.0
            CLOAD(AE, "ASSEM1.INSTNC1.$(ixxxx)", 2, F[ixxxx])
        end
    end
    END_STEP(AE)
    close(AE)

    nlines = 0
    open(AE.filename) do f
        s = readlines(f)
        nlines = length(s)
    end
    @test nlines == 2886
    rm(AE.filename)
end
end
using .mmmCookmstressisommm
mmmCookmstressisommm.test()

module mmLE10expimpmm
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()


    # Thick elliptical plate with an elliptical hole is clamped on its exterior
    # boundary and is loaded with transverse  pressure.
    # This is a NAFEMS Benchmark, Test No. LE10.
    # The plate is discretized with solid elements.
    # Because of the symmetries of the geometry and load, only quarter of the plate is modeled.
    # The $\sigma_y=\sigma_2$ at the point $P$ is to be determined. Since the
    # target point is on the boundary of the domain it will not be an
    # integration node as we use Gauss quadrature. The reference value is -5.38 MPa.

    # println("LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole."        )
    t0 = time()

    E = 210e3 * phun("MEGA*PA")# 210e3 MPa
    nu = 0.3
    qmagn = 1.0 * phun("MEGA*PA")# transverse pressure
    sigma_yP = -5.38 * phun("MEGA*PA")# tensile stress at [2.0, 0.0] meters
    Ae = 3.25 * phun("m") # Major radius of the exterior ellipse
    Be = 2.75 * phun("m") # Minor radius of the exterior ellipse
    Ai = 2.0 * phun("m") # Major radius of the interior ellipse
    Bi = 1.0 * phun("m") # Minor radius of the interior ellipse
    Thickness = 0.6 * phun("m")
    tolerance = Thickness / 1000.0 # Geometrical tolerance

    INP_file = """
    *HEADING
    NAFEMS TEST NLE10, COARSE MESH, C3D10M ELEMENTS
    *NODE
           1,     2.83277,       1.348
           2,     2.48116,     1.04837
           3,     2.12955,    0.748738
           4,     3.14146,    0.704779
           5,     2.33893,    0.399071
           6,     2.68977,    0.374369
           7,        3.25,          0.
           8,      2.8335,          0.
           9,       2.417,          0.
          10,     2.83277,       1.348,        0.15
          11,     2.48116,     1.04837,        0.15
          12,     2.12955,    0.748738,        0.15
          13,     2.33841,    0.400381,        0.15
          14,     3.14128,     0.70533,        0.15
          15,        3.25,          0.,        0.15
          16,      2.8335,          0.,        0.15
          17,       2.417,          0.,        0.15
          18,     2.83277,       1.348,         0.3
          19,     2.48116,     1.04837,         0.3
          20,     2.12955,    0.748738,         0.3
          21,     2.62488,       0.674,         0.3
          22,     2.33893,    0.399071,         0.3
          23,     3.14146,    0.704779,         0.3
          24,        3.25,          0.,         0.3
          25,      2.8335,          0.,         0.3
          26,       2.417,          0.,         0.3
          27,     2.83277,       1.348,        0.45
          28,     2.48116,     1.04837,        0.45
          29,     2.12955,    0.748738,        0.45
          30,     2.33841,    0.400381,        0.45
          31,     3.14128,     0.70533,        0.45
          32,        3.25,          0.,        0.45
          33,      2.8335,          0.,        0.45
          34,       2.417,          0.,        0.45
          35,     2.83277,       1.348,         0.6
          36,     2.48116,     1.04837,         0.6
          37,     2.12955,    0.748738,         0.6
          38,     3.14146,    0.704779,         0.6
          39,     2.33893,    0.399071,         0.6
          40,     2.68977,    0.374369,         0.6
          41,        3.25,          0.,         0.6
          42,      2.8335,          0.,         0.6
          43,       2.417,          0.,         0.6
          45,     1.95628,    0.600869
          46,     1.78302,       0.453
          47,     2.06477,    0.374369
          48,     1.93715,    0.248725
          51,      2.2085,          0.
          52,          2.,          0.
          54,     1.95628,    0.600869,        0.15
          55,     1.78302,       0.453,        0.15
          56,     1.93661,    0.249767,        0.15
          59,      2.2085,          0.,        0.15
          60,          2.,          0.,        0.15
          62,     1.95628,    0.600869,         0.3
          63,     1.78302,       0.453,         0.3
          64,     1.93715,    0.248725,         0.3
          66,     2.10001,      0.2265,         0.3
          68,      2.2085,          0.,         0.3
          69,          2.,          0.,         0.3
          71,     1.95628,    0.600869,        0.45
          72,     1.78302,       0.453,        0.45
          74,     1.93661,    0.249767,        0.45
          76,      2.2085,          0.,        0.45
          77,          2.,          0.,        0.45
          79,     1.95628,    0.600869,         0.6
          80,     1.78302,       0.453,         0.6
          81,     2.06477,    0.374369,         0.6
          83,     1.93715,    0.248725,         0.6
          85,      2.2085,          0.,         0.6
          86,          2.,          0.,         0.6
          87,       1.783,     2.29921
          88,     1.57618,     1.80182
          89,     1.36937,     1.30443
          90,     1.95627,     1.52397
          91,     2.36495,     1.88628
          92,     1.78146,     1.06985
          96,       1.783,     2.29921,        0.15
          97,     1.57618,     1.80182,        0.15
          98,     1.36937,     1.30443,        0.15
          99,     1.78071,     1.07038,        0.15
         100,     2.36449,     1.88669,        0.15
         104,       1.783,     2.29921,         0.3
         105,     1.57618,     1.80182,         0.3
         106,     1.36937,     1.30443,         0.3
         107,     2.36495,     1.88628,         0.3
         108,     1.78146,     1.06985,         0.3
         109,     2.10107,     1.32621,         0.3
         113,       1.783,     2.29921,        0.45
         114,     1.57618,     1.80182,        0.45
         115,     1.36937,     1.30443,        0.45
         116,     2.36449,     1.88669,        0.45
         117,     1.78071,     1.07038,        0.45
         121,       1.783,     2.29921,         0.6
         122,     1.57618,     1.80182,         0.6
         123,     1.36937,     1.30443,         0.6
         124,     1.95627,     1.52397,         0.6
         125,     2.36495,     1.88628,         0.6
         126,     1.78146,     1.06985,         0.6
         131,     1.26718,     1.05863
         132,       1.165,    0.812831
         134,     1.49321,    0.665266
         135,     1.64727,    0.780784
         140,     1.26718,     1.05863,        0.15
         141,       1.165,    0.812831,        0.15
         143,      1.4924,    0.665723,        0.15
         148,     1.26718,     1.05863,         0.3
         149,       1.165,    0.812831,         0.3
         150,     1.57619,    0.878714,         0.3
         152,     1.49321,    0.665266,         0.3
         157,     1.26718,     1.05863,        0.45
         158,       1.165,    0.812831,        0.45
         160,      1.4924,    0.665723,        0.45
         165,     1.26718,     1.05863,         0.6
         166,       1.165,    0.812831,         0.6
         168,     1.49321,    0.665266,         0.6
         169,     1.64727,    0.780784,         0.6
         173,          0.,       1.583
         174,          0.,      1.2915
         175,          0.,          1.
         176,      0.5825,     1.19792
         177,    0.699442,     1.51527
         178,    0.590531,    0.955415
         182,          0.,       1.583,        0.15
         183,          0.,      1.2915,        0.15
         184,          0.,          1.,        0.15
         185,    0.698861,     1.51538,        0.15
         186,    0.590076,    0.955486,        0.15
         190,          0.,       1.583,         0.3
         191,          0.,      1.2915,         0.3
         192,          0.,          1.,         0.3
         193,    0.699442,     1.51527,         0.3
         194,    0.590531,    0.955415,         0.3
         195,    0.684684,     1.15221,         0.3
         199,          0.,       1.583,        0.45
         200,          0.,      1.2915,        0.45
         201,          0.,          1.,        0.45
         202,    0.698861,     1.51538,        0.45
         203,    0.590076,    0.955486,        0.45
         207,          0.,       1.583,         0.6
         208,          0.,      1.2915,         0.6
         209,          0.,          1.,         0.6
         210,      0.5825,     1.19792,         0.6
         211,    0.699442,     1.51527,         0.6
         212,    0.590531,    0.955415,         0.6
         216,          0.,        2.75
         217,          0.,      2.1665
         219,    0.920251,     2.63745
         221,      0.8915,      1.9411
         225,          0.,        2.75,        0.15
         226,          0.,      2.1665,        0.15
         228,    0.919707,     2.63759,        0.15
         233,          0.,        2.75,         0.3
         234,          0.,      2.1665,         0.3
         236,    0.684684,     2.02721,         0.3
         237,    0.920251,     2.63745,         0.3
         242,          0.,        2.75,        0.45
         243,          0.,      2.1665,        0.45
         245,    0.919707,     2.63759,        0.45
         250,          0.,        2.75,         0.6
         251,          0.,      2.1665,         0.6
         253,    0.920251,     2.63745,         0.6
         255,      0.8915,      1.9411,         0.6
    **
    **
    *ELEMENT, TYPE=C3D10M, ELSET=EALL
           1,       1,       7,      18,       3,       4,      14,      10,
           2,       6,      11
           2,      24,      18,       7,      26,      23,      14,      15,
          25,      21,      16
           3,       9,       3,      26,       7,       5,      13,      17,
           8,       6,      16
           4,      20,      26,       3,      18,      22,      13,      12,
          19,      21,      11
           5,       3,       7,      18,      26,       6,      14,      11,
          13,      16,      21
           6,      24,      41,      18,      26,      32,      31,      23,
          25,      33,      21
           7,      35,      18,      41,      37,      27,      31,      38,
          36,      28,      40
           8,      20,      37,      26,      18,      29,      30,      22,
          19,      28,      21
           9,      43,      26,      37,      41,      34,      30,      39,
          42,      33,      40
          10,      18,      26,      41,      37,      21,      33,      31,
          28,      30,      40
          11,       9,      26,       3,      52,      17,      13,       5,
          51,      59,      47
          12,      20,       3,      26,      63,      12,      13,      22,
          62,      54,      66
          13,      46,      63,      52,       3,      55,      56,      48,
          45,      54,      47
          14,      69,      52,      63,      26,      60,      56,      64,
          68,      59,      66
          15,       3,      52,      26,      63,      47,      59,      13,
          54,      56,      66
          16,      20,      26,      37,      63,      22,      30,      29,
          62,      66,      71
          17,      43,      37,      26,      86,      39,      30,      34,
          85,      81,      76
          18,      69,      63,      86,      26,      64,      74,      77,
          68,      66,      76
          19,      80,      86,      63,      37,      83,      74,      72,
          79,      81,      71
          20,      63,      26,      37,      86,      66,      30,      71,
          74,      76,      81
          21,       1,      18,      87,       3,      10,     100,      91,
           2,      11,      90
          22,     104,      87,      18,     106,      96,     100,     107,
         105,      97,     109
          23,      89,     106,       3,      87,      98,      99,      92,
          88,      97,      90
          24,      20,       3,     106,      18,      12,      99,     108,
          19,      11,     109
          25,      87,       3,      18,     106,      90,      11,     100,
          97,      99,     109
          26,     104,      18,     121,     106,     107,     116,     113,
         105,     109,     114
          27,      35,     121,      18,      37,     125,     116,      27,
          36,     124,      28
          28,      20,     106,      37,      18,     108,     117,      29,
          19,     109,      28
          29,     123,      37,     106,     121,     126,     117,     115,
         122,     124,     114
          30,     106,      18,     121,      37,     109,     116,     114,
         117,      28,     124
          31,      89,       3,     106,     132,      92,      99,      98,
         131,     135,     140
          32,      20,     106,       3,      63,     108,      99,      12,
          62,     150,      54
          33,      46,     132,      63,       3,     134,     143,      55,
          45,     135,      54
          34,     149,      63,     132,     106,     152,     143,     141,
         148,     150,     140
          35,     132,       3,     106,      63,     135,      99,     140,
         143,      54,     150
          36,      20,      37,     106,      63,      29,     117,     108,
          62,      71,     150
          37,     123,     106,      37,     166,     115,     117,     126,
         165,     157,     169
          38,     149,     166,      63,     106,     158,     160,     152,
         148,     157,     150
          39,      80,      63,     166,      37,      72,     160,     168,
          79,      71,     169
          40,     106,      63,      37,     166,     150,      71,     117,
         157,     160,     169
          41,      89,     106,     173,     132,      98,     185,     177,
         131,     140,     176
          42,     190,     173,     106,     192,     182,     185,     193,
         191,     183,     195
          43,     175,     192,     132,     173,     184,     186,     178,
         174,     183,     176
          44,     149,     132,     192,     106,     141,     186,     194,
         148,     140,     195
          45,     173,     132,     106,     192,     176,     140,     185,
         183,     186,     195
          46,     190,     106,     207,     192,     193,     202,     199,
         191,     195,     200
          47,     123,     207,     106,     166,     211,     202,     115,
         165,     210,     157
          48,     149,     192,     166,     106,     194,     203,     158,
         148,     195,     157
          49,     209,     166,     192,     207,     212,     203,     201,
         208,     210,     200
          50,     192,     106,     207,     166,     195,     202,     200,
         203,     157,     210
          51,     216,      87,     233,     173,     219,     228,     225,
         217,     221,     226
          52,     104,     233,      87,     106,     237,     228,      96,
         105,     236,      97
          53,      89,     173,     106,      87,     177,     185,      98,
          88,     221,      97
          54,     190,     106,     173,     233,     193,     185,     182,
         234,     236,     226
          55,     173,      87,     233,     106,     221,     228,     226,
         185,      97,     236
          56,     104,     121,     233,     106,     113,     245,     237,
         105,     114,     236
          57,     250,     233,     121,     207,     242,     245,     253,
         251,     243,     255
          58,     190,     207,     106,     233,     199,     202,     193,
         234,     243,     236
          59,     123,     106,     207,     121,     115,     202,     211,
         122,     114,     255
          60,     233,     106,     121,     207,     236,     114,     245,
         243,     202,     255
    """
    write("NLE10.inp", INP_file)

    output = MeshImportModule.import_ABAQUS("NLE10.inp")
    fens, fes = output["fens"], output["fesets"][1]

    try
        rm("NLE10.inp")
    catch
    end

    # Select the  boundary faces, on the boundary that is clamped,  and on the part
    # of the boundary that is loaded with the transverse pressure
    bdryfes = meshboundary(fes)
    exteriorbfl = selectelem(
        fens,
        bdryfes,
        facing = true,
        direction = [1.0, 1.0, 0.0],
        dotmin = 0.001,
    )
    topbfl = selectelem(
        fens,
        bdryfes,
        box = [0.0, Inf, 0.0, Inf, Thickness, Thickness],
        inflate = tolerance,
    )

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    L12 = connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
    setebc!(u, L12, true, 1, 0.0)
    setebc!(u, L12, true, 2, 0.0)
    LL = selectnode(
        fens;
        box = [0.0, Inf, 0.0, Inf, Thickness / 2.0, Thickness / 2.0],
        inflate = tolerance,
    )
    L3 = intersect(LL, connectednodes(subset(bdryfes, exteriorbfl)))
    setebc!(u, L3, true, 3, 0.0)
    L1 = selectnode(fens; box = [0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u, L1, true, 1, 0.0) # symmetry plane X = 0
    L2 = selectnode(fens; box = [0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u, L2, true, 2, 0.0) # symmetry plane Y = 0

    applyebc!(u)
    numberdofs!(u)

    eL1femm = FEMMBase(IntegDomain(subset(bdryfes, topbfl), TriRule(3)))
    function pfun(forceout::Vector{T}, XYZ, tangents, feid, qpid) where {T}
        forceout .= [0.0, 0.0, -qmagn]
        return forceout
    end
    fi = ForceIntensity(Float64, 3, pfun)
    F2 = distribloads(eL1femm, geom, u, fi, 2)

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)

    u = solve_blocked!(u, K, F2)

    nl = selectnode(fens, box = [Ai, Ai, 0, 0, Thickness, Thickness], inflate = tolerance)
    thecorneru = zeros(Float64, 1, 3)
    gathervalues_asmat!(u, thecorneru, nl)
    thecorneru = thecorneru / phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")
    @test norm(thecorneru - [-0.0268854 0.0 -0.0919955]) < 1.0e-5

    fld = fieldfromintegpoints(
        femm,
        geom,
        u,
        :Cauchy,
        2;
        nodevalmethod = :invdistance,
        reportat = :meanonly,
    )#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl, 1][1] / phun("MPa") - -2.2616980060965024) < 1.0e-3

    fld = fieldfromintegpoints(
        femm,
        geom,
        u,
        :Cauchy,
        2;
        nodevalmethod = :averaging,
        reportat = :extrapmean,
    )#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl, 1][1] / phun("MPa") - -2.382478776709117) < 1.0e-3

    fld = fieldfromintegpoints(
        femm,
        geom,
        u,
        :Cauchy,
        2;
        nodevalmethod = :averaging,
        reportat = :extraptrend,
    )#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    # println("$(fld.values[nl,1][1]/phun("MPa"))")
    @test abs(fld.values[nl, 1][1] / phun("MPa") - -5.470291697493607) < 1.0e-3

    File = "LE10NAFEMS_MST10_sigmay.vtk"
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.T10;
        vectors = [("u", u.values)],
        scalars = [("sigmay", fld.values)],
    )
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    AE = AbaqusExporter("LE10NAFEMS_MST10")
    HEADING(
        AE,
        "LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole.",
    )
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "c3d10", "AllElements", 1, connasarray(femm.integdomain.fes))
    ELEMENT(
        AE,
        "SFM3D6",
        "TractionElements",
        1 + count(femm.integdomain.fes),
        connasarray(eL1femm.integdomain.fes),
    )
    NSET_NSET(AE, "L1", L1)
    NSET_NSET(AE, "L2", L2)
    NSET_NSET(AE, "L3", L3)
    NSET_NSET(AE, "L12", L12)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl")
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L1", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L2", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L3", 3)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L12", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.L12", 2)
    DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec([0.0, 0.0, -qmagn]))
    END_STEP(AE)
    close(AE)

    output = MeshImportModule.import_ABAQUS(AE.filename)
    fens, fes = output["fens"], output["fesets"][1]
    try
        rm(AE.filename)
    catch
    end

    # Select the  boundary faces, on the boundary that is clamped,  and on the part
    # of the boundary that is loaded with the transverse pressure
    bdryfes = meshboundary(fes)
    exteriorbfl = selectelem(
        fens,
        bdryfes,
        facing = true,
        direction = [1.0, 1.0, 0.0],
        dotmin = 0.001,
    )
    topbfl = selectelem(
        fens,
        bdryfes,
        box = [0.0, Inf, 0.0, Inf, Thickness, Thickness],
        inflate = tolerance,
    )

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    L12 = connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
    setebc!(u, L12, true, 1, 0.0)
    setebc!(u, L12, true, 2, 0.0)
    LL = selectnode(
        fens;
        box = [0.0, Inf, 0.0, Inf, Thickness / 2.0, Thickness / 2.0],
        inflate = tolerance,
    )
    L3 = intersect(LL, connectednodes(subset(bdryfes, exteriorbfl)))
    setebc!(u, L3, true, 3, 0.0)
    L1 = selectnode(fens; box = [0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u, L1, true, 1, 0.0) # symmetry plane X = 0
    L2 = selectnode(fens; box = [0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u, L2, true, 2, 0.0) # symmetry plane Y = 0

    applyebc!(u)
    numberdofs!(u)

    eL1femm = FEMMBase(IntegDomain(subset(bdryfes, topbfl), TriRule(3)))
    # function pfun(forceout::FVec{T}, XYZ, tangents, feid) where {T}
    #     forceout .=  [0.0, 0.0, -qmagn]
    #     return forceout
    # end
    fi = ForceIntensity(Float64, 3, pfun)
    F2 = distribloads(eL1femm, geom, u, fi, 2)

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)

    u = solve_blocked!(u, K, F2)

    nl = selectnode(fens, box = [Ai, Ai, 0, 0, Thickness, Thickness], inflate = tolerance)
    thecorneru = zeros(Float64, 1, 3)
    gathervalues_asmat!(u, thecorneru, nl)
    thecorneru = thecorneru / phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")
    @test norm(thecorneru - [-0.0268854 0.0 -0.0919955]) < 1.0e-5

    fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 2; reportat = :meanonly)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl, 1][1] / phun("MPa") - -2.2616980060965024) < 1.0e-3

    fld = fieldfromintegpoints(
        femm,
        geom,
        u,
        :Cauchy,
        2;
        nodevalmethod = :averaging,
        reportat = :extrapmean,
    )#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    @test abs(fld.values[nl, 1][1] / phun("MPa") - -2.382478776709117) < 1.0e-3

    fld = fieldfromintegpoints(
        femm,
        geom,
        u,
        :Cauchy,
        2;
        nodevalmethod = :averaging,
        reportat = :extraptrend,
    )#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")
    # println("$(fld.values[nl,1][1]/phun("MPa"))")
    @test abs(fld.values[nl, 1][1] / phun("MPa") - -5.470291697493607) < 1.0e-3

    File = "LE10NAFEMS_MST10_sigmay.vtk"
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.T10;
        vectors = [("u", u.values)],
        scalars = [("sigmay", fld.values)],
    )
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end
end
end
using .mmLE10expimpmm
mmLE10expimpmm.test()


# module mmtruncatedmfreem1
# using FinEtools
# using FinEtoolsDeforLinear
# using FinEtoolsDeforLinear.AlgoDeforLinearModule
# using Test
# import Arpack: eigs
# import LinearAlgebra: norm, cholesky, cross
# function test()
#     # println("""
#     # Vibration modes of truncated cylindrical shell.
#     # """)

#     # t0 = time()

#     E = 205000*phun("MPa");# Young's modulus
#     nu = 0.3;# Poisson ratio
#     rho = 7850*phun("KG*M^-3");# mass density
#     OmegaShift = (2*pi*100) ^ 2; # to resolve rigid body modes
#     h = 0.05*phun("M");
#     l = 10*h;
#     Rmed = h/0.2;
#     psi   = 0;    # Cylinder
#     nh = 5; nl  = 12; nc = 40;
#     tolerance = h/nh/100;
#     neigvs = 20;

#     MR = DeforModelRed3D
#     fens,fes  = H8block(h,l,2*pi,nh,nl,nc)
#     # Shape into a cylinder
#     R = zeros(3, 3)
#     for i = 1:count(fens)
#         x, y, z = fens.xyz[i,:];
#         rotmat3!(R, [0, z, 0])
#         Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
#             -sin(psi*pi/180) cos(psi*pi/180) 0;
#             0 0 1]
#         fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
#     end
#     candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h/1000)
#     fens,fes = mergenodes(fens, fes,  tolerance, candidates);

#     geom = NodalField(fens.xyz)
#     u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

#     numberdofs!(u)

#     material=MatDeforElastIso(MR, rho, E, nu, 0.0)

#     femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)), material)
#     femm = associategeometry!(femm, geom)
#     K =stiffness(femm, geom, u)
#     femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
#     M =mass(femm, geom, u)


#     # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
#     # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
#     # of iterations niter and the number of matrix vector multiplications nmult, as
#     # well as the final residual vector resid.

#     if true
#         d,v,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
#         d[:] = d .- OmegaShift;
#         fs = real(sqrt.(complex(d)))/(2*pi)
#         # println("Eigenvalues: $fs [Hz]")
#         @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 0.000166835, 0.000182134, 517.147, 517.147]) < 2.0e-2
#         # mode = 7
#         # scattersysvec!(u, v[:,mode])
#         # File =  "unit_cube_modes.vtk"
#         # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
#         # @async run(`"paraview.exe" $File`)
#     end


#     v0 = [i==j ? one(Float64) : zero(Float64) for i=1:size(K,1), j=1:2*neigvs]
#     tol = 1.0e-2
#     maxiter = 20
#     lamb, v, nconv, niter, lamberr =
#     AlgoDeforLinearModule.ssit(K+OmegaShift*M, M; nev=neigvs, v0=v0, tol=tol, maxiter=maxiter)
#     @test nconv == neigvs
#         # if nconv < neigvs
#             # println("NOT converged")
#         # end
#     broadcast!(+, lamb, lamb, - OmegaShift);
#     fs = real(sqrt.(complex(lamb)))/(2*pi)
#         # println("Eigenvalues: $fs [Hz]")
#         # println("$(sort(fs))")
#     @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 7.9048e-5, 0.0, 517.147, 517.147]) < 2.0e-2
#         # println("Eigenvalue errors: $lamberr [ND]")
#         # mode = 7
#         # scattersysvec!(u, v[:,mode])
#         # File =  "unit_cube_modes.vtk"
#         # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
#         # @async run(`"paraview.exe" $File`)


# end
# end
# using .mmtruncatedmfreem1
# mmtruncatedmfreem1.test()

module mmFV32mm1
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    # println("""
    # FV32: Cantilevered tapered membrane
    # This is a test recommended by the National Agency for Finite Element Methods and
    # Standards (U.K.): Test FV32 from NAFEMS publication TNSB, Rev. 3, “The Standard
    # NAFEMS Benchmarks,” October 1990.
    #
    # Reference solution: 44.623	130.03	162.70	246.05	379.90	391.44 for the first
    # six modes.
    # """)

    t0 = time()


    E = 200 * phun("GPA")
    nu = 0.3
    rho = 8000 * phun("KG/M^3")
    L = 10 * phun("M")
    W0 = 5 * phun("M")
    WL = 1 * phun("M")
    H = 0.05 * phun("M")
    nL, nW, nH = 8, 4, 1# How many element edges per side?
    neigvs = 20                   # how many eigenvalues
    Reffs = [44.623 130.03 162.70 246.05 379.90 391.44]

    fens, fes = H20block(1.0, 2.0, 1.0, nL, nW, nH)
    for i = 1:count(fens)
        xi, eta, theta = fens.xyz[i, :]
        eta = eta - 1.0
        fens.xyz[i, :] = [xi * L eta * (1.0 - 0.8 * xi) * W0 / 2 theta * H / 2]
    end
    # File =  "mesh.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, rho, E, nu, 0.0)
    region1 = FDataDict(
        "femm" => FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material),
        "femm_mass" => FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 3)), material),
    )

    nl1 = selectnode(fens; plane = [1.0 0.0 0.0 0.0], thickness = H / 1.0e4)
    ebc1 = FDataDict("node_list" => nl1, "component" => 1, "displacement" => 0.0)
    ebc2 = FDataDict("node_list" => nl1, "component" => 2, "displacement" => 0.0)
    ebc3 = FDataDict("node_list" => nl1, "component" => 3, "displacement" => 0.0)

    nl4 = selectnode(fens; plane = [0.0 0.0 1.0 0.0], thickness = H / 1.0e4)
    ebc4 = FDataDict("node_list" => nl4, "component" => 3, "displacement" => 0.0)

    # Make model data
    modeldata = FDataDict(
        "fens" => fens,
        "regions" => [region1],
        "essential_bcs" => [ebc1 ebc2 ebc3 ebc4],
        "neigvs" => neigvs,
    )

    # Solve
    modeldata = FinEtoolsDeforLinear.AlgoDeforLinearModule.modal(modeldata)

    fs = modeldata["omega"] / (2 * pi)
    #println("Eigenvalues: $fs [Hz]")
    # println("Percentage frequency errors: $((vec(fs[1:6]) - vec(Reffs))./vec(Reffs)*100)")
    @test norm(
        (vec(fs[1:6]) - vec(Reffs)) ./ vec(Reffs) * 100 -
        [0.0162775, 0.0623384, 0.00799148, 0.151669, 0.376663, 0.0191388],
    ) < 1.0e-6

    modeldata["postprocessing"] = FDataDict("file" => "FV32-modes", "mode" => 1:10)
    modeldata = FinEtoolsDeforLinear.AlgoDeforLinearModule.exportmode(modeldata)
    # @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
    try
        rm(modeldata["postprocessing"]["file"] * "1.vtk")
    catch
    end
end
end
using .mmFV32mm1
mmFV32mm1.test()

# module mmtruncatedmfreem2
# using FinEtools
# using FinEtoolsDeforLinear
# using FinEtoolsDeforLinear.AlgoDeforLinearModule
# using Test
# import Arpack: eigs
# import LinearAlgebra: norm, cholesky, cross
# function test()
#     # println("""
#     # Vibration modes of truncated cylindrical shell.
#     # """)

#     # t0 = time()

#     E = 205000*phun("MPa");# Young's modulus
#     nu = 0.3;# Poisson ratio
#     rho = 7850*phun("KG*M^-3");# mass density
#     OmegaShift = (2*pi*100) ^ 2; # to resolve rigid body modes
#     h = 0.05*phun("M");
#     l = 10*h;
#     Rmed = h/0.2;
#     psi   = 0;    # Cylinder
#     nh = 5; nl  = 12; nc = 40;
#     tolerance = h/nh/100;
#     neigvs = 20;

#     MR = DeforModelRed3D
#     fens,fes  = H8block(h,l,2*pi,nh,nl,nc)
#     # Shape into a cylinder
#     R = zeros(3, 3)
#     for i = 1:count(fens)
#         x, y, z = fens.xyz[i,:];
#         rotmat3!(R, [0, z, 0])
#         Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
#             -sin(psi*pi/180) cos(psi*pi/180) 0;
#             0 0 1]
#         fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
#     end
#     candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h/1000)
#     fens,fes = mergenodes(fens, fes,  tolerance, candidates);

#     geom = NodalField(fens.xyz)
#     u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

#     numberdofs!(u)

#     material=MatDeforElastIso(MR, rho, E, nu, 0.0)

#     femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)), material)
#     femm = associategeometry!(femm, geom)
#     K =stiffness(femm, geom, u)
#     femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
#     M =mass(femm, geom, u)


#     # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
#     # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
#     # of iterations niter and the number of matrix vector multiplications nmult, as
#     # well as the final residual vector resid.

#     if true
#         d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
#         broadcast!(+, d, d, - OmegaShift);
#         fs = real(sqrt.(complex(d)))/(2*pi)
#         # println("Eigenvalues: $fs [Hz]")
#         @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 0.000166835, 0.000182134, 517.147, 517.147]) < 2.0e-2
#         # mode = 7
#         # scattersysvec!(u, v[:,mode])
#         # File =  "unit_cube_modes.vtk"
#         # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
#         # @async run(`"paraview.exe" $File`)
#     end

#     v0 = [i==j ? one(Float64) : zero(Float64) for i=1:size(K,1), j=1:2*neigvs]
#     tol = 1.0e-2
#     maxiter = 20
#     lamb, v, nconv, niter, lamberr =
#     AlgoDeforLinearModule.ssit(K+OmegaShift*M, M; nev=neigvs, v0=v0, tol=tol, maxiter=maxiter)
#     @test nconv == neigvs
#         # if nconv < neigvs
#             # println("NOT converged")
#         # end
#     lamb = lamb .- OmegaShift;
#     fs = real(sqrt.(complex(lamb)))/(2*pi)
#         # println("Eigenvalues: $fs [Hz]")
#         # println("$(sort(fs))")
#     @test norm(sort(fs)[1:8] - [0.0, 0.0, 0.0, 0.0, 7.9048e-5, 0.0, 517.147, 517.147]) < 2.0e-2
#         # println("Eigenvalue errors: $lamberr [ND]")
#         # mode = 7
#         # scattersysvec!(u, v[:,mode])
#         # File =  "unit_cube_modes.vtk"
#         # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
#         # @async run(`"paraview.exe" $File`)

# end
# end
# using .mmtruncatedmfreem2
# mmtruncatedmfreem2.test()
