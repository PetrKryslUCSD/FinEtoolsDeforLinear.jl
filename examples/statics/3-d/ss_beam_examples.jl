module ss_beam_examples

using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy, solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule: linearstatics,
    exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky

# Isotropic material
E = 209000.0
nu = 0.3 #Taylor data
W = 30.0
H = 16.25
L = 200.0
htol = minimum([L, H, W]) / 1000
magn = 1.0
uzex = -5 / 384 * magn * W * L^4 / (E * W * H^3 / 12)
n = 8 #

function getfrcL!(forceout::FFltVec,
    XYZ::FFltMat,
    tangents::FFltMat,
    feid::FInt,
    qpid::FInt)
    copyto!(forceout, [0.0; 0.0; -magn])
    forceout
end

function test_h8()
    elementtag = "H8"
    println("""
    SS beam example. Element: $(elementtag)
    """)

    fens, fes = H8block(L, W, H, n, 2, 1)
    bfes = meshboundary(fes)
    # Simple support
    sectionL = selectelem(fens, bfes; facing = true, direction = [+1.0 0.0 0.0])
    # Simple support
    section0 = selectelem(fens, bfes; facing = true, direction = [-1.0 0.0 0.0])
    # Top loaded surface
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [0.0 0.0 1.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

    # Material orientation matrix
    csmat = [i == j ? one(FFlt) : zero(FFlt) for i in 1:3, j in 1:3]

    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, feid::FInt)
        copyto!(csmatout, csmat)
    end

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx0 = connectednodes(subset(bfes, sectionL))
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    # Fix one node in the X direction
    lx1 = connectednodes(subset(bfes, sectionlateral))[1]
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionlateral), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    u = solve_blocked!(u, K, F2)

    utip = minimum(u.values[:, 3])
    println("Deflection: $(utip), compared to $(uzex)")

    File = "test_h8.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # test_h8

function test_h20r()
    elementtag = "H20R"
    println("""
    SS beam example. Element: $(elementtag)
    """)

    fens, fes = H20block(L, W, H, n, 2, 1)
    bfes = meshboundary(fes)
    # Simple support
    sectionL = selectelem(fens, bfes; facing = true, direction = [+1.0 0.0 0.0])
    # Simple support
    section0 = selectelem(fens, bfes; facing = true, direction = [-1.0 0.0 0.0])
    # Top loaded surface
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [0.0 0.0 1.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)

    # Material orientation matrix
    csmat = [i == j ? one(FFlt) : zero(FFlt) for i in 1:3, j in 1:3]

    function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, feid::FInt)
        copyto!(csmatout, csmat)
    end

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx0 = connectednodes(subset(bfes, sectionL))
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    # Fix one node in the X direction
    lx1 = connectednodes(subset(bfes, sectionlateral))[1]
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionlateral), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    u = solve_blocked!(u, K, F2)

    utip = minimum(u.values[:, 3])
    println("Deflection: $(utip), compared to $(uzex)")

    File = "test_h20r.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # test_h8

function allrun()
    println("#####################################################")
    println("# test_h8 ")
    test_h8()
    println("#####################################################")
    println("# test_h20r ")
    test_h20r()

    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
