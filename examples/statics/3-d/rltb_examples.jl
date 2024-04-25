module rltb_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy, solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule:
    linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky

# Isotropic material
E = 1000.0
nu = 0.4999 # Taylor data
nu = 0.3 # Taylor data
W = 2.5
H = 5.0
L = 50.0
htol = minimum([L, H, W]) / 1000
uzex = -12.6
magn = 0.2 * uzex / 4
Force = magn * W * H * 2
CTE = 0.0
mult = 8
n = 3 #

function getfrcL!(
    forceout,
    XYZ,
    tangents,
    feid,
    qpid,
)
    copyto!(forceout, [0.0; 0.0; magn])
    return forceout
end

function rltb_H8_by_hand()
    elementtag = "H8"
    println("""
    Taylor Cantilever example. Element: $(elementtag)
    """)

    fens, fes = H8block(W, L, H, n, mult * n, 2 * n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(
        csmatout,
        XYZ,
        tangents,
        feid,
        qpid,
    )
        copyto!(csmatout, csmat)
        return csmatout
    end

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2)
    
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    u = solve_blocked!(u, K, F2)

    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "rltb_H8_by_hand.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # rltb_H8_by_hand

function rltb_H8_sri_by_hand()
    elementtag = "H8-SRI"
    println("""
    Taylor Cantilever example. Element: $(elementtag)
    """)

    fens, fes = H8block(W, L, H, n, mult * n, 2 * n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(
        csmatout,
        XYZ,
        tangents,
        feid,
        qpid,
    )
        copyto!(csmatout, csmat)
        return csmatout
    end

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2)
    

    # First compute the bulk portion of the stiffness matrix
    material.D .= FinEtoolsDeforLinear.MatDeforElastIsoModule.bulk_split_tangent_moduli_bulk(E, nu)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    # Next add the shear portion of the stiffness matrix
    material.D .= FinEtoolsDeforLinear.MatDeforElastIsoModule.bulk_split_tangent_moduli_shear(E, nu)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    associategeometry!(femm, geom)
    K += stiffness(femm, geom, u)

    u = solve_blocked!(u, K, F2)

    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "rltb_H8_by_hand.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # rltb_H8_sri_by_hand

function rltb_H20_by_hand()
    elementtag = "H20"
    println("""
    Taylor Cantilever example. Element: $(elementtag)
    """)

    fens, fes = H20block(W, L, H, n, mult * n, 2 * n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid)
        copyto!(csmatout, csmat)
    end

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2)
    
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    u = solve_blocked!(u, K, F2)

    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "rltb_H20_by_hand.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    # @async run(`"paraview.exe" $File`)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)", "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy, "component"=>[5])
    # modeldata = exportstresselementwise(modeldata)

    # modeldata["postprocessing"] = FDataDict("file"=>"hughes_cantilever_stresses_$(elementtag)",
    # "outputcsys"=>CSys(3, 3, updatecs!), "quantity"=>:Cauchy,
    # "component"=>collect(1:6))
    # modeldata = exportstresselementwise(modeldata)
    # stressfields = ElementalField[modeldata["postprocessing"]["exported"][1]["field"]]

    true
end # rltb_H20_by_hand

function rltb_H8_abaqus_export()
    elementtag = "H8"
    println("""
    Taylor Cantilever example. Element: $(elementtag)
    """)

    fens, fes = H8block(W, L, H, n, mult * n, 2 * n)
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(
        csmatout,
        XYZ,
        tangents,
        feid,
        qpid,
    )
        copyto!(csmatout, csmat)
        return csmatout
    end

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F2 = distribloads(el2femm, geom, u, fi, 2)
    
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    u = solve_blocked!(u, K, F2)

    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    AE = AbaqusExporter("rltb_H8_abaqus_export")
    HEADING(
        AE,
        "Taylor Cantilever example. Element: $(elementtag)",
    )
    PART(AE, "PART1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "C3D8", "ALLELEMENTS", 1, finite_elements(femm).conn)
    ELEMENT(
        AE,
        "SFM3D4",
        "TRACTIONELEMENTS",
        1 + count(finite_elements(femm)),
        finite_elements(el2femm).conn,
    )
    NSET_NSET(AE, "LX0", lx0)
    NSET_NSET(AE, "LX1", lx1)
    ORIENTATION(AE, "GLOBALORIENTATION", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "ELASTICITY", "GLOBALORIENTATION", "ALLELEMENTS", "HOURGLASSCTL")
    SURFACE_SECTION(AE, "TRACTIONELEMENTS")
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "ELASTICITY")
    ELASTIC(AE, E, nu)
    SECTION_CONTROLS(AE, "HOURGLASSCTL", "HOURGLASS=ENHANCED")
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.LX0", 1)
    BOUNDARY(AE, "ASSEM1.INSTNC1.LX0", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.LX0", 3)
    BOUNDARY(AE, "ASSEM1.INSTNC1.LX1", 1)
    DLOAD(AE, "ASSEM1.INSTNC1.TRACTIONELEMENTS", vec([0.0, 0.0, -magn]))
    END_STEP(AE)
    close(AE)

    true
end # rltb_H8_by_hand

function allrun()
    println("#####################################################")
    println("# rltb_H8_by_hand ")
    rltb_H8_by_hand()
    println("#####################################################")
    println("# rltb_H8_sri_by_hand ")
    rltb_H8_sri_by_hand()
    println("#####################################################")
    println("# rltb_H20_by_hand ")
    rltb_H20_by_hand()
    println("#####################################################")
    println("# rltb_H8_abaqus_export ")
    rltb_H8_abaqus_export()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
