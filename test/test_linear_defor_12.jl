"""
Composite tests
R0031(1): Laminated strip under three-point bending
R0031(2): Wrapped thick cylinder under pressure and thermal loading
R0031(3): Three-layer sandwich shell under normal pressure loading
"""

module NAFEMS_R0031_1_test_msh8
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Statistics: mean
using LinearAlgebra: norm, cross
using Test

function NAFEMS_R0031_1_msh8()
    # println("""
    # Laminated Strip Under Three-Point Bending. Mean-strain h8 elements.
    # """)

    # Determine the central transverse displacement in a simply-supported seven
    # layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
    # material lay-up is specified with the center ply being four times as
    # thick as the others.
    # Reference: NAFEMS Report R0031, Test No.1, 17-Dec-1998.

    # Because of the symmetries of the geometry and load, only the
    # first-quadrant   (in XY) quarter of the plate is modeled.

    # The coordinate system is centered at point E (at the difference with
    # respect to the original benchmark definition).  The  load is applied
    # along a curve passing through point C. The simple support is applied
    # along the curve passing through point B.

    # We realize the simple supports along the lines  A, B and the line load at
    # point C  are illegal from the point of view of convergence.  No
    # convergence can be hoped for as the stress underneath the load and above
    # the simple supports  is infinite in the limit (these locations are stress
    # singularities).   However, for relatively coarse meshes the results away
    # from the singularities are still meaningful.

    # The target quantities are displacement at the bottom surface at point E,
    # the tensile axial stress at the same point,  and of the transverse shear
    # stress at point D  in between the bottommost two layers (See figure 1).

    t0 = time()
    # Orthotropic material parameters of the material of the layers
    E1s = 100.0 * phun("GPa")
    E2s = E3s = 5.0 * phun("GPa")
    nu12s = 0.4
    nu13s = 0.3
    nu23s = 0.3
    G12s = 3.0 * phun("GPa")
    G13s = G23s = 2.0 * phun("GPa")
    CTE1 = 3.0e-6
    CTE2 = 2.0e-5
    CTE3 = 2.0e-5

    AB = 30.0 * phun("mm") # span between simple supports
    CD = 4.0 * phun("mm") # distance between the point D and the center
    OH = 10.0 * phun("mm") # overhang
    W = 10.0 * phun("mm") # width of the strip

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0 90 0 90 0 90 0])
    ts = vec([0.1 0.1 0.1 0.4 0.1 0.1 0.1]) * phun("mm") # layer thicknesses
    TH = sum(ts) # total thickness of the plate

    tolerance = 0.0001 * TH

    # The line load is in the negative Z direction.
    q0 = 10 * phun("N/mm") #    line load

    # Reference deflection under the load is
    wEref = -1.06 * phun("mm")

    # The reference tensile stress at the bottom of the lowest layer is
    sigma11Eref = 684 * phun("MPa")

    # Because we model the first-quadrant quarter of the plate using
    # coordinate axes centered  at the point E  the shear at the point D is
    # positive instead of negative as in the benchmark where the coordinate
    # system is located at the outer corner of the strip.
    sigma13Dref = 4.1 * phun("MPa")

    Refinement = 4
    # We select 8 elements spanwise and 2 elements widthwise.  The overhang
    # of the plate is given one element.
    nL = Refinement * 4
    nO = Refinement * 1
    nW = Refinement * 1

    # Each layer is modeled with a single element times the refinement.
    nts = Refinement * ones(Int, length(angles))# number of elements per layer

    xs = unique(vcat(collect(linearspace(0.0, AB / 2, nL + 1)), [CD],
        collect(linearspace(AB / 2, AB / 2 + OH, nO + 1))))
    xs = xs[sortperm(xs)]
    ys = collect(linearspace(0.0, W / 2, nW + 1))

    fens, fes = H8layeredplatex(xs, ys, ts, nts)

    # This is the material  model
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
        0.0, E1s, E2s, E3s,
        nu12s, nu13s, nu23s,
        G12s, G13s, G23s,
        CTE1, CTE2, CTE3)

    # The material coordinate system function is defined as:
    function _updatecs!(csmatout, feid, labels)
        rotmat3!(csmatout, angles[labels[feid]] / 180.0 * pi * [0.0; 0.0; 1.0])
        csmatout
    end

    # The vvolume integrals are evaluated using this rule
    gr = GaussRule(3, 2)

    # We will create two regions, one for the layers with 0°  orientation,
    # and one for the layers with 90° orientation.
    rl1 = vcat(selectelem(fens, fes, label = 1),
        selectelem(fens, fes, label = 3),
        selectelem(fens, fes, label = 5),
        selectelem(fens, fes, label = 7))
    rfes1 = subset(fes, rl1)
    region1 = FDataDict("femm" => FEMMDeforLinearMSH8(MR,
        IntegDomain(rfes1, gr),
        CSys(3,
            3,
            (csmatout, XYZ, tangents, feid, qpid) -> _updatecs!(csmatout,
                feid,
                rfes1.label)),
        material))
    rl2 = vcat(selectelem(fens, fes, label = 2),
        selectelem(fens, fes, label = 4),
        selectelem(fens, fes, label = 6))
    rfes2 = subset(fes, rl2)
    region2 = FDataDict("femm" => FEMMDeforLinearMSH8(MR,
        IntegDomain(rfes2, gr),
        CSys(3,
            3,
            (csmatout, XYZ, tangents, feid, qpid) -> _updatecs!(csmatout,
                feid,
                rfes2.label)),
        material))

    # File =  "NAFEMS-R0031-1-plate-r1.vtk"
    # vtkexportmesh(File, region1["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "NAFEMS-R0031-1-plate-r2.vtk"
    # vtkexportmesh(File, region2["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)

    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate = tolerance)
    ex0 = FDataDict("displacement" => 0.0, "component" => 1, "node_list" => lx0)
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box = [-Inf Inf 0.0 0.0 -Inf Inf], inflate = tolerance)
    ey0 = FDataDict("displacement" => 0.0, "component" => 2, "node_list" => ly0)
    # The transverse displacement is fixed along the line  passing through
    # point B. The nodes are fixed in the box along this line in the Z
    # direction.
    lz0 = selectnode(fens, box = [AB / 2 AB / 2 -Inf Inf -Inf Inf], inflate = tolerance)
    ez0 = FDataDict("displacement" => 0.0, "component" => 3, "node_list" => lz0)

    # The traction boundary condition is applied  along  the edge of the
    # mesh passing through point C at the top surface of the strip.   First
    # we extract the boundary of the hexahedral mesh.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # X = 0
    xl = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate = tolerance)
    # Now we extract the boundary  of these selected quadrilaterals
    bbfes = meshboundary(subset(bfes, xl))
    # …  And from these  we extract the ones at the top
    zl = selectelem(fens, bbfes, box = [0.0 0.0 -Inf Inf TH TH], inflate = tolerance)
    # Note that  we have to apply only half of the line load given that
    # were modeling  just one quarter of the geometry and we are splitting
    # the line load  with the symmetry plane X=0. Also note that the
    # quadrature rule is one-dimensional  since we are integrating along
    # a curve.
    Trac = FDataDict("traction_vector" => vec([0.0; 0.0; -q0 / 2]),
        "femm" => FEMMBase(IntegDomain(subset(bbfes, zl), GaussRule(1, 3))))

    modeldata = FDataDict("fens" => fens,
        "regions" => [region1, region2],
        "essential_bcs" => [ex0, ey0, ez0],
        "traction_bcs" => [Trac])
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS_R0031_1_msh8-plate")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

    u = modeldata["u"]
    geom = modeldata["geom"]

    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nE = selectnode(fens, box = [0.0 0.0 0.0 0.0 0.0 0.0], inflate = tolerance)
    nC = selectnode(fens, box = [0.0 0.0 0.0 0.0 TH TH], inflate = tolerance)
    nD = selectnode(fens, box = [CD CD 0.0 0.0 ts[1] ts[1]], inflate = tolerance)
    n0z = selectnode(fens, box = [0.0 0.0 0.0 0.0 0.0 TH], inflate = tolerance)
    ix = sortperm(geom.values[n0z, 3])
    # # println("ix = $(ix)")

    cdis = mean(u.values[nE, 3])
    # println("")
    # println("Normalized Center deflection: $(cdis/wEref)")

    extrap = :extraptrend
    # # extrap = :extrapmean
    inspectormeth = :averaging
    # extrap = :default
    # inspectormeth = :invdistance

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS_R0031_1_msh8-plate-sx",
        "quantity" => :Cauchy, "component" => 1, "outputcsys" => CSys(3),
        "nodevalmethod" => inspectormeth, "reportat" => extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    # println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")
    # println("Reference sx@E = $(sigma11Eref/phun("MPa")) [MPa]")
    @test norm(s.values[nE] .- 654.948 * phun("MPa")) / sigma11Eref < 1.0e-3

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS_R0031_1_msh8-plate-sxz",
        "quantity" => :Cauchy, "component" => 5, "outputcsys" => CSys(3),
        "nodevalmethod" => inspectormeth, "reportat" => extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    # println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")
    @test norm(s.values[nD] .- 3.65590 * phun("MPa")) / sigma13Dref < 1.0e-3
    s = modeldata["postprocessing"]["exported"][2]["field"]
    # println("sxz@D_2 = $(s.values[nD]/phun("MPa")) [MPa]")
    @test norm(s.values[nD] .- 4.14278 * phun("MPa")) / sigma13Dref < 1.0e-3
    # println("Reference sxz@D = $(sigma13Dref/phun("MPa")) [MPa]")

    # println("Done")
    true
end # NAFEMS_R0031_1

    NAFEMS_R0031_1_msh8()
    
end # module 
nothing


module NAFEMS_R0031_1_test_esniceh8
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Statistics: mean
using LinearAlgebra: norm, cross
using Test

function NAFEMS_R0031_1_esnice_h8()
    # println("""
    # Laminated Strip Under Three-Point Bending. ESNICE h8 elements. 
    # """)

    # Determine the central transverse displacement in a simply-supported seven
    # layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
    # material lay-up is specified with the center ply being four times as
    # thick as the others.
    # Reference: NAFEMS Report R0031, Test No.1, 17-Dec-1998.

    # Because of the symmetries of the geometry and load, only the
    # first-quadrant   (in XY) quarter of the plate is modeled.

    # The coordinate system is centered at point E (at the difference with
    # respect to the original benchmark definition).  The  load is applied
    # along a curve passing through point C. The simple support is applied
    # along the curve passing through point B.

    # We realize the simple supports along the lines  A, B and the line load at
    # point C  are illegal from the point of view of convergence.  No
    # convergence can be hoped for as the stress underneath the load and above
    # the simple supports  is infinite in the limit (these locations are stress
    # singularities).   However, for relatively coarse meshes the results away
    # from the singularities are still meaningful.

    # The target quantities are displacement at the bottom surface at point E,
    # the tensile axial stress at the same point,  and of the transverse shear
    # stress at point D  in between the bottommost two layers (See figure 1).

    t0 = time()
    # Orthotropic material parameters of the material of the layers
    E1s = 100.0 * phun("GPa")
    E2s = E3s = 5.0 * phun("GPa")
    nu12s = 0.4
    nu13s = 0.3
    nu23s = 0.3
    G12s = 3.0 * phun("GPa")
    G13s = G23s = 2.0 * phun("GPa")
    CTE1 = 3.0e-6
    CTE2 = 2.0e-5
    CTE3 = 2.0e-5

    AB = 30.0 * phun("mm") # span between simple supports
    CD = 4.0 * phun("mm") # distance between the point D and the center
    OH = 10.0 * phun("mm") # overhang
    W = 10.0 * phun("mm") # width of the strip

    # Here we define the layout and the thicknesses of the layers. Region 1 has
    # orientation 0 degrees, and region 2 has orientation 90 degrees.
    angles = vec([0 90])
    ts = vec([0.1 0.1 0.1 0.4 0.1 0.1 0.1]) * phun("mm") # layer thicknesses
    TH = sum(ts) # total thickness of the plate

    tolerance = 0.0001 * TH

    # The line load is in the negative Z direction.
    q0 = 10 * phun("N/mm") #    line load

    # Reference deflection under the load is
    wEref = -1.06 * phun("mm")

    # The reference tensile stress at the bottom of the lowest layer is
    sigma11Eref = 684 * phun("MPa")

    # Because we model the first-quadrant quarter of the plate using
    # coordinate axes centered  at the point E  the shear at the point D is
    # positive instead of negative as in the benchmark where the coordinate
    # system is located at the outer corner of the strip.
    sigma13Dref = 4.1 * phun("MPa")

    Refinement = 4
    # We select 8 elements spanwise and 2 elements widthwise.  The overhang
    # of the plate is given one element.
    nL = Refinement * 4
    nO = Refinement * 1
    nW = Refinement * 1

    # Each layer is modeled with a single element times the refinement.
    nts = Refinement * ones(Int, length(ts))# number of elements per layer

    xs = unique(vcat(collect(linearspace(0.0, AB / 2, nL + 1)), [CD],
        collect(linearspace(AB / 2, AB / 2 + OH, nO + 1))))
    xs = xs[sortperm(xs)]
    ys = collect(linearspace(0.0, W / 2, nW + 1))

    fens, fes = H8layeredplatex(xs, ys, ts, nts)

    # This is the material  model
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
        0.0, E1s, E2s, E3s,
        nu12s, nu13s, nu23s,
        G12s, G13s, G23s,
        CTE1, CTE2, CTE3)

    # The material coordinate system function is defined as:
    function _updatecs!(csmatout, feid, labels)
        rotmat3!(csmatout, angles[feid] / 180.0 * pi * [0.0; 0.0; 1.0])
        csmatout
    end

    # The vvolume integrals are evaluated using this rule
    gr = TrapezoidalRule(3)

    # We will create two regions, one for the layers with 0°  orientation,
    # and one for the layers with 90° orientation.
    rl1 = vcat(selectelem(fens, fes, label = 1),
        selectelem(fens, fes, label = 3),
        selectelem(fens, fes, label = 5),
        selectelem(fens, fes, label = 7))
    rfes1 = subset(fes, rl1)
    region1 = FDataDict("femm" => FEMMDeforLinearESNICEH8(MR,
        IntegDomain(rfes1, gr),
        CSys(3,
            3,
            (csmatout, XYZ, tangents, feid, qpid) -> _updatecs!(csmatout, 1, [])),
        material))
    rl2 = vcat(selectelem(fens, fes, label = 2),
        selectelem(fens, fes, label = 4),
        selectelem(fens, fes, label = 6))
    rfes2 = subset(fes, rl2)
    region2 = FDataDict("femm" => FEMMDeforLinearESNICEH8(MR,
        IntegDomain(rfes2, gr),
        CSys(3,
            3,
            (csmatout, XYZ, tangents, feid, qpid) -> _updatecs!(csmatout, 2, [])),
        material))

    # File =  "NAFEMS-R0031-1-plate-r1.vtk"
    # vtkexportmesh(File, region1["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "NAFEMS-R0031-1-plate-r2.vtk"
    # vtkexportmesh(File, region2["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)

    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate = tolerance)
    ex0 = FDataDict("displacement" => 0.0, "component" => 1, "node_list" => lx0)
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box = [-Inf Inf 0.0 0.0 -Inf Inf], inflate = tolerance)
    ey0 = FDataDict("displacement" => 0.0, "component" => 2, "node_list" => ly0)
    # The transverse displacement is fixed along the line  passing through
    # point B. The nodes are fixed in the box along this line in the Z
    # direction.
    lz0 = selectnode(fens, box = [AB / 2 AB / 2 -Inf Inf -Inf Inf], inflate = tolerance)
    ez0 = FDataDict("displacement" => 0.0, "component" => 3, "node_list" => lz0)

    # The traction boundary condition is applied  along  the edge of the
    # mesh passing through point C at the top surface of the strip.   First
    # we extract the boundary of the hexahedral mesh.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # X = 0
    xl = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate = tolerance)
    # Now we extract the boundary  of these selected quadrilaterals
    bbfes = meshboundary(subset(bfes, xl))
    # …  And from these  we extract the ones at the top
    zl = selectelem(fens, bbfes, box = [0.0 0.0 -Inf Inf TH TH], inflate = tolerance)
    # Note that  we have to apply only half of the line load given that
    # were modeling  just one quarter of the geometry and we are splitting
    # the line load  with the symmetry plane X=0. Also note that the
    # quadrature rule is one-dimensional  since we are integrating along
    # a curve.
    Trac = FDataDict("traction_vector" => vec([0.0; 0.0; -q0 / 2]),
        "femm" => FEMMBase(IntegDomain(subset(bbfes, zl), TrapezoidalRule(1))))

    modeldata = FDataDict("fens" => fens,
        "regions" => [region1, region2],
        "essential_bcs" => [ex0, ey0, ez0],
        "traction_bcs" => [Trac])
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS-R0031-1-plate-esnice-h8")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

    u = modeldata["u"]
    geom = modeldata["geom"]

    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nE = selectnode(fens, box = [0.0 0.0 0.0 0.0 0.0 0.0], inflate = tolerance)
    nC = selectnode(fens, box = [0.0 0.0 0.0 0.0 TH TH], inflate = tolerance)
    nD = selectnode(fens, box = [CD CD 0.0 0.0 ts[1] ts[1]], inflate = tolerance)
    n0z = selectnode(fens, box = [0.0 0.0 0.0 0.0 0.0 TH], inflate = tolerance)
    ix = sortperm(geom.values[n0z, 3])
    # # println("ix = $(ix)")

    cdis = mean(u.values[nE, 3])
    # println("")
    # println("Normalized Center deflection: $(cdis/wEref)")

    extrap = :extraptrend
    # # extrap = :extrapmean
    inspectormeth = :averaging
    # extrap = :default
    # inspectormeth = :invdistance

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS-R0031-1-plate-esnice-h8-sx",
        "quantity" => :Cauchy, "component" => 1, "outputcsys" => CSys(3),
        "nodevalmethod" => inspectormeth, "reportat" => extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    # println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")
    # println("Reference sx@E = $(sigma11Eref/phun("MPa")) [MPa]")
    @test norm(s.values[nE] .- 667.6985 * phun("MPa")) / sigma11Eref < 1.0e-3

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS-R0031-1-plate-esnice-h8-sxz",
        "quantity" => :Cauchy, "component" => 5, "outputcsys" => CSys(3),
        "nodevalmethod" => inspectormeth, "reportat" => extrap)
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    # println("sxz@D_1 = $(s.values[nD]/phun("MPa")) [MPa]")
    @test norm(s.values[nD] .- 3.75795 * phun("MPa")) / sigma13Dref < 1.0e-3
    s = modeldata["postprocessing"]["exported"][2]["field"]
    # println("sxz@D_2 = $(s.values[nD]/phun("MPa")) [MPa]")
    @test norm(s.values[nD] .- 4.36560 * phun("MPa")) / sigma13Dref < 1.0e-3
    # println("Reference sxz@D = $(sigma13Dref/phun("MPa")) [MPa]")

    # println("Done")
    true
end # NAFEMS_R0031_1

    NAFEMS_R0031_1_esnice_h8()
 
end # module 
nothing


module NAFEMS_R0031_1_test_h20
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Statistics: mean
using LinearAlgebra: norm, cross
using Test

function NAFEMS_R0031_1_H20()
    # println("""
    # Laminated Strip Under Three-Point Bending. H20 hexahedral elements.
    # """)

    # Determine the central transverse displacement in a simply-supported seven
    # layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
    # material lay-up is specified with the center ply being four times as
    # thick as the others.
    # Reference: NAFEMS Report R0031, Test No.1, 17-Dec-1998.

    # Because of the symmetries of the geometry and load, only the
    # first-quadrant   (in XY) quarter of the plate is modeled.

    # The coordinate system is centered at point E (at the difference with
    # respect to the original benchmark definition).  The  load is applied
    # along a curve passing through point C. The simple support is applied
    # along the curve passing through point B.

    # We realize the simple supports along the lines  A, B and the line load at
    # point C  are illegal from the point of view of convergence.  No
    # convergence can be hoped for as the stress underneath the load and above
    # the simple supports  is infinite in the limit (these locations are stress
    # singularities).   However, for relatively coarse meshes the results away
    # from the singularities are still meaningful.

    # The target quantities are displacement at the bottom surface at point E,
    # the tensile axial stress at the same point,  and of the transverse shear
    # stress at point D  in between the bottommost two layers (See figure 1).

    t0 = time()
    # Orthotropic material parameters of the material of the layers
    E1s = 100.0 * phun("GPa")
    E2s = E3s = 5.0 * phun("GPa")
    nu12s = 0.4
    nu13s = 0.3
    nu23s = 0.3
    G12s = 3.0 * phun("GPa")
    G13s = G23s = 2.0 * phun("GPa")
    CTE1 = 3.0e-6
    CTE2 = 2.0e-5
    CTE3 = 2.0e-5

    AB = 30.0 * phun("mm") # span between simple supports
    CD = 4.0 * phun("mm") # distance between the point D and the center
    OH = 10.0 * phun("mm") # overhang
    W = 10.0 * phun("mm") # width of the strip

    # Here we define the layout and the thicknesses of the layers.
    angles = vec([0 90 0 90 0 90 0])
    ts = vec([0.1 0.1 0.1 0.4 0.1 0.1 0.1]) * phun("mm") # layer thicknesses
    TH = sum(ts) # total thickness of the plate

    tolerance = 0.0001 * TH

    # The line load is in the negative Z direction.
    q0 = 10 * phun("N/mm") #    line load

    # Reference deflection under the load is
    wEref = -1.06 * phun("mm")

    # The reference tensile stress at the bottom of the lowest layer is
    sigma11Eref = 684 * phun("MPa")

    # Because we model the first-quadrant quarter of the plate using
    # coordinate axes centered  at the point E  the shear at the point D is
    # positive instead of negative as in the benchmark where the coordinate
    # system is located at the outer corner of the strip.
    sigma13Dref = 4.1 * phun("MPa")

    Refinement = 4
    # We select 8 elements spanwise and 2 elements widthwise.  The overhang
    # of the plate is given one element.
    nL = Refinement * 4
    nO = Refinement * 1
    nW = Refinement * 1

    # Each layer is modeled with a single element times the refinementp.
    nts = Refinement * ones(Int, length(angles))# number of elements per layer

    xs = unique(vcat(collect(linearspace(0.0, AB / 2, nL + 1)), [CD],
        collect(linearspace(AB / 2, AB / 2 + OH, nO + 1))))
    xs = xs[sortperm(xs)]
    ys = collect(linearspace(0.0, W / 2, nW + 1))

    fens, fes = H8layeredplatex(xs, ys, ts, nts)
    fens, fes = H8toH20(fens, fes)

    # This is the material  model
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
        0.0, E1s, E2s, E3s,
        nu12s, nu13s, nu23s,
        G12s, G13s, G23s,
        CTE1, CTE2, CTE3)

    # The material coordinate system function is defined as:
    function _updatecs!(csmatout, feid, labels)
        rotmat3!(csmatout, angles[labels[feid]] / 180.0 * pi * [0.0; 0.0; 1.0])
        csmatout
    end

    # The vvolume integrals are evaluated using this rule
    gr = GaussRule(3, 3)

    # We will create two regions, one for the layers with 0°  orientation,
    # and one for the layers with 90° orientation.
    rl1 = vcat(selectelem(fens, fes, label = 1),
        selectelem(fens, fes, label = 3),
        selectelem(fens, fes, label = 5),
        selectelem(fens, fes, label = 7))
    rfes1 = subset(fes, rl1)
    region1 = FDataDict("femm" => FEMMDeforLinear(MR,
        IntegDomain(rfes1, gr),
        CSys(3,
            3,
            (csmatout, XYZ, tangents, feid, qpid) -> _updatecs!(csmatout,
                feid,
                rfes1.label)),
        material))
    rl2 = vcat(selectelem(fens, fes, label = 2),
        selectelem(fens, fes, label = 4),
        selectelem(fens, fes, label = 6))
    rfes2 = subset(fes, rl2)
    region2 = FDataDict("femm" => FEMMDeforLinear(MR,
        IntegDomain(rfes2, gr),
        CSys(3,
            3,
            (csmatout, XYZ, tangents, feid, qpid) -> _updatecs!(csmatout,
                feid,
                rfes2.label)),
        material))

    # File =  "NAFEMS-R0031-1-plate-r1.vtk"
    # vtkexportmesh(File, region1["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # # @async run(`"paraview.exe" $File`)
    # File =  "NAFEMS-R0031-1-plate-r2.vtk"
    # vtkexportmesh(File, region2["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.H8)
    # @async run(`"paraview.exe" $File`)

    # The essential boundary conditions are applied on the symmetry planes.
    # First the plane X=0;...
    lx0 = selectnode(fens, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate = tolerance)
    ex0 = FDataDict("displacement" => 0.0, "component" => 1, "node_list" => lx0)
    # ... and then the plane Y=0.
    ly0 = selectnode(fens, box = [-Inf Inf 0.0 0.0 -Inf Inf], inflate = tolerance)
    ey0 = FDataDict("displacement" => 0.0, "component" => 2, "node_list" => ly0)
    # The transverse displacement is fixed along the line  passing through
    # point B. The nodes are fixed in the box along this line in the Z
    # direction.
    lz0 = selectnode(fens, box = [AB / 2 AB / 2 -Inf Inf 0.0 0.0], inflate = tolerance)
    ez0 = FDataDict("displacement" => 0.0, "component" => 3, "node_list" => lz0)

    # The traction boundary condition is applied  along  the edge of the
    # mesh passing through point C at the top surface of the strip.   First
    # we extract the boundary of the hexahedral mesh.
    bfes = meshboundary(fes)
    # From  the entire boundary we select those quadrilaterals that lie on the plane
    # X = 0
    xl = selectelem(fens, bfes, box = [0.0 0.0 -Inf Inf -Inf Inf], inflate = tolerance)
    # Now we extract the boundary  of these selected quadrilaterals
    bbfes = meshboundary(subset(bfes, xl))
    # …  And from these  we extract the ones at the top
    zl = selectelem(fens, bbfes, box = [0.0 0.0 -Inf Inf TH TH], inflate = tolerance)
    # Note that  we have to apply only half of the line load given that
    # were modeling  just one quarter of the geometry and we are splitting
    # the line load  with the symmetry plane X=0. Also note that the
    # quadrature rule is one-dimensional  since we are integrating along
    # a curve.
    Trac = FDataDict("traction_vector" => vec([0.0; 0.0; -q0 / 2]),
        "femm" => FEMMBase(IntegDomain(subset(bbfes, zl), GaussRule(1, 3))))

    modeldata = FDataDict("fens" => fens,
        "regions" => [region1, region2],
        "essential_bcs" => [ex0, ey0, ez0],
        "traction_bcs" => [Trac])
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS_R0031_1_H20-plate")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

    u = modeldata["u"]
    geom = modeldata["geom"]

    # The results of the displacement and stresses will be reported at
    # nodes located at the appropriate points.
    nE = selectnode(fens, box = [0.0 0.0 0.0 0.0 0.0 0.0], inflate = tolerance)
    nC = selectnode(fens, box = [0.0 0.0 0.0 0.0 TH TH], inflate = tolerance)
    nD = selectnode(fens, box = [CD CD 0.0 0.0 ts[1] ts[1]], inflate = tolerance)

    cdis = mean(u.values[nE, 3])
    # println("")
    # println("Normalized Center deflection: $(cdis/wEref)")

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS_R0031_1_H20-plate-sx",
        "quantity" => :Cauchy, "component" => 1, "outputcsys" => CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    # println("sx@E = $(s.values[nE]/phun("MPa")) [MPa]")
    # println("Reference sx@E = $(sigma11Eref/phun("MPa")) [MPa]")
    @test norm(s.values[nE] .- 668.16037 * phun("MPa")) / sigma11Eref < 1.0e-3

    modeldata["postprocessing"] = FDataDict("file" => "NAFEMS_R0031_1_H20-plate-sxz",
        "quantity" => :Cauchy, "component" => 5, "outputcsys" => CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    s = modeldata["postprocessing"]["exported"][1]["field"]
    # println("sxz@D = $(s.values[nD]/phun("MPa")) [MPa]")
    @test norm(s.values[nD] .- 3.619757 * phun("MPa")) / sigma13Dref < 1.0e-3
    s = modeldata["postprocessing"]["exported"][2]["field"]
    # println("sxz@D = $(s.values[nD]/phun("MPa")) [MPa]")
    @test norm(s.values[nD] .- 4.100860 * phun("MPa")) / sigma13Dref < 1.0e-3
    # println("Reference sxz@D = $(sigma13Dref/phun("MPa")) [MPa]")

    # println("Done")
    true
end # NAFEMS_R0031_1_H20

    NAFEMS_R0031_1_H20()

end # module 
nothing
