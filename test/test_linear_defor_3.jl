
module scratch1_06092017_ortho
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using Test

mutable struct MyIData
    c::Int
    r::Vector{Float64}
    s::Vector{Float64}
end

function test()

    # println("Thick pipe with internal pressure: axially symmetric model")
    #=
    This is a simple modification of the full three-dimensional simulation of
    the tutorial pub_thick_pipe that implements the axially-symmetric model
    reduction procedure.

    An infinitely long thick walled cylindrical pipe
    with inner boundary radius of 3 mm and outer boundary radius of 9 mm is
    subjected to an internal pressure of 1.0 MPa. A wedge   with thickness of
    2 mm and a 90-degree angle sector is considered for the finite element
    analysis. The material properties are taken as  isotropic linear elastic
    with $E=1000$ MPa and $\nu=0.4999$ to represent nearly incompressible
    behavior. This problem has been proposed to by MacNeal and Harder as a
    test of an element's ability to represent the  response of a nearly
    incompressible material. The plane-strain condition is assumed in the
    axial direction of the pipe which together with the radial symmetry
    confines the material in all but the radial direction and therefore
    amplifies the numerical difficulties associated with the confinement of
    the nearly incompressible material.

    There is an analytical solution to this problem. Timoshenko and Goodier
    presented the original solution of Lame in their textbook. We are going
    to compare with  both the stress distribution (radial and hoop stresses)
    and the displacement of the inner  cylindrical surface.

    References:
    - Macneal RH, Harder RL (1985) A proposed standard set of problems to test
    finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
    - Timoshenko S. and Goodier J. N., Theory of Elasticity, McGraw-Hill, 2nd ed., 1951.

    =#

    # Internal radius of the pipe.
    a = 3 * phun("MM")
    ##
    # External radius of the pipe.
    b = 9 * phun("MM")
    ##
    # Thickness of the slice.
    t = 2 * phun("MM")

    ##
    # Geometrical tolerance.
    tolerance = a / 10000.0
    ##
    # Young's modulus and Poisson's ratio.
    E = 1000 * phun("MEGA*PA")
    nu = 0.499
    ##
    # Applied pressure on the internal surface.
    press = 1.0 * phun("MEGA*PA")

    ##
    # Analytical solutions.   Radial stress:
    radial_stress(r) = press * a .^ 2 / (b^2 - a^2) .* (1 - (b^2) ./ r .^ 2)
    ##
    # Circumferential (hoop) stress:
    hoop_stress(r) = press * a .^ 2 / (b^2 - a^2) .* (1 + (b^2) ./ r .^ 2)

    ##
    # Radial displacement:
    radial_displacement(r) =
        press * a^2 * (1 + nu) * (b^2 + r .^ 2 * (1 - 2 * nu)) / (E * (b^2 - a^2) .* r)

    ##
    # Therefore the radial displacement of the loaded surface will be:
    urex = radial_displacement(a)


    ##
    # The mesh parameters: The numbers of element edges axially,
    # and through the thickness of the pipe wall (radially).

    na = 1
    nt = 10

    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in effect.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm
    axisymmetric = true

    # Create the mesh and initialize the geometry.  First we are going
    # to construct the block of elements with the first coordinate
    # corresponding to the thickness in the radial direction, and the second
    # coordinate is the thickness in the axial direction.
    fens, fes = Q8block(b - a, t, nt, na)

    # Extract the boundary  and mark the finite elements on the
    # interior surface.
    bdryfes = meshboundary(fes)

    bcl = selectelem(fens, bdryfes, box = [0.0, 0.0, -Inf, Inf], inflate = tolerance)
    internal_fenids = connectednodes(subset(bdryfes, bcl))
    # Now  shape the block  into  the actual wedge piece of the pipe.
    for i = 1:count(fens)
        fens.xyz[i, :] = fens.xyz[i, :] + [a; 0.0]
    end

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # The plane-strain condition in the axial direction  is specified by selecting nodes
    # on the plane y=0 and y=t.
    l1 = selectnode(fens; box = [-Inf Inf 0.0 0.0], inflate = tolerance)
    setebc!(u, l1, true, 2, 0.0)
    l1 = selectnode(fens; box = [-Inf Inf t t], inflate = tolerance)
    setebc!(u, l1, true, 2, 0.0)

    applyebc!(u)
    numberdofs!(u)

    # The traction boundary condition is applied in the radial
    # direction.

    el1femm = FEMMBase(IntegDomain(subset(bdryfes, bcl), GaussRule(1, 3), axisymmetric))
    fi = ForceIntensity([press; 0.0])
    F2 = distribloads(el1femm, geom, u, fi, 2)

    # Property and material
    material = MatDeforElastOrtho(MR, E, nu)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), axisymmetric), material)

    K = stiffness(femm, geom, u)

    u = solve_blocked!(u, K, F2)

    # Transfer the solution of the displacement to the nodes on the
    # internal cylindrical surface and convert to
    # cylindrical-coordinate displacements there.
    uv = u.values[internal_fenids, :]
    # Report the  relative displacement on the internal surface:
    # println("(Approximate/true displacement) at the internal surface: $( mean(uv[:,1])/urex*100  ) %")

    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)


    # File =  "thick_pipe_sigmax.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])

    # Produce a plot of the solution components in the cylindrical
    # coordinate system.

    function inspector(idat::MyIData, elnum, conn, xe, out, xq)
        push!(idat.r, xq[1])
        push!(idat.s, out[idat.c])
        return idat
    end

    idat = MyIData(1, Int[], Int[])
    idat =
        inspectintegpoints(femm, geom, u, collect(1:count(fes)), inspector, idat, :Cauchy)

    # using Plots
    # plotly()
    #
    # # Plot the analytical solution.
    # r = linearspace(a,b,100);
    # plot(r, radial_stress(r))
    # # Plot the computed  integration-point data
    # plot!(idat.r, idat.s, m=:circle, color=:red)
    # gui()

    @test abs(idat.r[1] - 0.003126794919243112) < 1.0e-9
    @test abs(idat.s[1] - -910911.9777008593) < 1.0e-2

    ## Discussion
    #
    ##
    # The axially symmetric model is clearly very effective
    # computationally, as the size is much reduced compared to the 3-D
    # model.  In conjunction with  uniform or selective reduced integration
    # it can be very accurate as well.
end
end
using .scratch1_06092017_ortho
scratch1_06092017_ortho.test()

module mmLE11Q8mm
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea = 210000 * phun("MEGA*Pa")
    nua = 0.3
    alphaa = 2.3e-4              # thermal expansion coefficient
    sigmaA = -105 * phun("MEGA*Pa")
    nref = 1                        # how many times should we refine the mesh?
    X =
        [
            1.0 0.0#A
            1.4 0.0#B
            0.995184726672197 0.098017140329561
            1.393258617341076 0.137223996461385
            0.980785 0.195090#
            1.37309939 0.27312645
            0.956940335732209 0.290284677254462
            1.339716470025092 0.406398548156247
            0.9238795 0.38268#C
            1.2124 0.7#D
            0.7071 0.7071#E
            1.1062 1.045#F
            0.7071 (0.7071+1.79)/2#(E+H)/2
            1.0 1.39#G
            0.7071 1.79#H
            1.0 1.79#I
        ] * phun("M")
    tolerance = 1.e-6 * phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm



    fens = FENodeSet(X)
    fes =
        FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15])
    for ref = 1:nref
        fens, fes = Q4refine(fens, fes)
        list = selectnode(
            fens,
            distance = 1.0 + 0.1 / 2^nref,
            from = [0.0 0.0],
            inflate = tolerance,
        )
        fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)
    end
    fens, fes = Q4toQ8(fens, fes)
    list = selectnode(
        fens,
        distance = 1.0 + 0.1 / 2^nref,
        from = [0.0 0.0],
        inflate = tolerance,
    )
    fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)

    #     File  =   "mesh.vtk"
    # vtkexportmesh(File, fens, fes)

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply EBC's
    l1 = selectnode(fens, box = [-Inf Inf 0 0], inflate = tolerance)
    setebc!(u, l1, true, 2, 00.0)
    l1 = selectnode(fens, box = [-Inf Inf 1.79 1.79], inflate = tolerance)
    setebc!(u, l1, true, 2, 00.0)
    applyebc!(u)
    numberdofs!(u)

    # Temperature field
    dT = NodalField(reshape(fens.xyz[:, 1] + fens.xyz[:, 2], size(fens.xyz, 1), 1))


    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

    K = stiffness(femm, geom, u)
    F = thermalstrainloads(femm, geom, u, dT)

    u = solve_blocked!(u, K, F)

    nA = selectnode(fens, box = Float64[1.0 1.0 0.0 0.0], inflate = tolerance)

    fld = fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


    File = "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("sigmay", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-1.443052182185006e8, -1.4106181545272605e7],
    ) < 1.0e-2
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    sA = fld.values[nA] / phun("MEGA*Pa")
    sAn = fld.values[nA] / sigmaA
    # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    @test abs(sA[1] - (-93.8569)) < 1.0e-3

    fen2fe = FENodeToFEMap(fes, nnodes(geom))
    function inspector(idat, elnum, conn, xe, out, xq)
        println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end

    # inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
    # inspector, []; quantity = :Cauchy)

end
end
using .mmLE11Q8mm
mmLE11Q8mm.test()

module mmLE11Q8mmortho
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea = 210000 * phun("MEGA*Pa")
    nua = 0.3
    alphaa = 2.3e-4              # thermal expansion coefficient
    sigmaA = -105 * phun("MEGA*Pa")
    nref = 1                        # how many times should we refine the mesh?
    X =
        [
            1.0 0.0#A
            1.4 0.0#B
            0.995184726672197 0.098017140329561
            1.393258617341076 0.137223996461385
            0.980785 0.195090#
            1.37309939 0.27312645
            0.956940335732209 0.290284677254462
            1.339716470025092 0.406398548156247
            0.9238795 0.38268#C
            1.2124 0.7#D
            0.7071 0.7071#E
            1.1062 1.045#F
            0.7071 (0.7071+1.79)/2#(E+H)/2
            1.0 1.39#G
            0.7071 1.79#H
            1.0 1.79#I
        ] * phun("M")
    tolerance = 1.e-6 * phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm



    fens = FENodeSet(X)
    fes =
        FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15])
    for ref = 1:nref
        fens, fes = Q4refine(fens, fes)
        list = selectnode(
            fens,
            distance = 1.0 + 0.1 / 2^nref,
            from = [0.0 0.0],
            inflate = tolerance,
        )
        fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)
    end
    fens, fes = Q4toQ8(fens, fes)
    list = selectnode(
        fens,
        distance = 1.0 + 0.1 / 2^nref,
        from = [0.0 0.0],
        inflate = tolerance,
    )
    fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)

    #     File  =   "mesh.vtk"
    # vtkexportmesh(File, fens, fes)

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply EBC's
    l1 = selectnode(fens, box = [-Inf Inf 0 0], inflate = tolerance)
    setebc!(u, l1, true, 2, 00.0)
    l1 = selectnode(fens, box = [-Inf Inf 1.79 1.79], inflate = tolerance)
    setebc!(u, l1, true, 2, 00.0)
    applyebc!(u)
    numberdofs!(u)

    # Temperature field
    dT = NodalField(reshape(fens.xyz[:, 1] + fens.xyz[:, 2], size(fens.xyz, 1), 1))


    # Property and material
    material = MatDeforElastOrtho(MR, 0.0, Ea, nua, alphaa)
    # display(material )

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

    K = stiffness(femm, geom, u)
    F = thermalstrainloads(femm, geom, u, dT)

    u = solve_blocked!(u, K, F)

    nA = selectnode(fens, box = Float64[1.0 1.0 0.0 0.0], inflate = tolerance)

    fld = fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


    File = "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("sigmay", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-1.443052182185006e8, -1.4106181545272605e7],
    ) < 1.0e-2
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    sA = fld.values[nA] / phun("MEGA*Pa")
    sAn = fld.values[nA] / sigmaA
    # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    @test abs(sA[1] - (-93.8569)) < 1.0e-3

    fen2fe = FENodeToFEMap(fes, nnodes(geom))
    function inspector(idat, elnum, conn, xe, out, xq)
        println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end

    # inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
    # inspector, []; quantity = :Cauchy)

end
end
using .mmLE11Q8mmortho
mmLE11Q8mmortho.test()

module mLE11Q8aximmm
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea = 210000 * phun("MEGA*Pa")
    nua = 0.3
    alphaa = 2.3e-4              # thermal expansion coefficient
    sigmaA = -105 * phun("MEGA*Pa")
    nref = 2                        # how many times should we refine the mesh?
    X =
        [
            1.0 0.0#A
            1.4 0.0#B
            0.995184726672197 0.098017140329561
            1.393258617341076 0.137223996461385
            0.980785 0.195090#
            1.37309939 0.27312645
            0.956940335732209 0.290284677254462
            1.339716470025092 0.406398548156247
            0.9238795 0.38268#C
            1.2124 0.7#D
            0.7071 0.7071#E
            1.1062 1.045#F
            0.7071 (0.7071+1.79)/2#(E+H)/2
            1.0 1.39#G
            0.7071 1.79#H
            1.0 1.79#I
        ] * phun("M")
    tolerance = 1.e-6 * phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm



    fens = FENodeSet(X)
    fes =
        FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15])
    for ref = 1:nref
        fens, fes = Q4refine(fens, fes)
        list = selectnode(
            fens,
            distance = 1.0 + 0.1 / 2^nref,
            from = [0.0 0.0],
            inflate = tolerance,
        )
        fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)
    end
    fens, fes = Q4toQ8(fens, fes)
    list = selectnode(
        fens,
        distance = 1.0 + 0.1 / 2^nref,
        from = [0.0 0.0],
        inflate = tolerance,
    )
    fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)

    #     File  =   "mesh.vtk"
    # vtkexportmesh(File, fens, fes)

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply EBC's
    lbottom = selectnode(fens, box = [-Inf Inf 0 0], inflate = tolerance)
    setebc!(u, lbottom, true, 2, 00.0)
    ltop = selectnode(fens, box = [-Inf Inf 1.79 1.79], inflate = tolerance)
    setebc!(u, ltop, true, 2, 00.0)
    applyebc!(u)
    numberdofs!(u)

    # Temperature field
    dT = NodalField(reshape(fens.xyz[:, 1] + fens.xyz[:, 2], size(fens.xyz, 1), 1))


    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

    K = stiffness(femm, geom, u)
    F = thermalstrainloads(femm, geom, u, dT)

    u = solve_blocked!(u, K, F)

    nA = selectnode(fens, box = Float64[1.0 1.0 0.0 0.0], inflate = tolerance)

    fld = fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


    File = "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("sigmay", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-1.6338426447540134e8, -4.961956343464769e6],
    ) < 1.e-1
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end


    sA = fld.values[nA] / phun("MEGA*Pa")
    sAn = fld.values[nA] / sigmaA
    # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")

    fen2fe = FENodeToFEMap(fes, nnodes(geom))
    function inspector(idat, elnum, conn, xe, out, xq)
        # println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end

    inspectintegpoints(
        femm,
        geom,
        u,
        dT,
        fen2fe.map[nA[1]],
        inspector,
        [];
        quantity = :Cauchy,
    )

    fld = fieldfromintegpoints(femm, geom, u, dT, :Pressure, 1)
    File = "LE11NAFEMS_Q8_pressure.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("pressure", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of  pressure = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-1.1881819144904878e7, 7.555030948761216e7],
    ) < 1.e-1
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    fld = fieldfromintegpoints(femm, geom, u, dT, :vm, 1)
    File = "LE11NAFEMS_Q8_vm.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("pressure", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of von Mises = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [3.221370699152578e7, 1.4437590830351183e8],
    ) < 1.e-1
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    AE = AbaqusExporter("LE11NAFEMS_Q8_export_stress")
    HEADING(AE, "NAFEMS LE11 benchmark with Q8 elements.")
    COMMENT(AE, "sigmaA = -105 MPa ")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "cax8", "AllElements", 1, connasarray(fes))
    NSET_NSET(AE, "ltop", ltop)
    NSET_NSET(AE, "lbottom", lbottom)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements")
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, Ea, nua)
    EXPANSION(AE, alphaa)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.ltop", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.lbottom", 2)
    TEMPERATURE(AE, "ASSEM1.INSTNC1.", collect(1:count(fens)), vec(dT.values))
    END_STEP(AE)
    close(AE)
    nlines = 0
    open(AE.filename) do f
        s = readlines(f)
        nlines = length(s)
    end
    @test nlines == 963
    try
        rm(AE.filename)
    catch
    end

end
end
using .mLE11Q8aximmm
mLE11Q8aximmm.test()


module mLE11Q8aximorthom
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    # NAFEMS LE11 benchmark with Q8 elements.
    # # This is a test recommended by the National Agency for Finite Element
    # # Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
    # # Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
    # #
    # # Target solution: Direct stress,   =  –105 MPa at point A.
    #function  LE11NAFEMS()
    # Parameters:
    Ea = 210000 * phun("MEGA*Pa")
    nua = 0.3
    alphaa = 2.3e-4              # thermal expansion coefficient
    sigmaA = -105 * phun("MEGA*Pa")
    nref = 2                        # how many times should we refine the mesh?
    X =
        [
            1.0 0.0#A
            1.4 0.0#B
            0.995184726672197 0.098017140329561
            1.393258617341076 0.137223996461385
            0.980785 0.195090#
            1.37309939 0.27312645
            0.956940335732209 0.290284677254462
            1.339716470025092 0.406398548156247
            0.9238795 0.38268#C
            1.2124 0.7#D
            0.7071 0.7071#E
            1.1062 1.045#F
            0.7071 (0.7071+1.79)/2#(E+H)/2
            1.0 1.39#G
            0.7071 1.79#H
            1.0 1.79#I
        ] * phun("M")
    tolerance = 1.e-6 * phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm



    fens = FENodeSet(X)
    fes =
        FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15])
    for ref = 1:nref
        fens, fes = Q4refine(fens, fes)
        list = selectnode(
            fens,
            distance = 1.0 + 0.1 / 2^nref,
            from = [0.0 0.0],
            inflate = tolerance,
        )
        fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)
    end
    fens, fes = Q4toQ8(fens, fes)
    list = selectnode(
        fens,
        distance = 1.0 + 0.1 / 2^nref,
        from = [0.0 0.0],
        inflate = tolerance,
    )
    fens.xyz[list, :] = FinEtools.MeshUtilModule.ontosphere(fens.xyz[list, :], 1.0)

    #     File  =   "mesh.vtk"
    # vtkexportmesh(File, fens, fes)

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    # Apply EBC's
    lbottom = selectnode(fens, box = [-Inf Inf 0 0], inflate = tolerance)
    setebc!(u, lbottom, true, 2, 00.0)
    ltop = selectnode(fens, box = [-Inf Inf 1.79 1.79], inflate = tolerance)
    setebc!(u, ltop, true, 2, 00.0)
    applyebc!(u)
    numberdofs!(u)

    # Temperature field
    dT = NodalField(reshape(fens.xyz[:, 1] + fens.xyz[:, 2], size(fens.xyz, 1), 1))


    # Property and material
    material = MatDeforElastOrtho(MR, 0.0, Ea, nua, alphaa)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

    K = stiffness(femm, geom, u)
    F = thermalstrainloads(femm, geom, u, dT)

    u = solve_blocked!(u, K, F)

    nA = selectnode(fens, box = Float64[1.0 1.0 0.0 0.0], inflate = tolerance)

    fld = fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


    File = "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("sigmay", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-1.6338426447540134e8, -4.961956343464769e6],
    ) < 1.e-1
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end


    sA = fld.values[nA] / phun("MEGA*Pa")
    sAn = fld.values[nA] / sigmaA
    # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")

    fen2fe = FENodeToFEMap(fes, nnodes(geom))
    function inspector(idat, elnum, conn, xe, out, xq)
        # println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
        return idat
    end

    inspectintegpoints(
        femm,
        geom,
        u,
        dT,
        fen2fe.map[nA[1]],
        inspector,
        [];
        quantity = :Cauchy,
    )

    fld = fieldfromintegpoints(femm, geom, u, dT, :Pressure, 1)
    File = "LE11NAFEMS_Q8_pressure.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("pressure", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of  pressure = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-1.1881819144904878e7, 7.555030948761216e7],
    ) < 1.e-1
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    fld = fieldfromintegpoints(femm, geom, u, dT, :vm, 1)
    File = "LE11NAFEMS_Q8_vm.vtk"
    vtkexportmesh(
        File,
        fens,
        fes;
        scalars = [("pressure", fld.values)],
        vectors = [("u", u.values)],
    )
    # println("range of von Mises = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [3.221370699152578e7, 1.4437590830351183e8],
    ) < 1.e-1
    # @async run(`"paraview.exe" $File`)
    try
        rm(File)
    catch
    end

    AE = AbaqusExporter("LE11NAFEMS_Q8_export_stress")
    HEADING(AE, "NAFEMS LE11 benchmark with Q8 elements.")
    COMMENT(AE, "sigmaA = -105 MPa ")
    PART(AE, "part1")
    END_PART(AE)
    ASSEMBLY(AE, "ASSEM1")
    INSTANCE(AE, "INSTNC1", "PART1")
    NODE(AE, fens.xyz)
    ELEMENT(AE, "cax8", "AllElements", 1, connasarray(fes))
    NSET_NSET(AE, "ltop", ltop)
    NSET_NSET(AE, "lbottom", lbottom)
    ORIENTATION(AE, "GlobalOrientation", vec([1.0 0 0]), vec([0 1.0 0]))
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements")
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE)
    END_ASSEMBLY(AE)
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, Ea, nua)
    EXPANSION(AE, alphaa)
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1.ltop", 2)
    BOUNDARY(AE, "ASSEM1.INSTNC1.lbottom", 2)
    TEMPERATURE(AE, "ASSEM1.INSTNC1.", collect(1:count(fens)), vec(dT.values))
    END_STEP(AE)
    close(AE)
    nlines = 0
    open(AE.filename) do f
        s = readlines(f)
        nlines = length(s)
    end
    @test nlines == 963
    try
        rm(AE.filename)
    catch
    end

end
end
using .mLE11Q8aximorthom
mLE11Q8aximorthom.test()

module mmmCookmmstrainmmisommm
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
    MR = DeforModelRed2DStrain
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
    @test abs(theutip[2] - 21.35955642390279) < 1.0e-3

    modeldata["postprocessing"] =
        FDataDict("file" => "cookstress-ew", "quantity" => :Cauchy, "component" => :xy)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.07028875155067116, 0.1301698279821655],
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
        [0.007503804468283987, 0.33798754356331173],
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
        [-0.11777749217431245, 0.23457099031101358],
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
        [-0.16098150217425994, 0.24838761904231466],
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
        "component" => 3,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.4523049370669106, 0.02980811337406548],
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
    COMMENT(AE, "CPE3 are pretty poor-accuracy elements, but here we don't care about it.")
    @test nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
    ELEMENT(
        AE,
        "CPE3",
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
    for ix  in eachindex(bfes)
        for jx = 1:2
            F[bfes.conn[ix][jx]] += 1.0 / n / 2 / thickness
        end
    end
    for ixxxx in eachindex(F)
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
using .mmmCookmmstrainmmisommm
mmmCookmmstrainmmisommm.test()

module mmmCookmmstrainmorthommm
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
    MR = DeforModelRed2DStrain
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
    @test abs(theutip[2] - 21.35955642390279) < 1.0e-3

    modeldata["postprocessing"] =
        FDataDict("file" => "cookstress-ew", "quantity" => :Cauchy, "component" => :xy)
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of Cauchy_xy = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.07028875155067116, 0.1301698279821655],
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
        [0.007503804468283987, 0.33798754356331173],
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
        [-0.11777749217431245, 0.23457099031101358],
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
        [-0.16098150217425994, 0.24838761904231466],
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
        "component" => 3,
    )
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    fld = modeldata["postprocessing"]["exported"][1]["field"]
    # println("range of princCauchy Min = $((minimum(fld.values), maximum(fld.values)))")
    @test norm(
        [minimum(fld.values), maximum(fld.values)] -
        [-0.4523049370669106, 0.02980811337406548],
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
    COMMENT(AE, "CPE3 are pretty poor-accuracy elements, but here we don't care about it.")
    @test nodesperelem(modeldata["regions"][1]["femm"].integdomain.fes) == 3
    ELEMENT(
        AE,
        "CPE3",
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
    for ixxxx in eachindex(F)
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
using .mmmCookmmstrainmorthommm
mmmCookmmstrainmorthommm.test()
