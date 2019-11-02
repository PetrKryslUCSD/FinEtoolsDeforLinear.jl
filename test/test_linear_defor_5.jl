
module mfiber_reinf_cant_yn_strong_Abaqus
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Test
import LinearAlgebra: Symmetric, cholesky
import Statistics: mean
function test()


# println("""
# Cantilever example.  Strongly orthotropic material. Orientation "y".
# @article{
# author = {Krysl, P.},
# title = {Mean-strain 8-node hexahedron with optimized energy-sampling stabilization},
# journal = {Finite Elements in Analysis and Design},
# volume = {108}, pages = {41-53}, DOI = {10.1016/j.finel.2015.09.008}, year = {2016}
# }
# """)

t0 = time()
pu = ustring -> phun(ustring; system_of_units = :SIMM)

# # Orthotropic material
E1s = 100000.0*pu("GPa")
E2s = 1.0*pu("GPa")
E3s = E2s
nu23s = nu12s = nu13s = 0.25
G12s = 0.2*pu("GPa")
G23s = G13s = G12s
CTE1 = 0.0
CTE2 = 0.0
CTE3 = 0.0
# # Isotropic material
# E = 1.0e9*pu("Pa")
# nu = 0.25
# CTE = 0.0

# Reference value for  the vertical deflection of the tip
uz_ref = -1.027498445054843e-05*pu("m");

a = 90.0*pu("mm") # length of the cantilever
b = 10.0*pu("mm") # width of the cross-section
t = 20.0*pu("mm") # height of the cross-section
q0 = -1000.0*pu("Pa") # shear traction
dT = 0*pu("K") # temperature rise

tolerance = 0.00001*t

# Generate mesh
n = 2
na = 4*n # number of elements lengthwise
nb = 2*n # number of elements through the depth
nt = n # number of elements through the thickness
xs = collect(linearspace(0.0, a, na+1))
ys = collect(linearspace(0.0, b, nb+1))
ts = collect(linearspace(0.0, t, nt+1))
# println("Mesh generation")
fens,fes = H8blockx(xs, ys, ts)
fens,fes = H8toH20(fens,fes)

bfes = meshboundary(fes)
# end cross-section surface  for the shear loading
sshearl = selectelem(fens, bfes; facing=true, direction = [+1.0 0.0 0.0])

MR = DeforModelRed3D
material = MatDeforElastOrtho(MR,
  0.0, E1s, E2s, E3s,
  nu12s, nu13s, nu23s,
  G12s, G13s, G23s,
  CTE1, CTE2, CTE3)
# material = MatDeforElastIso(MR,
#   0.0, E, nu, CTE)

# Material orientation matrix
csmat = zeros(3, 3)
rotmat3!(csmat, -45.0/180.0*pi*[0,1,0])

function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copyto!(csmatout, csmat)
end

femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), CSys(3, 3, updatecs!), material)

lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
nnodes(geom)

setebc!(u, lx0, true, 1, zeros(size(lx0)))
setebc!(u, lx0, true, 2, zeros(size(lx0)))
setebc!(u, lx0, true, 3, zeros(size(lx0)))
applyebc!(u)

# S = connectionmatrix(femm.femmbase, nnodes(geom))
numberdofs!(u)

function getshr!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copyto!(forceout, q0*[0.0; 0.0; 1.0])
end

Tracfemm = FEMMBase(IntegDomain(subset(bfes, sshearl), GaussRule(2, 3)))

# println("Stiffness")
K = stiffness(femm, geom, u)
fi = ForceIntensity(FFlt, 3, getshr!);
# println("Traction loads")
F =  distribloads(Tracfemm, geom, u, fi, 2);

# println("Factorization")
K = (K + K')/2;
K = cholesky(Symmetric(K))
# println("U = K\\F")
U = K\F
# # println("U = cg(K, F; tol=1e-3, maxiter=2000)")
# U = cg(K, F; tol=1e-3, maxiter=2000)
scattersysvec!(u, U[:])

Tipl = selectnode(fens, box=[a a b b 0. 0.], inflate=tolerance)
utip = mean(u.values[Tipl, 3])
# println("Deflection $utip, normalized: $(utip/uz_ref)")
@test abs(utip - -0.00653888266072445)/abs(utip) < 1.0e-6
# println("Solution: $(  time()-t0 )")

AE = AbaqusExporter("fiber_reinf_cant_yn_strong_Abaqus");
HEADING(AE, "fiber_reinf_cant_yn_strong_Abaqus.");
PART(AE, "part1");
END_PART(AE);
ASSEMBLY(AE, "ASSEM1");
INSTANCE(AE, "INSTNC1", "PART1");
NODE(AE, fens.xyz);
ELEMENT(AE, "c3d20r", "AllElements", connasarray(femm.integdomain.fes))
ELEMENT(AE, "SFM3D8", "TractionElements", connasarray(Tracfemm.integdomain.fes))
NSET_NSET(AE, "L1", lx0)
NSET_NSET(AE, "L2", lx0)
NSET_NSET(AE, "L3", lx0)
NSET_NSET(AE, "tip", Tipl)
ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
ORIENTATION(AE, "MaterialOrientation", vec(csmat[:,1]), vec(csmat[:,2]));
SOLID_SECTION(AE, "elasticity", "MaterialOrientation", "AllElements", "Hourglassctl");
SURFACE_SECTION(AE, "TractionElements")
END_INSTANCE(AE);
END_ASSEMBLY(AE);
MATERIAL(AE, "elasticity")
ELASTIC(AE, E1s, E2s, E3s, nu12s, nu13s, nu23s, G12s, G13s, G23s)
SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
STEP_PERTURBATION_STATIC(AE)
BOUNDARY(AE, "ASSEM1.INSTNC1.L1", 1)
BOUNDARY(AE, "ASSEM1.INSTNC1.L2", 2)
BOUNDARY(AE, "ASSEM1.INSTNC1.L3", 3)
DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec([0.0, 0.0, q0]))
NODE_PRINT(AE, "ASSEM1.INSTNC1.tip")
ENERGY_PRINT(AE)
END_STEP(AE)
close(AE)
try rm(AE.filename); catch end

# println("Done: $(  time()-t0 )")
true
end
end
using .mfiber_reinf_cant_yn_strong_Abaqus
mfiber_reinf_cant_yn_strong_Abaqus.test()

module mmorthoballoonpenaltymm
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

    # Orthotropic balloon inflation, axially symmetric model

    # Parameters:
    E1=1.0;
    E2=1.0;
    E3=3.0;
    nu12=0.29;
    nu13=0.29;
    nu23=0.19;
    G12=0.3;
    G13=0.3;
    G23=0.3;
    p= 0.15;
    rin=1.;
    rex =1.2;
    tolerance=rin/1000.

    MR = DeforModelRed2DAxisymm

    fens,fes = Q4block(rex-rin,pi/2,5,20);
    bdryfes = meshboundary(fes);

    icl = selectelem(fens, bdryfes, box=[0.,0.,0.,pi/2],inflate=tolerance);
    for i=1:count(fens)
        r=rin+fens.xyz[i,1]; a=fens.xyz[i,2];
        fens.xyz[i,:]=[r*cos(a) r*sin(a)];
    end

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

    # the symmetry plane
    ly =selectelem(fens, bdryfes; box=[0 rex 0 0], inflate = tolerance)
    # the axis of symmetry
    lx =selectelem(fens, bdryfes; box=[0 0 0 rex], inflate = tolerance)

    # No EBC
    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")

    # The traction boundary condition is applied in the radial
    # direction.

    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(1, 3), true))
    function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
      copyto!(forceout, XYZ/norm(XYZ)*p)
      return forceout
    end
    fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
    F2= distribloads(el1femm, geom, u, fi, 2);

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

    ##
    # The restraints of the nodes on the bounding cross-sections in the direction
    # of the normal to the plane of the cross-section  in the
    # circumferential direction are introduced using a penalty formulation.
    # For that purpose we introduce  a finite element model machine for the
    # surface  finite elements on the cross-sections.
    springcoefficient =1.0e9/ (abs(p)/E1)
    xsfemm = FEMMDeforWinkler(IntegDomain(subset(bdryfes, lx), GaussRule(1, 3), true))
    ysfemm = FEMMDeforWinkler(IntegDomain(subset(bdryfes, ly), GaussRule(1, 3), true))
    H = surfacenormalspringstiffness(xsfemm,  geom, u, springcoefficient, SurfaceNormal(3)) +
        surfacenormalspringstiffness(ysfemm,  geom, u, springcoefficient, SurfaceNormal(3))
    K =stiffness(femm, geom, u)
    U=  (K + H)\(F2)
    scattersysvec!(u,U[:])

    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 3)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
@test abs(minimum(fld.values) - (-0.04635309320688638)) < 1.0e-5
@test abs(maximum(fld.values) - (0.5708149883384825)) < 1.0e-5
    # File =  "orthoballoon_penalty_sigmaz.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)

end
end
using .mmorthoballoonpenaltymm
mmorthoballoonpenaltymm.test()

module mbar1
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    Area = 2.0*phun("in^2")
    E = 30e6*phun("psi") # Young's modulus
    nu = 0.0
    alpha = 5.5e-6*phun("in")/phun("in")/phun("F")
    fens = FENodeSetModule.FENodeSet([10. -5 20;
    30 25 -15]*phun("in") )
    fes = FESetL2(reshape([1,2],1,2))

    # function otherdimensionfu(loc::FFltMat,
    #   conn::CC, N::FFltMat)::FFlt where {CC<:AbstractArray{FInt}}
    #   return otherdimension::FFlt
    # end
    integdomain = IntegDomain(fes, GaussRule(1, 2), (loc, conn, N) -> Area, false)
    # display(IntegDomain)

    MR = DeforModelRed1D
    material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
    # display(material )

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, integdomain, CSys(3, 1), material)
    K = stiffness(femm,  geom,  u)
    # println("K = $(K/(phun("lbf")/phun("in")))")
    ref_K=   1.0e+05 *[
    1.8916    2.8373   -3.3102   -1.8916   -2.8373    3.3102
    2.8373    4.2560   -4.9653   -2.8373   -4.2560    4.9653
   -3.3102   -4.9653    5.7929    3.3102    4.9653   -5.7929
   -1.8916   -2.8373    3.3102    1.8916    2.8373   -3.3102
   -2.8373   -4.2560    4.9653    2.8373    4.2560   -4.9653
    3.3102    4.9653   -5.7929   -3.3102   -4.9653    5.7929];
    @test norm(K/(phun("lbf")/phun("in")) - ref_K)/1.0e5 < 1.0e-3

    dT = NodalField(broadcast(+, zeros(size(fens.xyz, 1), 1),
        100*phun("F"))) # temperature field
    # display(dT)

    F2 = thermalstrainloads(femm, geom, u, dT)
    # println("F2 = $(F2/(phun("lbf")))")
    ref_F=  1.0e+04 *[
    -1.313449091077187
    -1.970173636615779
    2.298535909385076
    1.313449091077187
    1.970173636615779
    -2.298535909385076]
    @test norm(F2/(phun("lbf")) - ref_F) < 1.0e-2

   # K = cholesky(K)
   # U=  K\(F2)
   # scattersysvec!(u, U[:])

    # File = "playground.vtk"
    # MeshExportModule.VTK.vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
    # try rm(File) catch end

end
end
using .mbar1
mbar1.test()

module mbar2
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.FENodeSetModule
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    Area = 1.5
    E = 1.0e7 # Young's modulus
    nu = 0.0
    alpha = 0.0
    fens = FENodeSetModule.FENodeSet([0.0 0;
    0 40;
    40 0;
    40 40;
    80 0;
    80 40] )
    fes = FESetL2([1     3
     1     4
     2     4
     3     4
     3     5
     5     4
     6     4
     5     6])


    MR = DeforModelRed1D
    material = MatDeforElastIso(MR,  0.0, E, nu, alpha)
    # display(material )

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field
    setebc!(u, 1)
    setebc!(u, 2)
    applyebc!(u)
    numberdofs!(u)
    # display(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(1, 1),  (loc, conn, N) -> Area, false), CSys(2, 1), material)
    K = stiffness(femm,  geom,  u)

    fi = ForceIntensity(vec([0 -2000.0]))
    lfemm = FEMMBase(IntegDomain(FESetP1(reshape([3], 1,1)), PointRule()))
    F = distribloads(lfemm,  geom,  u,  fi,  3);
    fi = ForceIntensity(vec([+2000.0 0]))
    lfemm = FEMMBase(IntegDomain(FESetP1(reshape([5], 1,1)), PointRule()))
    F = F + distribloads(lfemm,  geom,  u,  fi,  3);
    fi = ForceIntensity(vec([+4000.0 +6000.0]))
    lfemm = FEMMBase(IntegDomain(FESetP1(reshape([6], 1,1)), PointRule()))
    F = F + distribloads(lfemm,  geom,  u,  fi,  3);

    K = cholesky(K)
    U=  K\F
    scattersysvec!(u, U[:])
    @test norm(u.values  - [ 0.0         0.0
                              0.0         0.0
                              0.0213333   0.0408366
                             -0.016       0.0461699
                              0.0426667   0.150091
                             -0.00533333  0.166091 ]) < 1.0e-4

    sfld =  elemfieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    # display(sfld)
    # println("Cauchy = $(sfld.values)")
    @test norm(sfld.values - [5333.33; 3771.24; -4000.0; 1333.33; 5333.33; -5656.85; 2666.67; 4000.0]) < 1.0e-2
    vfld =  elemfieldfromintegpoints(femm, geom, u, :vm, 1)
    # display(vfld)

    File = "Planar_truss.vtk"
    MeshExportModule.VTK.vtkexportmesh(File, fens, fes;
    scalars=[("sx", sfld.values), ("vm", vfld.values)])
    # @async run(`"paraview.exe" $File`)
    try rm(File) catch end

end
end
using .mbar2
mbar2.test()

module mmmLE10expiAbaqus2mmmm
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
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

    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    qmagn = 1.0*phun("MEGA*PA");# transverse pressure
    sigma_yP = -5.38*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Ae =3.25*phun("m"); # Major radius of the exterior ellipse
    Be =2.75*phun("m"); # Minor radius of the exterior ellipse
    Ai =2.0*phun("m"); # Major radius of the interior ellipse
    Bi =1.0*phun("m"); # Minor radius of the interior ellipse
    Thickness = 0.6*phun("m")
    nc = 6; # number of elements per side
    nr = 5; # number of elements per side
    nt = 2; # number of elements through the thickness
    # nc = 26; # number of elements per side
    # nr = 25; # number of elements per side
    # nt = 18; # number of elements through the thickness
    tolerance = Thickness/nt/1000.; # Geometrical tolerance

    fens,fes = Q4block(1.0, pi/2, nr, nc)
    #
    @test  nt % 2 == 0
    fens,fes  = H8extrudeQ4(fens, fes,
      nt, (xyz, layer)->[xyz[1], xyz[2], (layer)/nt*Thickness]);

    # Select the  boundary faces, on the boundary that is clamped,  and on the part
    # of the boundary that is loaded with the transverse pressure
    bdryfes = meshboundary(fes);
    exteriorbfl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    topbfl = selectelem(fens, bdryfes, box=[0.0, 1.0, 0.0, pi/2, Thickness, Thickness], inflate=tolerance);

    # Reshape the generated block into the elliptical plate
    for i=1:count(fens)
        r=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(r*Ae+(1-r)*Ai)*cos(a) (r*Be+(1-r)*Bi)*sin(a) z];
    end


    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    l12 =connectednodes(subset(bdryfes, exteriorbfl)) # external boundary
    setebc!(u, l12, true, 1, 0.0)
    setebc!(u, l12, true, 2, 0.0)
    ll = selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness/2.0, Thickness/2.0], inflate = tolerance)
    l3 = intersect(ll, connectednodes(subset(bdryfes, exteriorbfl)))
    setebc!(u, l3, true, 3, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0) # symmetry plane X = 0
    l2 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l2,true, 2, 0.0) # symmetry plane Y = 0

    applyebc!(u)
    numberdofs!(u)

    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,topbfl), GaussRule(2, 2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        forceout .=  [0.0, 0.0, -qmagn]
        return forceout
    end
    fi = ForceIntensity(FFlt, 3, pfun);
    F2 = distribloads(el1femm, geom, u, fi, 2);

    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the fully three-dimensional solid.
    MR = DeforModelRed3D

    material = MatDeforElastIso(MR, E, nu)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    # The geometry field now needs to be associated with the FEMM
    femm = associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)
    K = cholesky(K)
    U = K\(F2)
    scattersysvec!(u, U[:])

    nl = selectnode(fens, box=[Ai,Ai,0,0,Thickness,Thickness],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.030939, 0, -0.10488] [MM]")


    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2; nodevalmethod = :averaging, reportat = :extraptrend)#
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yP = $(sigma_yP/phun("MPa")) [MPa]")

    # println("$((nc, nr, nt)), $(fld.values[nl,1][1]/phun("MPa"))")
# println("$(fld.values[nl,1][1]/phun("MPa"))")
    @test abs(fld.values[nl,1][1]/phun("MPa") - -4.627214556813842) < 1.0e-3

    # File =  "LE10NAFEMS_sigmay.vtk"
    # vtkexportmesh(File, fes.conn, geom.values,
    #                FinEtools.MeshExportModule.VTK.H8; vectors=[("u", u.values)],
    #                scalars=[("sigmay", fld.values)])
    # @async run(`"paraview.exe" $File`)
    # true


    AE = AbaqusExporter("LE10NAFEMS_H8");
    HEADING(AE, "LE10NAFEMS: Transverse deflection of elliptical plate with elliptical hole.");
    PART(AE, "part1");
    END_PART(AE);
    ASSEMBLY(AE, "ASSEM1");
    INSTANCE(AE, "INSTNC1", "PART1");
    NODE(AE, fens.xyz);
    ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(femm.integdomain.fes))
    ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(femm.integdomain.fes), connasarray(el1femm.integdomain.fes))
    NSET_NSET(AE, "l1", l1)
    NSET_NSET(AE, "l2", l2)
    NSET_NSET(AE, "l3", l3)
    NSET_NSET(AE, "l12", l12)
    ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
    SURFACE_SECTION(AE, "TractionElements")
    END_INSTANCE(AE);
    END_ASSEMBLY(AE);
    MATERIAL(AE, "elasticity")
    ELASTIC(AE, E, nu)
    SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
    STEP_PERTURBATION_STATIC(AE)
    BOUNDARY(AE, "ASSEM1.INSTNC1", u.is_fixed, u.fixed_values)
    DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec([0.0, 0.0, -qmagn]))
    END_STEP(AE)
    close(AE)
    lines = read(AE.filename)
    @test length(lines) - 10270 == 0
    # try rm(AE.filename) catch end

end
end
using .mmmLE10expiAbaqus2mmmm
mmmLE10expiAbaqus2mmmm.test()

module mplate_w_hole_RECT_MSH8m
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
# using DataFrames
# using CSV
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
    E = 210000*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Ri= 0.1*phun("M"); # hole radius
    Re = 2*Ri; # outer radius
    H = 0.01*phun("M") # thickness of the plate
    nRadial, nCircumferential=6, 3;
    sigma0=1*phun("MEGA*PA");

    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    sigyderrs = Dict{Symbol, FFltVec}()

    nelems = []
    for extrapolation in [:extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nelems = []
        for ref in [1]
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance

            fens,fes = H8elliphole(Ri, Ri, Re, Re, Thickness,
                2^ref*nCircumferential, 2^ref*nCircumferential, 2^ref*nRadial, 1)
            @test count(fes) == 144
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.VTK.H8)
            # @async run(`"paraview.exe" $File`)

            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)

            applyebc!(u)
            numberdofs!(u)


            bdryfes = meshboundary(fes);
            ixl = selectelem(fens, bdryfes, plane=[1.0, 0.0, 0.0, Re], thickness=tolerance);
            elxfemm =  FEMMBase(IntegDomain(subset(bdryfes,ixl), GaussRule(2, 2)))
            function pfunx(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxx(XYZ)
                forceout[2] = sigmaxy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfunx);
            Fx = distribloads(elxfemm, geom, u, fi, 2);
            iyl = selectelem(fens, bdryfes, plane=[0.0, 1.0, 0.0, Re], thickness=tolerance);
            elyfemm =  FEMMBase(IntegDomain(subset(bdryfes,iyl), GaussRule(2, 2)))
            function pfuny(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = -sigmaxy(XYZ)
                forceout[2] = sigmayy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfuny);
            Fy = distribloads(elyfemm, geom, u, fi, 2);

            MR = DeforModelRed3D

            material = MatDeforElastIso(MR, E, nu)

            femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(Fx + Fy)
            scattersysvec!(u, U[:])

            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")

            # println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], dims = 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            # println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            @test abs(sigyA/phun("MPa") - -0.8521990950600441) < 1.0e-3
            sigxB = mean(sigx.values[nlB,1], dims = 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            # println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            @test abs(sigxB/phun("MPa") - 2.7749827820003374) < 1.0e-3
            # println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            # push!(nelems, count(fes))
            # push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.VTK.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA")),
            # ("sigmay", sigy.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end

    # df = DataFrame(nelems=vec(nelems),
    #     sigyderrtrend=vec(sigyderrs[:extraptrend]),
    #     sigyderrdefault=vec(sigyderrs[:extrapmean]))
    # File = "LE1NAFEMS_MSH8_convergence.CSV"
    # CSV.write(File, df)
    # @async run(`"paraview.exe" $File`)

end
end
using .mplate_w_hole_RECT_MSH8m
mplate_w_hole_RECT_MSH8m.test()

module mplate_w_hole_RECT_H20m
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
using FinEtools.MeshImportModule: import_ABAQUS
# using DataFrames
# using CSV
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
    E = 210000*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    Ri= 0.15*phun("M"); # hole radius
    Re = 2*Ri; # outer radius
    H = 0.01*phun("M") # thickness of the plate
    nRadial, nCircumferential=6, 3;
    sigma0=1*phun("MEGA*PA");

    function sigmaxx(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(vec(x[1:2]));
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    sigyderrs = Dict{Symbol, FFltVec}()

    nelems = []
    for extrapolation in [:extrapmean]
        sigyderrs[extrapolation] = FFltVec[]
        nelems = []
        for ref in [1]
            Thickness = H
            # Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance

            fens,fes = H8elliphole(Ri, Ri, Re, Re, Thickness,
            2^ref*nCircumferential, 2^ref*nCircumferential, 2^ref*nRadial, 1)
            fens,fes = H8toH20(fens,fes)
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.VTK.H20)
            # @async run(`"paraview.cexe" $File`)

            # println("My mesh=>$((count(fens), count(fes)))")
            @test count(fens) == 1131
            @test count(fes) == 144
            #
            # output = import_ABAQUS("plane_w_hole_m_debug.inp")
            # fens1,fes1 = output["fens"], output["fesets"][1]
            # println("Matlab mesh=>$((count(fens1), count(fes1[1])))")
            #
            #  fens3, newfes1, fes2 = mergemeshes(fens,fes, fens1,fes1[1], tolerance)
            #  fes3 = cat(2, newfes1)
            #  println("Merged mesh=>$((count(fens3), count(fes3)))")

            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)

            applyebc!(u)
            numberdofs!(u)


            bdryfes = meshboundary(fes);
            # ixl = selectelem(fens, bdryfes, plane=[1.0, 0.0, 0.0, Re], thickness=tolerance);
            ixl = selectelem(fens, bdryfes, box=[Re, Re, -Inf, +Inf, -Inf, +Inf], inflate = tolerance);
            elxfemm =  FEMMBase(IntegDomain(subset(bdryfes,ixl), GaussRule(2, 2)))
            function pfunx(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxx(XYZ)
                forceout[2] = sigmaxy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfunx);
            Fx = distribloads(elxfemm, geom, u, fi, 2);
            # iyl = selectelem(fens, bdryfes, plane=[0.0, 1.0, 0.0, Re], thickness=tolerance);
            iyl = selectelem(fens, bdryfes, box=[-Inf, +Inf, Re, Re, -Inf, +Inf], inflate = tolerance);
            elyfemm =  FEMMBase(IntegDomain(subset(bdryfes,iyl), GaussRule(2, 2)))
            function pfuny(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                forceout[1] = sigmaxy(XYZ)
                forceout[2] = sigmayy(XYZ)
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfuny);
            Fy = distribloads(elyfemm, geom, u, fi, 2);

            MR = DeforModelRed3D

            material = MatDeforElastIso(MR, E, nu)

            femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(Fx + Fy)
            scattersysvec!(u, U[:])
            # println("oof load = $(norm(Fx + Fy, 2))")
            @test abs(norm(Fx + Fy, 2) - 883.437848042617) < 1.0e-2

            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, 00.0], inflate=tolerance);
            pointu = zeros(FFlt,length(nlA),3)
            gathervalues_asmat!(u, pointu, nlA);
            # println("disp@A = $(pointu/phun("mm")) [MM]")
            @test norm(pointu/phun("mm") - [0.00213238 0.0 0.0]) < 1.0e-4
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, 0.0], inflate=tolerance);
            pointu = zeros(FFlt,length(nlB),3)
            gathervalues_asmat!(u, pointu, nlB);
            # println("disp@B = $(pointu/phun("mm")) [MM]")
            @test norm(pointu/phun("mm") - [0.0 -0.000708141 0.0]) < 1.0e-4
            nlC = selectnode(fens, box=[Re, Re, Re, Re, Thickness, Thickness], inflate=tolerance);
            pointu = zeros(FFlt,length(nlC),3)
            gathervalues_asmat!(u, pointu, nlC);
            # println("disp@C = $(pointu/phun("mm")) [MM]")
            @test norm(pointu/phun("mm") - [0.00168556 -0.000455007 -1.4286e-5]) < 1.0e-4

            nlAallz = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlBallz = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlAallz,1], dims = 1)[1]
            sigyAtrue = sigmayy([Ri, 0.0, 0.0])
            # println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            @test abs(sigyA/phun("MPa") - -0.8513053526935438)/(sigyAtrue/phun("MPa")) < 1.0e-4
            sigxB = mean(sigx.values[nlBallz,1], dims = 1)[1]
            sigxBtrue = sigmaxx([0.0, Ri, 0.0])
            # println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            # @test abs(sigxB/phun("MPa") - 2.789413093796375)/3.0 < 1.0e-4
            # println("$extrapolation, $(count(fes)), $(sigyd/phun("MPa"))")
            # push!(nelems, count(fes))
            # push!(sigyderrs[extrapolation], abs(sigyd/sigma_yD - 1.0))
            File =  "a.vtk"
            vtkexportmesh(File, connasarray(fes), geom.values,
                FinEtools.MeshExportModule.VTK.H20; vectors=[("u", u.values)],
                scalars=[("sigmax", sigx.values/phun("MEGA*PA")),
                ("sigmay", sigy.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
            try rm(File); catch end
        end
    end

    # df = DataFrame(nelems=vec(nelems),
    #     sigyderrtrend=vec(sigyderrs[:extraptrend]),
    #     sigyderrdefault=vec(sigyderrs[:extrapmean]))
    # File = "LE1NAFEMS_MSH8_convergence.CSV"
    # CSV.write(File, df)
    # @async run(`"paraview.exe" $File`)

end
end
using .mplate_w_hole_RECT_H20m
mplate_w_hole_RECT_H20m.test()


module mplate_w_hole_MST10m
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential=3, 5;
    sigma0=1*phun("MEGA*PA");

    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    sigxderrs = Dict{Symbol, FFltVec}()
    sigyderrs = Dict{Symbol, FFltVec}()
    numelements = []
    numnodes = []
    for extrapolation in [:extrapmean] # :extraptrend
        sigxderrs[extrapolation] = FFltVec[]
        sigyderrs[extrapolation] = FFltVec[]
        numelements = []
        numnodes = []
        for ref in 1:1
            # Thickness = H
            Thickness = H/2^ref
            tolerance = Thickness/2^ref/1000.; # Geometrical tolerance

            fens,fes = T10block(1.0, pi/2, Thickness, 2^ref*nRadial, 2^ref*nCircumferential, 1)

            bdryfes = meshboundary(fes);
            icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);

            for i=1:count(fens)
                t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
                fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
            end

            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, fens.xyz,
            #     FinEtools.MeshExportModule.VTK.H8)
            # @async run(`"paraview.exe" $File`)

            geom = NodalField(fens.xyz)
            u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

            l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 2, 0.0)
            l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
            setebc!(u,l1,true, 1, 0.0)
            # Plane-stress constraint: assume the plane z=0 is the plane of symmetry of the plate
            l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
            setebc!(u,l1,true, 3, 0.0)
            # If this was enabled, the plane-strain  constraint would be enforced.
            # l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
            # setebc!(u,l1,true, 3, 0.0)

            applyebc!(u)
            numberdofs!(u)
            el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), SimplexRule(2, 3)))
            function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
                local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
                nx = XYZ[1]/r; ny = XYZ[2]/r
                # local sx, sy, txy
                # sx, sy, txy = sigmaxx(XYZ), sigmayy(XYZ), sigmaxy(XYZ)
                # sn = sx * nx^2 + sy * ny^2 + 2 * nx * ny * txy
                # tn = -(sx - sy) * nx * ny + (nx^2 - ny^2) * txy
                # forceout[1] = sn * nx - tn * ny
                # forceout[2] = sn * ny + tn * nx
                # forceout[3] = 0.0
                forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
                forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
                forceout[3] = 0.0
                return forceout
            end
            fi = ForceIntensity(FFlt, 3, pfun);
            F2 = distribloads(el1femm, geom, u, fi, 2);

            MR = DeforModelRed3D

            material = MatDeforElastIso(MR, E, nu)

            femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, SimplexRule(3, 4)), material)

            # The geometry field now needs to be associated with the FEMM
            femm = associategeometry!(femm, geom)

            K = stiffness(femm, geom, u)
            K = cholesky(K)
            U = K\(F2)
            scattersysvec!(u, U[:])

            nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
            nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
            # thecorneru = zeros(FFlt,length(nlA),3)
            # gathervalues_asmat!(u, thecorneru, nl);
            # thecorneru = mean(thecorneru, 1)[1]/phun("mm")
            # println("displacement = $(thecorneru) vs -0.10215 [MM]")

            # println("Extrapolation: $( extrapolation )")
            sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
                nodevalmethod = :averaging, reportat = extrapolation)
            sigyA = mean(sigy.values[nlA,1], dims = 1)[1]
            sigyAtrue = sigmatt([Ri, 0.0, 0.0])
            # println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
            @test  abs(sigyA/phun("MPa") - -0.6705333234697736)/(sigyAtrue/phun("MPa")) < 1.0e-4
            sigxB = mean(sigx.values[nlB,1], dims = 1)[1]
            sigxBtrue = sigmatt([0.0, Ri, 0.0])
            # println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
            @test  abs(sigxB/phun("MPa") - 2.301542874107758)/3.0 < 1.0e-4
            push!(numnodes, count(fens))
            push!(numelements, count(fes))
            push!(sigxderrs[extrapolation], abs(sigxB/sigxBtrue - 1.0))
            push!(sigyderrs[extrapolation], abs(sigyA/sigyAtrue - 1.0))
            # File =  "a.vtk"
            # vtkexportmesh(File, fes.conn, geom.values,
            # FinEtools.MeshExportModule.VTK.H8; vectors=[("u", u.values)],
            # scalars=[("sigmax", sigx.values/phun("MEGA*PA"))])
            # @async run(`"paraview.exe" $File`)
        end
    end

end
end
using .mplate_w_hole_MST10m
mplate_w_hole_MST10m.test()
