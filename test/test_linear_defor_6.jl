
module mmLE1NAFEMSsstress
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Radius = 1.0*phun("m")
    Thickness = 0.1*phun("m")
    n = 2; # number of elements per side
    tolerance = 1.0/n/1000.; # Geometrical tolerance

    fens,fes = Q4block(1.0, pi/2, n, n*2)
    fens,fes  = H8extrudeQ4(fens, fes,
      1, (xyz, layer)->[xyz[1], xyz[2], (layer)*Thickness]);

    bdryfes = meshboundary(fes);
    icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    for i=1:count(fens)
        t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
    end


    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
    setebc!(u,l1,true, 3, 0.0)

    applyebc!(u)
    numberdofs!(u)


    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(2, 2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
        forceout .=    vec(p*pt/norm(pt));
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

    nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, 0.0],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")
    @test norm(thecorneru - [-0.107276 0.0 0.0]) < 1.0e-4


    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :invdistance, reportat = :meanonly)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 42.54884174624546) < 1.0e-3
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :averaging, reportat = :extrapmean)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 42.54884174624546) < 1.0e-3
    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :averaging, reportat = :extraptrend)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 45.44562958746983) < 1.0e-3

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2;
        nodevalmethod = :invdistance, reportat = :meanonly)
    # println("Sigma_y =$(fld.values[nl,1][1]/phun("MPa")) as compared to reference sigma_yD = $(sigma_yD/phun("MPa")) [MPa]")
    @test abs(fld.values[nl,1][1]/phun("MPa") - 42.54884174624546) < 1.0e-3

end
end
using .mmLE1NAFEMSsstress
mmLE1NAFEMSsstress.test()

module mocylpullFun
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Test
function test()

    # Cylinder  pulled by enforced displacement, axially symmetric model

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
    Length = 1*rex
    ua = -0.05*Length
    tolerance=rin/1000.

    ##
    # Note that the FinEtools objects needs to be created with the proper
    # model-dimension reduction at hand.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm


    # Create the mesh and initialize the geometry.  First we are going
    # to construct the block of elements with the first coordinate
    # corresponding to the thickness in the radial direction, and the second
    # coordinate is the thickness in the axial direction.

    fens,fes = Q4block(rex-rin,Length,5,20);
    fens.xyz[:, 1] = fens.xyz[:, 1] .+ rin
    bdryfes = meshboundary(fes);

    # the symmetry plane
    la1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    # The other end
    la2 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)

    e1 = FDataDict("node_list"=>la1, "component"=>2, "displacement"=>x -> 0.0)
    e2 = FDataDict("node_list"=>la2, "component"=>2, "displacement"=>x -> ua)

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

    # Make region
    region = FDataDict("femm"=>femm)

    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region],
    "essential_bcs"=>[e1, e2])

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]

    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050318853446676574) < 1.0e-4
    @test abs(maximum(fld.values) - -0.04973951673608955) < 1.0e-4
    # File =  "orthoballoon_sigmaz.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .mocylpullFun
mocylpullFun.test()

module mmLE11malgo
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
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
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;
    alphaa = 2.3e-4;              # thermal expansion coefficient
    sigmaA = -105*phun("MEGA*Pa")
    nref =  1;                        # how many times should we refine the mesh?
    X = [1.     0.;#A
    1.4    0.;#B
    0.995184726672197   0.098017140329561;
    1.393258617341076 0.137223996461385;
    0.980785  0.195090;#
    1.37309939 0.27312645;
    0.956940335732209   0.290284677254462
    1.339716470025092 0.406398548156247
    0.9238795  0.38268;#C
    1.2124  0.7;#D
    0.7071  0.7071;#E
    1.1062  1.045;#F
    0.7071  (0.7071+1.79)/2;#(E+H)/2
    1.      1.39;#G
    0.7071  1.79;#H
    1.      1.79;#I
    ]*phun("M")
    tolerance  = 1.e-6*phun("M")
    ##
    # Note that the material object needs to be created with the proper
    # model-dimension reduction in mind.  In this case that is the axial symmetry
    # assumption.
    MR  =  DeforModelRed2DAxisymm



    fens = FENodeSet(X);
    fes = FESetQ4([1 2 4 3; 3 4 6 5; 5 6 8 7; 7 8 10 9; 9 10 12 11; 11 12 14 13; 13 14 16 15]);
    for ref = 1:nref
      fens,fes = Q4refine(fens,fes);
      list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
      fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);
    end
    fens,fes = Q4toQ8(fens,fes)
    list = selectnode(fens,distance = 1.0+0.1/2^nref, from = [0. 0.], inflate = tolerance);
    fens.xyz[list,:] =  FinEtools.MeshUtilModule.ontosphere(fens.xyz[list,:],1.0);


    # EBC's
    l1 = selectnode(fens,box = [-Inf Inf 0 0],inflate = tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)
    l1 = selectnode(fens,box = [-Inf Inf 1.79  1.79],inflate = tolerance)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>x -> 0.0)

    # Temperature field
    dtemp = FDataDict("temperature"=>x -> x[1] + x[2])

    # Property and material
    material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

    femm  =  FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 3), true), material)

    # Make region 1
    region = FDataDict("femm"=>femm)
    # Make model data
    modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region],
    "essential_bcs"=>[e1, e2], "temperature_change"=>dtemp)

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]
    dT = modeldata["temp"]

    nA  = selectnode(fens,box = FFlt[1.0  1.0 0.0 0.0], inflate = tolerance);

    fld =  fieldfromintegpoints(femm, geom, u, dT, :Cauchy, 2)


    File  =   "LE11NAFEMS_Q8_sigmay.vtk"
    vtkexportmesh(File, fens, fes; scalars = [("sigmay", fld.values)],
        vectors = [("u", u.values)])
    # println("range of  sigmay = $((minimum(fld.values), maximum(fld.values)))")
    @test norm([minimum(fld.values), maximum(fld.values)] - [-1.443052182185007e8, -1.4106181545272522e7]) < 1.0e-1
        # @async run(`"paraview.exe" $File`)
    try rm(File)   catch end

    sA  =  fld.values[nA]/phun("MEGA*Pa")
    sAn  =  fld.values[nA]/sigmaA
    # println("Stress at point A: $(sA) i. e.  $( sAn*100  )% of reference value")
    @test norm(sA .- -93.8569) < 1.0e-4

    fen2fe  = FENodeToFEMap(connasarray(fes), nnodes(geom))
    function inspector(idat, elnum, conn, xe,  out,  xq)
    #   println("loc = $(  xq  ) : $(  transpose(out)/phun("MEGA*Pa")  )")
      return idat
    end

    inspectintegpoints(femm, geom, u, dT,  fen2fe.map[nA[1]],
      inspector, []; quantity = :Cauchy)

end
end
using .mmLE11malgo
mmLE11malgo.test()

module mmtwistedmsh8ort
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()

  # println("""
  # The initially twisted cantilever beam is one of the standard test
  # problems for verifying the finite-element accuracy [1]. The beam is
  #   clamped at one end and loaded either with unit in-plane or
  #   unit out-of-plane force at the other. The centroidal axis of the beam is
  #   straight at the undeformed  configuration, while its cross-sections are
  #   twisted about the centroidal axis from 0 at the clamped end to pi/2 at
  #   the free end.
  #
  #   Reference:
  #   Zupan D, Saje M (2004) On "A proposed standard set of problems to test
  #   finite element accuracy": the twisted beam. Finite Elements in Analysis
  #   and Design 40: 1445-1451.
  #   """)
  E1s =   E2s =  E3s = 0.29e8
  nu12s = nu13s = nu23s = 0.22
  G12s = G13s = G23s = E1s/2/(1+nu12s)
  # E = 0.29e8;
  # nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 7;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

  # Reshape into a twisted beam shape
  for i = 1:count(fens)
    a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
    fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
  end

  # Clamped end of the beam
  l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
  e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
  e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
  e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

  # Traction on the opposite edge
  boundaryfes  =   meshboundary(fes);
  Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
  el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
  flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
  material))

  # Make model data
  modeldata =  FDataDict(
  "fens"=> fens, "regions"=>  [region1],
  "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

  # Call the solver
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
  geom = modeldata["geom"]
  u = modeldata["u"]

  # Extract the solution
  nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
  theutip = mean(u.values[nl,:], dims = 1)
  # println("displacement  = $(theutip[dir]) as compared to converged $uex")
  @test abs(theutip[dir]-uex)/uex < 0.0012

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of vm, nodal: $([minimum(vm.values),   maximum(vm.values)])")
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
  @test norm([minimum(vm.values),   maximum(vm.values)]-[4.78774, 522.126]) < 0.01

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of vm, elemental: $([minimum(vm.values),   maximum(vm.values)])")
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
  @test norm([minimum(vm.values),   maximum(vm.values)]-[1.85882, 522.126]) < 0.01

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedmsh8ort
mmtwistedmsh8ort.test()

module mmtwistedmsh9ort
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Test
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()

    # println("""
    # The initially twisted cantilever beam is one of the standard test
    # problems for verifying the finite-element accuracy [1]. The beam is
    #   clamped at one end and loaded either with unit in-plane or
    #   unit out-of-plane force at the other. The centroidal axis of the beam is
    #   straight at the undeformed  configuration, while its cross-sections are
    #   twisted about the centroidal axis from 0 at the clamped end to pi/2 at
    #   the free end.
    #
    #   Reference:
    #   Zupan D, Saje M (2004) On "A proposed standard set of problems to test
    #   finite element accuracy": the twisted beam. Finite Elements in Analysis
    #   and Design 40: 1445-1451.
    #   """)
    E1s =   E2s =  E3s = 0.29e8
    nu12s = nu13s = nu23s = 0.22
    G12s = G13s = G23s = E1s/2/(1+nu12s)
    # E = 0.29e8;
    # nu = 0.22;
    W = 1.1;
    L = 12.;
    t =  0.32;
    nl = 2; nt = 1; nw = 1; ref = 7;
    p =   1/W/t;
    #  Loading in the Z direction
    loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
    #   Loading in the Y direction
    #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
    tolerance  = t/1000;

    fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

    # Reshape into a twisted beam shape
    for i = 1:count(fens)
        a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
        fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end

    # Clamped end of the beam
    l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
    e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
    e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
    e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

    # Traction on the opposite edge
    boundaryfes  =   meshboundary(fes);
    Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
    el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
    flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


    # Make the region
    MR = DeforModelRed3D
    material = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
    region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
    material))

    # Make model data
    modeldata =  FDataDict( "fens"=> fens, "regions"=>  [region1],
    "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])

    # Call the solver
    modeldata = AlgoDeforLinearModule.linearstatics(modeldata)
    geom = modeldata["geom"]
    u = modeldata["u"]

    # Extract the solution
    nl = selectnode(fens, box = [L L -100*W 100*W -100*W 100*W],inflate = tolerance);
    theutip = mean(u.values[nl,:], dims = 1)
    # println("displacement  = $(theutip[dir]) as compared to converged $uex")
    @test abs(theutip[dir]-uex)/uex < 0.0012

    # Write out mesh with displacements
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8")
    modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :Cauchy, "component"=> :xy, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

    # Write out mesh with stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :Cauchy, "component"=> :xz, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

    # Write out mesh with von Mises stresses
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8",
    "quantity"=> :vm, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstress(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    # println("extremes of vm, nodal: $([minimum(vm.values),   maximum(vm.values)])")
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
    @test norm([minimum(vm.values),   maximum(vm.values)]-[4.78774, 522.126]) < 0.01

    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
    "quantity"=> :vm, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    vm  = modeldata["postprocessing"]["exported"][1]["field"]
    # println("extremes of vm, elemental: $([minimum(vm.values),   maximum(vm.values)])")
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end
    @test norm([minimum(vm.values),   maximum(vm.values)]-[1.85882, 522.126]) < 0.01

    # Write out mesh with von Mises stresses, elementwise
    modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam_msh8-ew",
    "quantity"=> :Cauchy, "component"=> :xz, "outputcsys"=>CSys(3))
    modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
    try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedmsh9ort
mmtwistedmsh9ort.test()

module mxRMSerror3a1
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshUtilModule
using FinEtools.AlgoBaseModule
using Test
function test()


    elementtag = "MSH8"
    # println("""
    # Pagano_3layer_cylindrical_bending: $(elementtag)
    # """)

    # This example provides three-dimensional finite element model for  the
    # transverse shear stress calculations. The problem consists of a one-, two- or
    # three-layer plate subjected to a sinusoidal distributed load, as
    # described by Pagano (1969). The resulting transverse shear and axial
    # stresses through the thickness of the plate are compared to two existing
    # analytical solutions by Pagano (1969). The first solution is derived from
    # classical laminated plate theory (CPT), while the second is an exact
    # solution from linear elasticity theory.

    filebase = "Pagano_3layer_cylindrical_bending_$(elementtag)_convergence"

    modeldatasequence = FDataDict[]
    for Refinement = [1, 2, 4]

        # Orthotropic material for the 3 layers
        E1 = 25e6*phun("PSI"); E2 = 1e6*phun("PSI"); E3 = E2;
        G12 = 0.5e6*phun("PSI");  G13 = G12; G23 = 0.2e6*phun("PSI")
        nu12 =  0.25; nu13 =  0.25; nu23 =  0.25;
        Span_to_thickness = 4.0;
        T = 2.5*phun("in"); # total thickness of the plate
        L = Span_to_thickness*T;
        h = 1*phun("in");  # depth of the plate
        q0 = 1*phun("PSI")
        CTE1 =  CTE2 =  CTE3 = 0.0

        # Here we define the layout and the thicknesses of the layers.
        angles = vec([0.0 90.0 0.0]);
        nLayers = length(angles)
        ts = T/nLayers * ones(nLayers); # layer thicknesses

        tolerance = 0.0001*T

        # Select how find the mesh should be
        nL, nh = Refinement * 2 * 4, Refinement*1;
        nts= Refinement * 2 * ones(Int, nLayers);# number of elements per layer

        xs = collect(linearspace(0.0, L, nL+1))
        ys = collect(linearspace(0.0, h, nh+1))

        fens,fes = H8layeredplatex(xs, ys, ts, nts)
        # println("count(fens) = $(count(fens))")

        # This is the material  model
        MR = DeforModelRed3D
        skinmaterial = MatDeforElastOrtho(MR,
        0.0, E1, E2, E3,
        nu12, nu13, nu23,
        G12, G13, G23,
        CTE1, CTE2, CTE3)

        # The material coordinate system function is defined as:
        function updatecs!(csmatout::FFltMat, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
            rotmat3!(csmatout, angles[fe_label]/180.0*pi* [0.0; 0.0; 1.0]);
        end

        # The vvolume integrals are evaluated using this rule
        gr = GaussRule(3, 2)

        # We will create 3 regions, one for each of the layers
        regions = FDataDict[]
        for layer = 1:nLayers
            rls = selectelem(fens, fes, label =  layer)
            push!(regions, FDataDict("femm"=>FEMMDeforLinearMSH8(MR,
            IntegDomain(subset(fes, rls), gr), CSys(3, 3, updatecs!), skinmaterial)))
        end

        # File =  "Meyer_Piening_sandwich-r1.vtk"
        # vtkexportmesh(File, skinregion["femm"].integdomain.fes.conn, fens.xyz, FinEtools.MeshExportModule.VTK.VTK.H8)
        # # @async run(`"paraview.exe" $File`)


        # The essential boundary conditions are applied to enforce the plane strain constraint.
        ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
        lyh = selectnode(fens, box=[-Inf Inf h h -Inf Inf], inflate=tolerance)
        ey = FDataDict("displacement"=>  0.0, "component"=> 2, "node_list"=>vcat(ly0, lyh))
        # The transverse displacement is fixed at the two ends.
        lz0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
        lzL = selectnode(fens, box=[L L -Inf Inf -Inf Inf], inflate=tolerance)
        ez = FDataDict("displacement"=>  0.0, "component"=> 3, "node_list"=>vcat(lz0, lzL))
        ex = FDataDict("displacement"=>  0.0, "component"=> 1, "node_list"=>[1])

        # The traction boundary condition is applied at the top of the plate.
        bfes = meshboundary(fes)
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            forceout[1] = 0.0
            forceout[2] = 0.0
            forceout[3] = -q0*sin(pi*XYZ[1]/L)
            return forceout
        end
        # From  the entire boundary we select those quadrilaterals that lie on the plane
        # Z = thickness
        tl = selectelem(fens, bfes, box = [-Inf Inf -Inf Inf T T], inflate=tolerance)
        Trac = FDataDict("traction_vector"=>pfun,
        "femm"=>FEMMBase(IntegDomain(subset(bfes, tl), GaussRule(2, 2))))

        modeldata = FDataDict("fens"=>fens,
        "regions"=>regions,
        "essential_bcs"=>[ex, ey, ez],
        "traction_bcs"=> [Trac]
        )
        modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

        modeldata["postprocessing"] = FDataDict("file"=>filebase * "-u")
        modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
        for e in modeldata["postprocessing"]["exported"]
            try rm(e["file"]) catch end
        end

        u = modeldata["u"]
        geom = modeldata["geom"]

        # The results of the displacement and stresses will be reported at
        # nodes located at the appropriate points.
        ntopcenter = selectnode(fens, box=[L/2 L/2 0.0 h T T], inflate=tolerance)
        ncenterline = selectnode(fens, box=[L/2 L/2 0.0 0.0 0.0 T], inflate=tolerance)
        nx0line = selectnode(fens, box=[0.0 0.0 0.0 0.0 0.0 T], inflate=tolerance)

        zclo = sortperm(vec(geom.values[ncenterline, 3]))
        ncenterline = ncenterline[zclo]
        centerz = geom.values[ncenterline, 3]

        # println("Top Center deflection: $(mean(u.values[ntopcenter, 3], 1)/phun("in")) [in]")

        # # extrap = :extrapmean
        extrap = :extraptrend
        nodevalmeth = :averaging
        # extrap = :default
        # nodevalmeth = :invdistance

        # Compute  all stresses
        modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s",
        "quantity"=>:Cauchy, "component"=>collect(1:6), "outputcsys"=>CSys(3),
        "nodevalmethod"=>nodevalmeth, "reportat"=>extrap)
        modeldata = AlgoDeforLinearModule.exportstress(modeldata)
        for e in modeldata["postprocessing"]["exported"]
            try rm(e["file"]) catch end
        end

        modeldata["elementsize"] = 1.0/Refinement
        modeldata["geometricaltolerance"] = tolerance
        modeldata["targetfields"] = [e["field"] for e in modeldata["postprocessing"]["exported"]]
        push!(modeldatasequence, modeldata)
    end # for refinement

    elementsizes, errornorms, p = AlgoBaseModule.evalconvergencestudy(modeldatasequence)

    # println("")
    # println("Normalized Approximate Error = $(errornorms)")
    @test abs(p[1] - 1.3347513854727369)/1.3347513854727369 < 1.0e-3
    # csvFile = filebase * "_errors" * ".CSV"
    # savecsv(csvFile,
    # elementsizes=vec(elementsizes[1:end-1]),
    # elementsizes2=vec(elementsizes[1:end-1].^2),
    # elementsizes3=vec(elementsizes[1:end-1].^3),
    # errornorms=vec(errornorms)
    # )

    # @async run(`"paraview.exe" $csvFile`)

    # println("Done")

end
end
using .mxRMSerror3a1
mxRMSerror3a1.test()

module munit_cube_modes_nice_t4
using FinEtools
using FinEtoolsDeforLinear
using Test
using Arpack
using LinearAlgebra
function test()
    # println("""
    # Vibration modes of unit cube  of almost incompressible material.
    # %
    # Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
    # tetrahedral. International Journal for Numerical Methods in
    # Engineering 67: 841-867.
    # """)
    t0 = time()

    E = 1*phun("PA");
    nu = 0.499;
    rho = 1*phun("KG/M^3");
    a = 1*phun("M"); b = a; h =  a;
    n1 = 10;# How many element edges per side?
    na =  n1; nb =  n1; nh  = n1;
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (0.01*2*pi)^2;
    stabfact = 0.015
    Eigenvalues = [0.0, 5.93656e-8, 7.54751e-8, 9.80131e-8, 1.14899e-7, 1.27725e-7, 0.264544, 0.266128, 0.350568, 0.352546, 0.355279, 0.357389, 0.357701, 0.359704, 0.402389, 0.402968, 0.404977, 0.45061, 0.450974, 0.452039]

    MR = DeforModelRed3D
    fens,fes  = T4block(a,b,h, na,nb,nh)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material, stabfact)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # # println("Eigenvalues: $fs [Hz]")
    @test norm(vec(fs) .- vec(Eigenvalues)) < 1.0e-4*maximum(vec(Eigenvalues))

    # mode = 17
    # scattersysvec!(u, v[:,mode])
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])

    true

end
end
using .munit_cube_modes_nice_t4
munit_cube_modes_nice_t4.test()

module malum_cyl_mode_nice_t4
using FinEtools
using FinEtoolsDeforLinear
using Test
using Arpack
using LinearAlgebra
function test()
    # Aluminum cylinder free vibration, mesh imported from Abaqus
    # Mesh converted from quadratic tetrahedra to linear tetrahedra
    # NICE tetrahedral elements used
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    stabfact = 0.005
    Eigenvalues = [4.54746e-5, 6.82231e-5, 8.7071e-5, 9.99708e-5, 0.000112778, 0.000116397, 2533.6, 2535.12, 2574.64, 4086.61, 4652.66, 4654.16, 5122.94, 6755.62, 6756.45, 6872.26, 6875.3, 6883.49, 6888.53, 6983.99]

    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material, stabfact)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # # println("Eigenvalues: $fs [Hz]")
    @test norm(vec(fs) .- vec(Eigenvalues)) < 1.0e-3*maximum(vec(Eigenvalues))

    true
end
end
using .malum_cyl_mode_nice_t4
malum_cyl_mode_nice_t4.test()

module malum_cyl_mode_esnice_t4
using FinEtools
using FinEtoolsDeforLinear
using Test
using Arpack
using LinearAlgebra
function test()
    # Aluminum cylinder free vibration, mesh imported from Abaqus
    # Mesh converted from quadratic tetrahedra to linear tetrahedra
    # NICE tetrahedral elements used
    E = 70000*phun("MPa");
    nu = 0.33;
    rho = 2700*phun("KG/M^3");
    radius = 0.5*phun("ft");
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0*2*pi)^2;
    stabfact = 0.005
    Eigenvalues =   [0.0, 0.0, 0.0, 1.8846e-5, 7.35917e-5, 0.000119445, 2498.15, 2498.88, 2513.31, 4082.65, 4585.99, 4586.42, 4987.01, 6648.02, 6648.48, 6679.04, 6682.16, 6777.89, 6780.59, 6799.36]
    MR = DeforModelRed3D
    output = import_ABAQUS("alum_cyl.inp")
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K  = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
    d = d .- OmegaShift;
    fs = real(sqrt.(complex(d)))/(2*pi)
    # println("Eigenvalues: $fs [Hz]")
    @test norm(vec(fs) .- vec(Eigenvalues)) < 1.0e-3*maximum(vec(Eigenvalues))

    true
end
end
using .malum_cyl_mode_esnice_t4
malum_cyl_mode_esnice_t4.test()


module mocylpull14 # From linear deformation
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    # Cylinder  compressed by enforced displacement, axially symmetric model


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
    Length = 1*rex
    ua = -0.05*Length
    tolerance=rin/1000.

    ##
    # Note that the FinEtools objects needs to be created with the proper
    # model-dimension reduction at hand.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm

    fens,fes = Q4block(rex-rin,Length,5,20);
    fens.xyz[:, 1] = fens.xyz[:, 1] .+ rin
    bdryfes = meshboundary(fes);

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # The other end
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, ua)

    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")
    @test u.nfreedofs == 240

    # Property and material
    material=MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    # display(material)
    # println("$(material.D)")

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])

    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 2)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - 0.0) < 1.0e-5
    @test abs(maximum(fld.values) - 0.0) < 1.0e-5
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 3)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050) < 1.0e-5
    @test abs(maximum(fld.values) - -0.04999999999999919) < 1.0e-5

    # File =  "mocylpull14.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .mocylpull14
mocylpull14.test()


module mocylpull1 # From deformation
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    # Cylinder  pulled by enforced displacement, axially symmetric model


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
    Length = 1*rex
    ua = -0.05*Length
    tolerance=rin/1000.

    ##
    # Note that the FinEtools objects needs to be created with the proper
    # model-dimension reduction at hand.  In this case that is the axial symmetry
    # assumption.
    MR = DeforModelRed2DAxisymm


    # Create the mesh and initialize the geometry.  First we are going
    # to construct the block of elements with the first coordinate
    # corresponding to the thickness in the radial direction, and the second
    # coordinate is the thickness in the axial direction.

    fens,fes = Q4block(rex-rin,Length,5,20);
    fens.xyz[:, 1] = fens.xyz[:, 1] .+ rin
    bdryfes = meshboundary(fes);

    # now we create the geometry and displacement fields
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

    # the symmetry plane
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    # The other end
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, ua)

    applyebc!(u)
    numberdofs!(u)
    # println("Number of degrees of freedom = $(u.nfreedofs)")
    @test u.nfreedofs == 240

    # Property and material
    material=MatDeforElastOrtho(MR, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

    K =stiffness(femm, geom, u)
    F = nzebcloadsstiffness(femm, geom, u)
    U=  K\(F)
    scattersysvec!(u,U[:])

    # Produce a plot of the radial stress component in the cylindrical
    # coordinate system. Note that this is the usual representation of
    # stress using nodal stress field.

    fld= fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
    # println("Minimum/maximum = $(minimum(fld.values))/$(maximum(fld.values))")
    @test abs(minimum(fld.values) - -0.050318853446676436) < 1.0e-5
    @test abs(maximum(fld.values) - -0.0497395167360893) < 1.0e-5
    # File =  "orthoballoon_sigmaz.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmaz", fld.values)],
    #               vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
end
end
using .mocylpull1
mocylpull1.test()


module mmLE1NAFEMSsstressx1
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 210e3*phun("MEGA*PA");# 210e3 MPa
    nu = 0.3;
    p = 10*phun("MEGA*PA");# 10 MPA Outward pressure on the outside ellipse
    sigma_yD = 92.7*phun("MEGA*PA");# tensile stress at [2.0, 0.0] meters
    Radius = 1.0*phun("m")
    Thickness = 0.1*phun("m")
    n = 20; # number of elements per side
    tolerance = 1.0/n/1000.; # Geometrical tolerance

    fens,fes = Q4block(1.0, pi/2, n, n*2)
    fens,fes  = H8extrudeQ4(fens, fes,
      1, (xyz, layer)->[xyz[1], xyz[2], (layer)*Thickness]);

    bdryfes = meshboundary(fes);
    icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);
    for i=1:count(fens)
        t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
        fens.xyz[i,:]=[(t*3.25+(1-t)*2)*cos(a), (t*2.75+(1-t)*1)*sin(a), z];
    end


    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
    setebc!(u,l1,true, 1, 0.0)
    l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
    setebc!(u,l1,true, 3, 0.0)

    applyebc!(u)
    numberdofs!(u)


    el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(2, 2)))
    function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
        pt= [2.75/3.25*XYZ[1], 3.25/2.75*XYZ[2], 0.0]
        forceout .=    vec(p*pt/norm(pt));
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

    nl = selectnode(fens, box=[2.0, 2.0, 0.0, 0.0, 0.0, 0.0],inflate=tolerance);
    thecorneru = zeros(FFlt,1,3)
    gathervalues_asmat!(u, thecorneru, nl);
    thecorneru = thecorneru/phun("mm")
    # println("displacement =$(thecorneru) [MM] as compared to reference [-0.10215,0] [MM]")
    @test norm(thecorneru - [-0.10082864119023721 0.0 0.0]) < 1.0e-4


    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    @test norm(maximum(fld.values) - 8.190847372888073e7)/8.190847372888073e7 <= 1.0e-4
    # File =  "mmLE1NAFEMSsstressx1-s1.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigma_1", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    fld= fieldfromintegpoints(femm, geom, u, :maxshear, 1)
    @test norm(maximum(fld.values) - 8.19559427861603e7)/8.190847372888073e7 <= 1.0e-4
    # File =  "mmLE1NAFEMSsstressx1-maxshear.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("maxshear", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    
end
end
using .mmLE1NAFEMSsstressx1
mmLE1NAFEMSsstressx1.test()

module mholestr1
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 100.0;
    nu = 1.0/3;
    cte = 0.15
    xradius = 1.0
    yradius = 1.0
    L = 3.0
    H = 3.0
    # nL = 50
    # nH = 50
    # nW = 70
    nL = 100
    nH = 100
    nW = 120
    tolerance = min(xradius, yradius, L, H)/min(nL, nH, nW)/1000.;#Geometrical tolerance

    fens, fes = Q4elliphole(xradius, yradius, L, H,    nL, nH, nW)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
    setebc!(u, l1, 1, 0.0)
    l2 = selectnode(fens; box=[-Inf,  Inf, 0, 0],  inflate = tolerance)
    setebc!(u, l2, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)

    boundaryfes =  meshboundary(fes);
    lce = selectelem(fens, boundaryfes, facing=true, direction = x -> -x,  inflate=  tolerance);
    lc = connectednodes(subset(boundaryfes, lce))
    
    le = selectelem(fens, boundaryfes,  box= [-Inf,  Inf, H, H],  inflate=  tolerance);
    el1femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([0.0, 1.0]);
    Fm = distribloads(el1femm,  geom,  u,  fi,  2);
    le = selectelem(fens, boundaryfes,  box= [L, L, -Inf,  Inf],  inflate=  tolerance);
    el2femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([1.0, 0.0]);
    Fm += distribloads(el2femm,  geom,  u,  fi,  2);

    MR = DeforModelRed2DStress
    material = MatDeforElastIso(MR,  0.0, E, nu, cte)

    femm = FEMMDeforLinear(MR, IntegDomain(fes,  GaussRule(2, 2)),  material)
    
    K = stiffness(femm,  geom,  u)
    U=  K\(Fm)
    scattersysvec!(u, U[:])

    
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # @test norm(maximum(fld.values) - 8.190847372888073e7)/8.190847372888073e7 <= 1.0e-4
    # File =  "mholestr1-s1.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigma_1", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    fld= fieldfromintegpoints(femm, geom, u, :maxshear, 1)
    @test norm(maximum(fld.values) - 2.210557410276065)/2.210557410276065 <= 1.0e-4
    # File =  "mholestr1-maxshear.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("maxshear", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    

end
end
using .mholestr1
mholestr1.test()


module mholestr2
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 100.0;
    nu = 1.0/3;
    cte = 0.15
    xradius = 1.0
    yradius = 1.0
    L = 3.0
    H = 3.0
    # nL = 50
    # nH = 50
    # nW = 70
    nL = 100
    nH = 100
    nW = 120
    tolerance = min(xradius, yradius, L, H)/min(nL, nH, nW)/1000.;#Geometrical tolerance

    fens, fes = Q4elliphole(xradius, yradius, L, H,    nL, nH, nW)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
    setebc!(u, l1, 1, 0.0)
    l2 = selectnode(fens; box=[-Inf,  Inf, 0, 0],  inflate = tolerance)
    setebc!(u, l2, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)

    boundaryfes =  meshboundary(fes);
    lce = selectelem(fens, boundaryfes, facing=true, direction = x -> -x,  inflate=  tolerance);
    lc = connectednodes(subset(boundaryfes, lce))
    
    le = selectelem(fens, boundaryfes,  box= [-Inf,  Inf, H, H],  inflate=  tolerance);
    el1femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([0.0, 1.0]);
    Fm = distribloads(el1femm,  geom,  u,  fi,  2);
    le = selectelem(fens, boundaryfes,  box= [L, L, -Inf,  Inf],  inflate=  tolerance);
    el2femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([1.0, 0.0]);
    Fm += distribloads(el2femm,  geom,  u,  fi,  2);

    MR = DeforModelRed2DStrain
    material = MatDeforElastIso(MR,  0.0, E, nu, cte)

    femm = FEMMDeforLinear(MR, IntegDomain(fes,  GaussRule(2, 2)),  material)
    
    K = stiffness(femm,  geom,  u)
    U=  K\(Fm)
    scattersysvec!(u, U[:])

    
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # @test norm(maximum(fld.values) - 8.190847372888073e7)/8.190847372888073e7 <= 1.0e-4
    # File =  "mholestr2-s1.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigma_1", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    fld= fieldfromintegpoints(femm, geom, u, :maxshear, 1)
    @test norm(maximum(fld.values) - 2.2102738887214257)/2.2102738887214257 <= 1.0e-4
    # File =  "mholestr2-maxshear.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("maxshear", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    

end
end
using .mholestr2
mholestr2.test()

module mholestr3
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()
    E = 100.0;
    nu = 1.0/3;
    cte = 0.15
    xradius = 1.0
    yradius = 1.0
    L = 3.0
    H = 3.0
    # nL = 50
    # nH = 50
    # nW = 70
    nL = 100
    nH = 100
    nW = 120
    tolerance = min(xradius, yradius, L, H)/min(nL, nH, nW)/1000.;#Geometrical tolerance

    fens, fes = Q4elliphole(xradius, yradius, L, H,    nL, nH, nW)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 2)) # displacement field

    l1 = selectnode(fens; box=[0, 0, -Inf,  Inf],  inflate = tolerance)
    setebc!(u, l1, 1, 0.0)
    l2 = selectnode(fens; box=[-Inf,  Inf, 0, 0],  inflate = tolerance)
    setebc!(u, l2, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)

    boundaryfes =  meshboundary(fes);
    lce = selectelem(fens, boundaryfes, facing=true, direction = x -> -x,  inflate=  tolerance);
    lc = connectednodes(subset(boundaryfes, lce))
    
    le = selectelem(fens, boundaryfes,  box= [-Inf,  Inf, H, H],  inflate=  tolerance);
    el1femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([0.0, -1.0]);
    Fm = distribloads(el1femm,  geom,  u,  fi,  2);
    le = selectelem(fens, boundaryfes,  box= [L, L, -Inf,  Inf],  inflate=  tolerance);
    el2femm =  FEMMBase(IntegDomain(subset(boundaryfes, le),  GaussRule(1, 2)))
    fi = ForceIntensity([1.0, 0.0]);
    Fm += distribloads(el2femm,  geom,  u,  fi,  2);

    MR = DeforModelRed2DStress
    material = MatDeforElastIso(MR,  0.0, E, nu, cte)

    femm = FEMMDeforLinear(MR, IntegDomain(fes,  GaussRule(2, 2)),  material)
    
    K = stiffness(femm,  geom,  u)
    U=  K\(Fm)
    scattersysvec!(u, U[:])

    
    fld= fieldfromintegpoints(femm, geom, u, :princCauchy, 1)
    # @test norm(maximum(fld.values) - 8.190847372888073e7)/8.190847372888073e7 <= 1.0e-4
    File =  "mholestr3-s1.vtk"
    vtkexportmesh(File, fens, fes; scalars=[("sigma_1", fld.values)], vectors=[("u", u.values)])
    @async run(`"paraview.exe" $File`)
    fld= fieldfromintegpoints(femm, geom, u, :maxshear, 1)
    @test norm(maximum(fld.values) - 5.921999943843146)/5.921999943843146 <= 1.0e-4
    # File =  "mholestr3-maxshear.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("maxshear", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    

end
end
using .mholestr3
mholestr3.test()



