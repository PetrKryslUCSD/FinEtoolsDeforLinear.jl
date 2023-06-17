
module mmmNAFEMS_R0031_3m
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using Test
import Statistics: mean
function test()
  # println("""
  # NAFEMS publication R0031/3 Composite plate test.
  # Simply supported on all four edges.  Uniform transverse  loading.
  # The modeled part is one quarter of the full plate here.
  # """)

  # This is a test recommended by the National Agency for Finite Element Methods
  # and Standards (U.K.): Test R0031/3 from NAFEMS publication R0031, “Composites
  # Benchmarks,” February 1995.
  t0 = time()
  # Skin (face) material parameters
  E1s = 1.0e7*phun("psi")
  E2s = 0.4e7*phun("psi")
  E3s = 0.4e7*phun("psi")
  nu12s = 0.3
  nu13s = 0.3
  nu23s = 0.3
  G12s = 0.1875e7*phun("psi")
  G13s = 0.1875e7*phun("psi")
  G23s = 0.1875e7*phun("psi")
  # Core material parameters
  E1c = 10.0*phun("psi")
  E2c = 10.0*phun("psi")
  E3c = 10e4.*phun("psi")
  nu12c = 0.
  nu13c = 0.
  nu23c = 0.
  G12c = 10.0*phun("psi")
  G13c = 3.0e4*phun("psi")
  G23c = 1.2e4*phun("psi")
  L = 10.0*phun("in") # side of the square plate
  nL = 8 # number of elements along the side of the plate
  tolerance = 0.0001*phun("in")
  xs = collect(linearspace(0.0, L/2, nL+1))
  ys = collect(linearspace(0.0, L/2, nL+1))
  ts = [0.028; 0.75; 0.028]*phun("in")
  nts = [2; 3;  2; ] # number of elements through the thickness
  tmag = 100*phun("psi")

  # Generate mesh
  fens,fes = H8layeredplatex(xs, ys, ts, nts)
  fens,fes = H8toH20(fens,fes)

  MR = DeforModelRed3D
  skinmaterial = MatDeforElastOrtho(MR,
    0.0, E1s, E2s, E3s,
    nu12s, nu13s, nu23s,
    G12s, G13s, G23s,
    0.0, 0.0, 0.0)
  corematerial = MatDeforElastOrtho(MR,
    0.0, E1c, E2c, E3c,
    nu12c, nu13c, nu23c,
    G12c, G13c, G23c,
    0.0, 0.0, 0.0)

  gr = GaussRule(3, 3)

  rl1 = selectelem(fens, fes, label=1)
  skinbot = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(subset(fes, rl1), gr), skinmaterial))

  rl3 = selectelem(fens, fes, label=3)
  skintop = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(subset(fes, rl3), gr), skinmaterial))

  rl2 = selectelem(fens, fes, label=2)
  core = FDataDict("femm"=>FEMMDeforLinear(MR,
      IntegDomain(subset(fes, rl2), gr), corematerial))

  lx0 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
  lxL2 = selectnode(fens, box=[L/2 L/2 -Inf Inf -Inf Inf], inflate=tolerance)
  ly0 = selectnode(fens, box=[-Inf Inf 0.0 0.0 -Inf Inf], inflate=tolerance)
  lyL2 = selectnode(fens, box=[-Inf Inf L/2 L/2 -Inf Inf], inflate=tolerance)

  ex0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>lx0 )
  exL2 = FDataDict( "displacement"=>  0.0, "component"=> 1, "node_list"=>lxL2 )
  ey0 = FDataDict( "displacement"=>  0.0, "component"=> 3, "node_list"=>ly0 )
  eyL2 = FDataDict( "displacement"=>  0.0, "component"=> 2, "node_list"=>lyL2 )

  bfes = meshboundary(fes)
  ttopl = selectelem(fens, bfes; facing=true, direction = [0.0 0.0 1.0])
  Trac = FDataDict("traction_vector"=>[0.0; 0.0; -tmag],
      "femm"=>FEMMBase(IntegDomain(subset(bfes, ttopl), GaussRule(2, 3))))

  modeldata = FDataDict("fens"=>fens,
   "regions"=>[skinbot, core, skintop], "essential_bcs"=>[ex0, exL2, ey0, eyL2],
   "traction_bcs"=> [Trac]
   )
  modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

  u = modeldata["u"]
  geom = modeldata["geom"]
  lcenter = selectnode(fens, box=[L/2 L/2  L/2 L/2 -Inf Inf], inflate=tolerance)
  cdis = mean(u.values[lcenter, 3])/phun("in")
  # println("Center node displacements $(cdis) [in]; NAFEMS-R0031-3 lists –0.123	[in]")
  # println("")
  @test abs(cdis - (-0.13634800328800462)) < 1.0e-5

  File =  "NAFEMS-R0031-3-plate.vtk"
  vtkexportmesh(File, connasarray(fes), geom.values, FinEtools.MeshExportModule.VTK.H20;
      scalars = [("Layer", fes.label)], vectors = [("displacement", u.values)])
  # @async run(`"paraview.exe" $File`)
  try  rm(File); catch end

end
end
using .mmmNAFEMS_R0031_3m
mmmNAFEMS_R0031_3m.test()

module mmtwistedmsh8mmm
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
import LinearAlgebra: norm, cholesky, cross
using Test
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
  E = 0.29e8;
  nu = 0.22;
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
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
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
using .mmtwistedmsh8mmm
mmtwistedmsh8mmm.test()

module mmunitmmccubemmvibrationmmms
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
function test()

  # println("""
  # % Vibration modes of unit cube  of almost incompressible material.
  # % Mean-strain hexahedron.
  # % Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
  # % tetrahedral. International Journal for Numerical Methods in
  # % Engineering 67: 841-867.""")
  t0 = time()


  E = 1*phun("PA");
  nu = 0.499;
  rho= 1*phun("KG/M^3");
  a=1*phun("M"); b=a; h= a;
  n1=8 # How many element edges per side?
  na= n1; nb= n1; nh =n1;
  neigvs=20                   # how many eigenvalues
  omega_shift=(0.1*2*pi)^2;

  fens,fes = H8block(a,b,h, na,nb,nh)

  # Make the region
  MR = DeforModelRed3D
  material = MatDeforElastIso(MR, rho, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
    material), "femm_mass"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,3)),
    material))

  # Make model data
  modeldata =  FDataDict(
    "fens"=> fens, "regions"=>  [region1],
    "omega_shift"=>omega_shift, "neigvs"=>neigvs)

  # Solve
  modeldata = FinEtoolsDeforLinear.AlgoDeforLinearModule.modal(modeldata)

  fs = modeldata["omega"]/(2*pi)
  # println("Eigenvalues: $fs [Hz]")
  @test norm(fs-[1.92866e-7, 2.07497e-7, 2.16105e-7, 2.31656e-7, 2.35711e-7, 2.53067e-7, 0.266016, 0.266016, 0.364001, 0.364001, 0.364001, 0.366888, 0.366888, 0.366888, 0.415044, 0.415044, 0.41703, 0.467364, 0.467364, 0.467364]) < 0.0001

  modeldata["postprocessing"] = FDataDict("file"=>"unit_cube_mode",
    "mode"=>10)
  modeldata=FinEtoolsDeforLinear.AlgoDeforLinearModule.exportmode(modeldata)
  # @async run(`"paraview.exe" $(modeldata["postprocessing"]["file"]*"1.vtk")`)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  true

end
end
using .mmunitmmccubemmvibrationmmms
mmunitmmccubemmvibrationmmms.test()

module mmtwistedbeamisomm
using FinEtools
using FinEtoolsDeforLinear
using Test
using FinEtoolsDeforLinear.AlgoDeforLinearModule
import LinearAlgebra: norm, cholesky, cross
import Statistics: mean
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 4;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)

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
  material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
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
    # println("normalized displacement  = $(theutip[dir]/uex*100) %")
    @test abs(theutip[dir]/uex*100-99.85504856450584) < 1.0e-6

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
  @test norm([minimum(vm.values),   maximum(vm.values)] - [6.94796, 451.904]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-1-ew",
  "quantity"=> :princCauchy, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of first principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [0.493918, 459.106]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-3-ew",
  "quantity"=> :princCauchy, "component"=> 3)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of third principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-459.106, -0.493918]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end


  # Write out mesh with pressure, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-press-ew",
  "quantity"=> :pressure, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of pressure: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-160.396, 160.396]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedbeamisomm
mmtwistedbeamisomm.test()

module mmtwistedbeamoorthomm
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm, cholesky, cross
using FinEtoolsDeforLinear.AlgoDeforLinearModule
import Statistics: mean
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 4;
  p =   1/W/t;
  #  Loading in the Z direction
  loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
  #   Loading in the Y direction
  #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
  tolerance  = t/1000;

  fens,fes  = H20block(L,W,t, nl*ref,nw*ref,nt*ref)

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
  material = MatDeforElastOrtho(MR, E, nu)
  region1 = FDataDict("femm"=>FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)),
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
    # println("normalized displacement  = $(theutip[dir]/uex*100) %")
    @test abs(theutip[dir]/uex*100-99.85504856450584) < 1.0e-6

  # Write out mesh with displacements
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
  modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xy)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstress(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :vm)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  vm  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
  @test norm([minimum(vm.values),   maximum(vm.values)] - [6.94796, 451.904]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with von Mises stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
  "quantity"=> :Cauchy, "component"=> :xz)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-1-ew",
  "quantity"=> :princCauchy, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of first principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [0.493918, 459.106]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

  # Write out mesh with principal stresses, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-3-ew",
  "quantity"=> :princCauchy, "component"=> 3)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of third principal stress: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-459.106, -0.493918]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end


  # Write out mesh with pressure, elementwise
  modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-press-ew",
  "quantity"=> :pressure, "component"=> 1)
  modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
  ps  = modeldata["postprocessing"]["exported"][1]["field"]
  # println("extremes of pressure: $([minimum(ps.values),   maximum(ps.values)])")
  @test norm([minimum(ps.values),   maximum(ps.values)] - [-160.396, 160.396]) < 1.e-3
  try rm(modeldata["postprocessing"]["exported"][1]["file"]); catch end

end
end
using .mmtwistedbeamoorthomm
mmtwistedbeamoorthomm.test()

module muunit_cube_modes_exportmmm
using FinEtools
using FinEtoolsDeforLinear
using FinEtools.MeshExportModule
using Test
import Arpack: eigs
import LinearAlgebra: norm, cholesky, cross
function test()


  # println("""
  # Vibration modes of unit cube  of almost incompressible material.
  #
  # This example EXPORTS the model to Abaqus.
  #
  # Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
  # tetrahedral. International Journal for Numerical Methods in
  # Engineering 67: 841-867.
  # """)
  t0 = time()


  E = 1*phun("PA");
  nu = 0.499;
  rho = 1*phun("KG/M^3");
  a = 1*phun("M"); b = a; h =  a;
  n1 = 5;# How many element edges per side?
  na =  n1; nb =  n1; nh  = n1;
  neigvs = 20                   # how many eigenvalues
  OmegaShift = (0.01*2*pi)^2;

  MR = DeforModelRed3D
  fens,fes  = H20block(a,b,h, na,nb,nh)

  geom = NodalField(fens.xyz)
  u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

  numberdofs!(u)

  material=MatDeforElastIso(MR, rho, E, nu, 0.0)

  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,2)), material)

  K =stiffness(femm, geom, u)
  femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3,3)), material)
  M =mass(femm, geom, u)
  d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM, explicittransform=:none)
  d = broadcast(-, d, OmegaShift);
  fs = real(sqrt.(complex(d)))/(2*pi)
  # println("Eigenvalues: $fs [Hz]")
  @test norm(fs - [2.73674e-7, 3.00469e-7, 3.14245e-7, 3.19946e-7, 3.42634e-7, 3.56347e-7, 0.262723, 0.262723, 0.357791, 0.357791, 0.357791, 0.36088, 0.36088, 0.36088, 0.408199, 0.408397, 0.408397, 0.461756, 0.461756, 0.461756]) < 1.0e-3

  mode = 7
  scattersysvec!(u, v[:,mode])
  File =  "unit_cube_modes.vtk"
  vtkexportmesh(File, fens, fes; vectors=[("mode$mode", u.values)])
  # @async run(`"paraview.exe" $File`)
  try rm(File); catch end


  AE = AbaqusExporter("unit_cube_modes_h20");
  # AE.ios = STDOUT;
  HEADING(AE, "Vibration modes of unit cube  of almost incompressible material.");
  COMMENT(AE, "The  first six frequencies are rigid body modes.");
  COMMENT(AE, "The  first nonzero frequency (7) should be around 0.26 Hz");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  COMMENT(AE, "The hybrid form of the serendipity hexahedron is chosen because");
  COMMENT(AE, "the material is  nearly incompressible.");
  ELEMENT(AE, "c3d20rh", "AllElements", 1, connasarray(fes))
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements");
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  DENSITY(AE, rho)
  STEP_FREQUENCY(AE, neigvs)
  END_STEP(AE)
  close(AE)
  nlines = 0
  open("unit_cube_modes_h20.inp") do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 1409
  rm("unit_cube_modes_h20.inp")
end
end
using .muunit_cube_modes_exportmmm
muunit_cube_modes_exportmmm.test()

module mmpipemmPSmorthom
using FinEtools
using FinEtools.AlgoBaseModule: solve!
using FinEtoolsDeforLinear
using Test
using LinearAlgebra: norm, cholesky, cross, dot
using Statistics: mean

mutable struct MyIData
    c::FInt
    r::FFltVec
    s::FFltVec
end

function test()

## Thick pipe with internal pressure: plane strain
#

##
# Link to the  <matlab:edit('pub_thick_pipe_ps') m-file>.

## Description
##
# This is a simple modification of the full three-dimensional simulation of
# the tutorial pub_thick_pipe takes advantage of the plane-strain model
# reduction procedure.
##
# An infinitely long thick walled cylindrical pipe
# with inner boundary radius of 3 mm and outer boundary radius of 9 mm is
# subjected to an internal pressure of 1.0 MPa. A wedge   with thickness of
# 2 mm and a 90-degree angle sector is considered for the finite element
# analysis. The material properties are taken as  isotropic linear elastic
# with $E=1000$ MPa and $\nu=0.4999$ to represent nearly incompressible
# behavior. This problem has been proposed to by MacNeal and Harder as a
# test of an element's ability to represent the  response of a nearly
# incompressible material. The plane-strain condition is assumed in the
# axial direction of the pipe which together with the radial symmetry
# confines the material in all but the radial direction and therefore
# amplifies the numerical difficulties associated with the confinement of
# the nearly incompressible material.
##
# There is an analytical solution to this problem. Timoshenko and Goodier
# presented the original solution of Lame in their textbook. We are going
# to compare with  both the stress distribution (radial and hoop stresses)
# and the displacement of the inner  cylindrical surface.

##
#
# <html>
# <table border=0><tr><td>
# <img src="../docs/pub_thick_pipe_ps.png" width = "30#">
# </td></tr>
# <tr><td>Figure 1. Definition of the geometry of the internally pressurized thick pipe</td></tr>
# </table>
# </html>

##
# References:
#
# # Macneal RH, Harder RL (1985) A proposed standard set of problems to test
# finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
#
# # Timoshenko S. and Goodier J. N., Theory of Elasticity, McGraw-Hill, 2nd ed., 1951.

## Solution
#

##
# Internal radius of the pipe.
a = 3*phun("MM");
##
# External radius of the pipe.
b = 9*phun("MM");
##
# Thickness of the slice.
t = 2*phun("MM");

##
# Geometrical tolerance.
tolerance  =a/10000.;
##
# Young's modulus and Poisson's ratio.
E = 1000*phun("MEGA*PA");
nu = 0.499;
##
# Applied pressure on the internal surface.
press = 1.0*phun("MEGA*PA");

##
# Analytical solutions.   Radial stress:
radial_stress(r) =press*a.^2/(b^2-a^2).*(1-(b^2)./r.^2);
##
# Circumferential (hoop) stress:
hoop_stress(r)=press*a.^2/(b^2-a^2).*(1+(b^2)./r.^2);

##
# Radial displacement:
radial_displacement(r)=press*a^2*(1+nu)*(b^2+r.^2*(1-2*nu))/(E*(b^2-a^2).*r);;

##
# Therefore the radial displacement of the loaded surface will be:
urex = radial_displacement(a);


##
# The mesh parameters: The numbers of element edges axially,
# and through the thickness of the pipe wall (radially).

nc=3; nt=3;

##
# Note that the material object needs to be created with the proper
# model-dimension reduction in mind.  In this case that is the axial symmetry
# assumption.
MR = DeforModelRed2DStrain

# Create the mesh and initialize the geometry.  First we are going
# to construct the block of elements with the first coordinate
# corresponding to the angle, and the second
# coordinate is the thickness in the radial direction.
anglrrange = 90.0/180*pi;
fens,fes =  Q8block(anglrrange, b-a, nc, nt);

# Extract the boundary  and mark the finite elements on the
# interior surface.
bdryfes = meshboundary(fes);
bcl = selectelem(fens, bdryfes, box=[-Inf,Inf,0.,0.], inflate=tolerance);
internal_fenids = connectednodes(subset(bdryfes,bcl));
# Now  shape the block  into  the actual wedge piece of the pipe.
ayr=fens.xyz;
for i=1:count(fens)
    angl=ayr[i,1]; r=a+ayr[i,2];
    fens.xyz[i,:] = [r*sin(angl),(r*cos(angl))];
end

# now we create the geometry and displacement fields
geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),2)) # displacement field

# The symmetry boundary condition  is specified by selecting nodes
# on the plane x=0.
l1 = selectnode(fens; box=[0.0 0.0 -Inf Inf], inflate = tolerance)
setebc!(u, l1, true, 1, 0.0)
# The second symmetry boundary condition is specified by selecting
# nodes on the plane y=0.
l1 = selectnode(fens; box=[-Inf Inf 0.0 0.0], inflate = tolerance)
setebc!(u, l1, true, 2, 0.0)

applyebc!(u)
numberdofs!(u)

# The traction boundary condition is applied in the radial
# direction.

el1femm =  FEMMBase(IntegDomain(subset(bdryfes,bcl), GaussRule(1, 3)))
function pressureloading!(forceout::FFltVec, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt)
  copyto!(forceout, XYZ/norm(XYZ)*press)
  return forceout
end
fi = ForceIntensity(FFlt, 2, pressureloading!); # pressure normal to the internal cylindrical surface
F2 = distribloads(el1femm, geom, u, fi, 2);

# Property and material
material = MatDeforElastOrtho(MR, E, nu)

femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2)), material)

K = stiffness(femm, geom, u)
u = solve!(u, K, F2)

# Transfer the solution of the displacement to the nodes on the
# internal cylindrical surface and convert to
# cylindrical-coordinate displacements there.
uv = u.values[internal_fenids,:];
ur = zeros(FFlt,length(internal_fenids));

for j = 1:length(internal_fenids)
    n = fens.xyz[internal_fenids[j],:];
    n = n'/norm(n);# normal to the cylindrical internal surface
    ur[j] = dot(vec(uv[j,:]),vec(n));
end

# Report the  relative displacement on the internal surface:
# println("(Approximate/true displacement) at the internal surface: $( mean(ur)/urex*100  ) %")
@test abs(mean(ur)/urex*100 - 100) < 0.1

# Produce a plot of the radial stress component in the cylindrical
# coordinate system. Note that this is the usual representation of
# stress using nodal stress field.

fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)


File =  "thick_pipe_sigmax.vtk"
vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)])
try rm(File); catch end

# Produce a plot of the solution components in the cylindrical
# coordinate system.
# Plot the analytical solution.

function inspector(idat::MyIData, elnum, conn, xe,  out,  xq)
  function outputRm(c)
    theNormal=c;
    r=norm(theNormal);# distance from the axis of symmetry
    theNormal =theNormal/r;# compute the unit normal vector
    e1p=[theNormal';0.];# local cylind. coordinate  basis vectors
    e3p=[0.,0.,1.]';# this one points along the axis of the cylinder
    e2p=cross(vec(e3p),vec(e1p));# this one is along the hoop direction
    R= [vec(e1p) vec(e2p) vec(e3p)];# transformation matrix for the stress
    return R
  end
  Rm=outputRm(xq)
  tm=zeros(FFlt,3,3)
  stressvtot!(MR, tm, out);# stress in global XYZ
  tpm = Rm'*tm*Rm;#  stress matrix in cylindrical coordinates
  sp=zeros(FFlt,6)
  stressttov!(MR, sp, tpm);# stress vector in cylindr. coord.
  push!(idat.r,norm(xq))
  push!(idat.s,sp[idat.c])
  return idat
end

idat = MyIData(1, FFltVec[], FFltVec[])
idat = inspectintegpoints(femm, geom, u, collect(1:count(fes)),
 inspector, idat, :Cauchy)
# show(idat)

@test norm(idat.s - [-7.44858e5, -3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5,
-1.08612e5, -2.19961e5, -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6,
-7.44858e5, -3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5, -1.08612e5,
-2.19961e5, -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6, -7.44858e5,
-3.55143e5, -7.44858e5, -3.55143e5, -2.19961e5, -1.08612e5, -2.19961e5,
 -1.08612e5, -58910.9, -12517.6, -58910.9, -12517.6])/1.0e5 < 1.e-3
# using Plots
# plotly()
#
# # Plot the analytical solution.
# r = linearspace(a,b,100);
# plot(r, radial_stress(r))
# # Plot the computed  integration-point data
# scatter!(idat.r, idat.s, m=:circle, color=:red)
# gui()
end
end
using .mmpipemmPSmorthom
mmpipemmPSmorthom.test()
