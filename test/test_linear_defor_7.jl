
module m111ocylpull14nnn # From miscellaneous
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using LinearAlgebra
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
    # println("Number of degrees of freedom = $(nfreedofs(u))")
    @test nfreedofs(u) == 240

    # Property and material
    material = MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    # display(material)
    # println("$(material.D)")
    # @show MR

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

    K = stiffness(femm, geom, u)

    K_ff, K_fd = matrix_blocked(K, nfreedofs(u), nfreedofs(u))[(:ff, :fd)]

    U_d = gathersysvec(u, :d)

    factor = cholesky(Symmetric(K_ff))
    U_f = factor\(-K_fd * U_d)
    scattersysvec!(u, U_f)


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
using .m111ocylpull14nnn
m111ocylpull14nnn.test()

module mophysun13 # From miscellaneous
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    E1=1.0;
    nu23=0.19;
    rin=1.;
    rex =1.2;
    Length = 1*rex
    tolerance=rin/1000.

    MR = DeforModelRed2DAxisymm

    fens,fes = Q4block(rex-rin,Length,5,20);
    fens.xyz[:, 1] = fens.xyz[:, 1] .+ rin
    bdryfes = meshboundary(fes);

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),2))
    l1 =selectnode(fens; box=[0 rex 0 0], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    l1 =selectnode(fens; box=[0 rex Length Length], inflate = tolerance)
    setebc!(u,l1,true, 2, 0.0)
    applyebc!(u)
    numberdofs!(u)
    @test nfreedofs(u) == 240

    material=MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)
    # println("success? ")
    # @code_llvm FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material, true)
    # femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material, true)
    # println("failure? ")
    # @code_llvm FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)

    true
end
end
using .mophysun13
mophysun13.test()

module m111ocylpull14n1 # From miscellaneous
using FinEtools
using FinEtoolsDeforLinear
using Test

function test()
    E1=1.0;
    nu23=0.19;
    rin=1.;
    rex =1.2;
    Length = 1*rex

    MR = DeforModelRed2DAxisymm
    fens,fes = Q4block(rex-rin,Length,5,20);
    material = MatDeforElastIso(MR, 00.0, E1, nu23, 0.0)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(2, 2), true), material)
    femm.mr == MR

    true
end
end
using .m111ocylpull14n1
m111ocylpull14n1.test()

module mocylpull14a
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    MR = DeforModelRed2DAxisymm
    material = MatDeforElastIso(MR, 0.0, 1.0, 0.0, 0.0)
    @test MR == material.mr
    femm = FEMMDeforLinear(MR, IntegDomain(FESetP1(reshape([1],1,1)), GaussRule(2, 2), true), material)

end
end
using .mocylpull14a
mocylpull14a.test()


module mstressconversionm
using FinEtools
using FinEtoolsDeforLinear
using Test
import LinearAlgebra: norm
function test()
  symmtens(N) = begin t=rand(N, N); t = (t+t')/2.0; end
  t = symmtens(2)
  v = zeros(3)
  strainttov!(DeforModelRed2DStrain, v, t)
  to = zeros(2, 2)
  strainvtot!(DeforModelRed2DStrain, to, v)
  @test norm(t-to) < eps(1.0)

  t = symmtens(3)
  v = zeros(6)
  strainttov!(DeforModelRed3D, v, t)
  to = zeros(3, 3)
  strainvtot!(DeforModelRed3D, to, v)
  @test norm(t-to) < eps(1.0)

  t = symmtens(2)
  v = zeros(3)
  stressttov!(DeforModelRed2DStress, v, t)
  to = zeros(2, 2)
  stressvtot!(DeforModelRed2DStress, to, v)
  @test norm(t-to) < eps(1.0)

  v = vec([1. 2. 3.])
  t = zeros(3, 3)
  stressvtot!(DeforModelRed2DStrain, t, v)
  to = [1. 3. 0; 3. 2. 0; 0 0 0]
  @test norm(t-to) < eps(1.0)

  v = vec([1. 2 3 4])
  t = zeros(3, 3)
  stressvtot!(DeforModelRed2DStrain, t, v)
  to = [1. 3 0; 3 2 0; 0 0 4]
  @test norm(t-to) < eps(1.0)

  v = rand(6)
  t = zeros(3, 3)
  stressvtot!(DeforModelRed3D, t, v)
  vo = zeros(6)
  stressttov!(DeforModelRed3D, vo, t)
  @test norm(v-vo) < eps(1.0)

  # v = rand(9)
  # t = zeros(3, 3)
  # strain9vto3x3t!(t, v)
  # t = (t + t')/2.0 # symmetrize
  # strain3x3tto9v!(v, t)
  # v6 = zeros(6)
  # strain9vto6v!(v6, v)
  # v9 = zeros(9)
  # strain6vto9v!(v9, v6)
  # @test norm(v-v9) < eps(1.0)

  # v = vec([1. 2 3 4 4 5 5 6 6])
  # v6 = zeros(6)
  # stress9vto6v!(v6, v)
  # v9 = zeros(9)
  # stress6vto9v!(v9, v6)
  # @test norm(v-v9) < eps(1.0)

end
end
using .mstressconversionm
mstressconversionm.test()


module mmmtdt1a
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
function test()
	MR = DeforModelRed3D
	symmet(a) = a + transpose(a)
	
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(3, 3)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6

	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6

	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
	a = symmet(rand(3, 3))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
end
end
using .mmmtdt1a
mmmtdt1a.test()

module mmmtdt2
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
function test()
	MR = DeforModelRed2DStrain
	symmet(a) = a + transpose(a)
	
	a = rand(2, 2)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(2, 2)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(2, 2)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(2, 2)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(2, 2)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = rand(2, 2)
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	@test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6

	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 6)
	strainttov!(MR, av, a)
	@test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6

	MR = DeforModelRed2DStrain
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 3)
	strainttov!(MR, av, a)
	@test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 3)
	strainttov!(MR, av, a)
	@test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
	a = symmet(rand(2, 2))
	av = fill(zero(eltype(a)), 3)
	strainttov!(MR, av, a)
	@test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
end
end
using .mmmtdt2
mmmtdt2.test()

module mmtwistedeexportmm
using FinEtools
using FinEtoolsDeforLinear
using Test
using FinEtools.MeshExportModule
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 3;
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


  AE = AbaqusExporter("twisted_beam");
  HEADING(AE, "Twisted beam example");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(region1["femm"].integdomain.fes))
  ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(region1["femm"].integdomain.fes), connasarray(flux1["femm"].integdomain.fes))
  NSET_NSET(AE, "l1", l1)
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglass");
  SURFACE_SECTION(AE, "TractionElements")
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  SECTION_CONTROLS(AE, "section1", "HOURGLASS=ENHANCED")
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3)
  DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
  END_STEP(AE)
  close(AE)
  nlines = 0
  open("twisted_beam.inp") do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 223
  rm("twisted_beam.inp")

  true
end
end
using .mmtwistedeexportmm
mmtwistedeexportmm.test()


module mmtwistedeexport2mm
using FinEtools
using FinEtoolsDeforLinear
using Test
using FinEtools.MeshExportModule
function test()
  E = 0.29e8;
  nu = 0.22;
  W = 1.1;
  L = 12.;
  t =  0.32;
  nl = 2; nt = 1; nw = 1; ref = 3;
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


  AE = AbaqusExporter("twisted_beam");
  HEADING(AE, "Twisted beam example");
  PART(AE, "part1");
  END_PART(AE);
  ASSEMBLY(AE, "ASSEM1");
  INSTANCE(AE, "INSTNC1", "PART1");
  NODE(AE, fens.xyz);
  ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(region1["femm"].integdomain.fes))
  ELEMENT(AE, "SFM3D4", "TractionElements",
    1+count(region1["femm"].integdomain.fes), connasarray(flux1["femm"].integdomain.fes))
  NSET_NSET(AE, "l1", l1)
  ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
  SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglass");
  SURFACE_SECTION(AE, "TractionElements")
  END_INSTANCE(AE);
  END_ASSEMBLY(AE);
  MATERIAL(AE, "elasticity")
  ELASTIC(AE, E, nu)
  SECTION_CONTROLS(AE, "section1", "HOURGLASS=ENHANCED")
  STEP_PERTURBATION_STATIC(AE)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1, 0.0)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2, 0.0)
  BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3, 0.0)
  DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
  END_STEP(AE)
  close(AE)
  nlines = 0
  open("twisted_beam.inp") do f
    s = readlines(f)
    nlines = length(s)
  end
  @test nlines == 223
  rm("twisted_beam.inp")

  true
end
end
using .mmtwistedeexport2mm
mmtwistedeexport2mm.test()


module mstresscomponentmap
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    MR = DeforModelRed1D
    @test stresscomponentmap(MR)[:x] == 1
    MR = DeforModelRed2DAxisymm
    @test stresscomponentmap(MR)[:x] == 1
    @test stresscomponentmap(MR)[:zz] == 3
end
end
using .mstresscomponentmap
mstresscomponentmap.test()


module mmMeasurement_3a
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    W = 1.1;
    L = 12.;
    t =  3.32;
    nl, nt, nw = 5, 3, 4;
    Ea =  210000*phun("MEGA*Pa")
    nua =  0.3;

    # println("New segmentation fault?")
    for orientation in [:a :b :ca :cb]
        fens,fes  = T4block(L,W,t, nl,nw,nt, orientation)
        geom  =  NodalField(fens.xyz)
        MR  =  DeforModelRed3D
        material = MatDeforElastIso(MR, 0.0, Ea, nua, 0.0)

        femm  =  FEMMDeforLinearNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
        V = integratefunction(femm, geom, (x) ->  1.0)
        @test abs(V - W*L*t)/V < 1.0e-5
    end

end
end
using .mmMeasurement_3a
mmMeasurement_3a.test()

module mmvgradmat1
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: infsup_sh
using Test
import LinearAlgebra: norm, cholesky
function test()
	Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = ( 1.0, 1.0, 1.0, 1, 1, 1, :a)
	Ea, nua, alphaa = ( 1.0, 0.3, 0.0)
	fens = FENodeSet(Float64[0     0     0
     0     3     3
     0     0     3
     3     0     3]);
	fes = FESetT4(reshape([1 2 3 4], 1, 4));


	MR  =  DeforModelRed3D

	# Property and material
	material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

	femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
	numberdofs!(u)

	Sh = infsup_sh(femm, geom, u);

	Sh1 = [0.5   0   0   0   0   0   -0.5   0   0   -0.0000   0   0
		   0    0.5   0   0   0   0   0   -0.5   0   0   -0.0000   0
		   0   0    0.5   0   0   0   0   0   -0.5   0   0   -0.0000
		   0   0   0    0.5   0   0   -0.5   0   0   0   0   0
		   0   0   0   0    0.5   0   0   -0.5   0   0   0   0
		   0   0   0   0   0    0.5   0   0   -0.5   0   0   0
		   -0.5   0   0   -0.5   0   0    1.5000   0   0   -0.5   0   0
		   0   -0.5   0   0   -0.5   0   0    1.5000   0   0   -0.5   0
		   0   0   -0.5   0   0   -0.5   0   0    1.5000   0   0   -0.5
		   -0.0000   0   0   0   0   0   -0.5   0   0    0.5   0   0
		   0   -0.0000   0   0   0   0   0   -0.5   0   0    0.5   0
		   0   0   -0.0000   0   0   0   0   0   -0.5   0   0    0.5]
	@test norm(Sh - Sh1)  <= 1.0e-9
end
end
using .mmvgradmat1
mmvgradmat1.test()

module mmdivmat2
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: infsup_gh
using Test
import LinearAlgebra: norm, cholesky
function test()
	Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = ( 1.0, 1.0, 1.0, 1, 1, 1, :a)
	Ea, nua, alphaa = ( 1.0, 0.3, 0.0)
	fens = FENodeSet(Float64[           0            0            0
							  -1.6200e-01   3.0000e+00   2.9791e+00
							            0            0   3.0000e+00
							   3.0000e+00   1.6200e-01   3.0000e+00]);
	fes = FESetT4(reshape([1 2 3 4], 1, 4));


	MR  =  DeforModelRed3D

	# Property and material
	material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

	femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(4)), material)

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
	numberdofs!(u)

	Gh = infsup_gh(femm, geom, u);

	Gh1 = [  7.0319e-08  -1.3022e-06  -1.8778e-04  -1.0111e-05   1.8724e-04            0  -1.7720e-04  -1.9604e-04   1.8778e-04   1.8724e-04   1.0111e-05            0
  -1.3022e-06   2.4115e-05   3.4774e-03   1.8724e-04  -3.4673e-03            0   3.2814e-03   3.6305e-03  -3.4774e-03  -3.4673e-03  -1.8724e-04            0
  -1.8778e-04   3.4774e-03   5.0146e-01   2.7000e-02  -5.0000e-01            0   4.7319e-01   5.2352e-01  -5.0146e-01  -5.0000e-01  -2.7000e-02            0
  -1.0111e-05   1.8724e-04   2.7000e-02   1.4538e-03  -2.6921e-02            0   2.5478e-02   2.8188e-02  -2.7000e-02  -2.6921e-02  -1.4538e-03            0
   1.8724e-04  -3.4673e-03  -5.0000e-01  -2.6921e-02   4.9855e-01            0  -4.7181e-01  -5.2200e-01   5.0000e-01   4.9855e-01   2.6921e-02            0
            0            0            0            0            0            0            0            0            0            0            0            0
  -1.7720e-04   3.2814e-03   4.7319e-01   2.5478e-02  -4.7181e-01            0   4.4651e-01   4.9401e-01  -4.7319e-01  -4.7181e-01  -2.5478e-02            0
  -1.9604e-04   3.6305e-03   5.2352e-01   2.8188e-02  -5.2200e-01            0   4.9401e-01   5.4656e-01  -5.2352e-01  -5.2200e-01  -2.8188e-02            0
   1.8778e-04  -3.4774e-03  -5.0146e-01  -2.7000e-02   5.0000e-01            0  -4.7319e-01  -5.2352e-01   5.0146e-01   5.0000e-01   2.7000e-02            0
   1.8724e-04  -3.4673e-03  -5.0000e-01  -2.6921e-02   4.9855e-01            0  -4.7181e-01  -5.2200e-01   5.0000e-01   4.9855e-01   2.6921e-02            0
   1.0111e-05  -1.8724e-04  -2.7000e-02  -1.4538e-03   2.6921e-02            0  -2.5478e-02  -2.8188e-02   2.7000e-02   2.6921e-02   1.4538e-03            0
            0            0            0            0            0            0            0            0            0            0            0            0]
    # @show Matrix(Gh)
	@test norm(Gh - Gh1) / norm(Gh1)  <= 2.0e-5
end
end
using .mmdivmat2
mmdivmat2.test()

module mmvgradmat2
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: infsup_sh
using Test
import LinearAlgebra: norm, cholesky
function test()
	Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = ( 1.0, 1.0, 1.0, 1, 1, 1, :a)
	Ea, nua, alphaa = ( 1.0, 0.3, 0.0)
	fens = FENodeSet(Float64[        0         0         0
								   -0.1620    3.0000    2.9791
								         0         0    3.0000
								    3.0000    0.1620    3.0000]);
	fes = FESetT4(reshape([1 2 3 4], 1, 4));


	MR  =  DeforModelRed3D

	# Property and material
	material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

	femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
	numberdofs!(u)

	Sh = infsup_sh(femm, geom, u);

	Sh1 = [5.0148e-01            0            0  -3.4774e-03            0            0  -4.9800e-01            0            0  -2.6321e-17            0            0
            0   5.0148e-01            0            0  -3.4774e-03            0            0  -4.9800e-01            0            0  -2.6321e-17            0
            0            0   5.0148e-01            0            0  -3.4774e-03            0            0  -4.9800e-01            0            0  -2.6321e-17
  -3.4774e-03            0            0   5.0000e-01            0            0  -4.9652e-01            0            0  -4.7374e-18            0            0
            0  -3.4774e-03            0            0   5.0000e-01            0            0  -4.9652e-01            0            0  -4.7374e-18            0
            0            0  -3.4774e-03            0            0   5.0000e-01            0            0  -4.9652e-01            0            0  -4.7374e-18
  -4.9800e-01            0            0  -4.9652e-01            0            0   1.4945e+00            0            0  -5.0000e-01            0            0
            0  -4.9800e-01            0            0  -4.9652e-01            0            0   1.4945e+00            0            0  -5.0000e-01            0
            0            0  -4.9800e-01            0            0  -4.9652e-01            0            0   1.4945e+00            0            0  -5.0000e-01
  -2.6321e-17            0            0  -4.7374e-18            0            0  -5.0000e-01            0            0   5.0000e-01            0            0
            0  -2.6321e-17            0            0  -4.7374e-18            0            0  -5.0000e-01            0            0   5.0000e-01            0
            0            0  -2.6321e-17            0            0  -4.7374e-18            0            0  -5.0000e-01            0            0   5.0000e-01]
	@test norm(Sh - Sh1) / norm(Sh1)  <= 1.0e-4
end
end
using .mmvgradmat2
mmvgradmat2.test()

module mmdivmat1
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: infsup_gh
using Test
import LinearAlgebra: norm, cholesky
function test()
	Length::FFlt, Width::FFlt, Height::FFlt, nL::FInt, nW::FInt, nH::FInt, orientation::Symbol = ( 1.0, 1.0, 1.0, 1, 1, 1, :a)
	Ea, nua, alphaa = ( 1.0, 0.3, 0.0)
	fens = FENodeSet(Float64[0     0     0
     0     3     3
     0     0     3
     3     0     3]);
	fes = FESetT4(reshape([1 2 3 4], 1, 4));


	MR  =  DeforModelRed3D

	# Property and material
	material = MatDeforElastIso(MR, 0.0, Ea, nua, alphaa)

	femm  =  FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
	numberdofs!(u)

	Gh = infsup_gh(femm, geom, u);

	Gh1 = [0.0  0    0.0  0   -0.0  0    0.0    0.0   -0.0   -0.0  0  0
			0  0  0  0  0  0  0  0  0  0  0  0
			0.0  0    0.5  0   -0.5  0    0.5    0.5   -0.5   -0.5  0  0
			0  0  0  0  0  0  0  0  0  0  0  0
			-0.0  0   -0.5  0    0.5  0   -0.5   -0.5    0.5    0.5  0  0
			0  0  0  0  0  0  0  0  0  0  0  0
			0.0  0    0.5  0   -0.5  0    0.5    0.5   -0.5   -0.5  0  0
			0.0  0    0.5  0   -0.5  0    0.5    0.5   -0.5   -0.5  0  0
			-0.0  0   -0.5  0    0.5  0   -0.5   -0.5    0.5    0.5  0  0
			-0.0  0   -0.5  0    0.5  0   -0.5   -0.5    0.5    0.5  0  0
			0  0  0  0  0  0  0  0  0  0  0  0
			0  0  0  0  0  0  0  0  0  0  0  0]
	@test norm(Gh - Gh1)  <= 1.0e-9
end
end
using .mmdivmat1
mmdivmat1.test()


module mmmMergem1
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;

    t0 = time()
    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, box = boundingbox([Rmed-h -Inf 0.0; Rmed+h +Inf 0.0]), inflate = tolerance)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    @test count(fens) == 8820

    # fens,fes = mergenodes(fens, fes,  tolerance);
    # println("  after merging  = $(count(fens))")
    #
    # println("Mesh generation ($(time() - t0) sec)")
    #
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
end
using .mmmMergem1
 mmmMergem1.test()


module mmmMergem2
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;

    t0 = time()
    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], inflate = tolerance)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    @test count(fens) == 8820

    # fens,fes = mergenodes(fens, fes,  tolerance);
    # println("  after merging  = $(count(fens))")
    #
    # println("Mesh generation ($(time() - t0) sec)")
    #
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
end
using .mmmMergem2
 mmmMergem2.test()


module mmmMergem3
using FinEtools
using FinEtoolsDeforLinear
using Test
function test()
    h = 0.05*phun("M");
    l = 10*h;
    Rmed = h/0.2;
    psi   = 0;    # Cylinder
    nh = 12; nl  = 40; nc = 120;
    nh = 6; nl  = 20; nc = 60;
    # nh = 3; nl  = 8; nc = 30;
    tolerance = h/nh/100;

    t0 = time()
    MR = DeforModelRed3D
    fens,fes  = H8block(h,l,2.0*pi,nh,nl,nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i,:];
        rotmat3!(R, [0, z, 0])
        Q = [cos(psi*pi/180) sin(psi*pi/180) 0;
            -sin(psi*pi/180) cos(psi*pi/180) 0;
            0 0 1]
        fens.xyz[i,:] = reshape([x+Rmed-h/2, y-l/2, 0], 1, 3)*Q*R;
    end
    # println("  before merging  = $(count(fens))")
    @test count(fens) == 8967
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h/1000)
    fens,fes = mergenodes(fens, fes,  tolerance, candidates);
    @test count(fens) == 8820

    # fens,fes = mergenodes(fens, fes,  tolerance);
    # println("  after merging  = $(count(fens))")
    #
    # println("Mesh generation ($(time() - t0) sec)")
    #
    # File =  "unit_cube_modes.vtk"
    # vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)
end
end
using .mmmMergem3
 mmmMergem3.test()


 module mtwistedalgofin # From miscellaneous
 using FinEtools
 using FinEtoolsDeforLinear
 using FinEtoolsDeforLinear.AlgoDeforLinearModule
 using Statistics
 using LinearAlgebra: norm
 using Test
 function test()
     # println("""
     # The initially twisted cantilever beam is one of the standard test    problems for verifying the finite-element accuracy [1]. The beam is        clamped at one end and loaded either with unit in-plane or        unit out-of-plane force at the other. The centroidal axis of the beam s
     # straight at the undeformed  configuration, while its cross-sections are
     # twisted about the centroidal axis from 0 at the clamped end to pi/2 at
     # the free end.
     #
     # Reference:
     # Zupan D, Saje M (2004) On "A proposed standard set of problems to test
     # finite element accuracy": the twisted beam. Finite Elements in Analysis
     # and Design 40: 1445-1451.
     # """)
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
     flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>ForceIntensity(loadv))


     # Make the region
     MR = DeforModelRed3D
     material = MatDeforElastIso(MR, 0.0, E, nu, 0.0)
     region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),    material))

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
     @test (theutip[dir] - 0.005443006890338948) / 0.005443006890338948 <= 1.0e-6

     # Write out mesh with displacements
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam")
     modeldata = AlgoDeforLinearModule.exportdeformation(modeldata)

     # Write out mesh with stresses
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
     "quantity"=> :Cauchy, "component"=> :xy)
     modeldata = AlgoDeforLinearModule.exportstress(modeldata)

     # Write out mesh with stresses
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
     "quantity"=> :Cauchy, "component"=> :xz)
     modeldata = AlgoDeforLinearModule.exportstress(modeldata)

     # Write out mesh with von Mises stresses
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam",
     "quantity"=> :vm)
     modeldata = AlgoDeforLinearModule.exportstress(modeldata)

     # Write out mesh with von Mises stresses, elementwise
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
     "quantity"=> :vm)
     modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
     vm  = modeldata["postprocessing"]["exported"][1]["field"]
     # println("extremes of von Mises: $([minimum(vm.values),   maximum(vm.values)])")
     @test norm([minimum(vm.values),   maximum(vm.values)] - [5.21266051216496, 452.7873821727328]) <= 1.0e-3

     # Write out mesh with von Mises stresses, elementwise
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-ew",
     "quantity"=> :Cauchy, "component"=> :xz)
     modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

     # Write out mesh with principal stresses, elementwise
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-1-ew",
     "quantity"=> :princCauchy, "component"=> 1)
     modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
     ps  = modeldata["postprocessing"]["exported"][1]["field"]
     # println("extremes of first principal stress: $([minimum(ps.values),   maximum(ps.values)])")
     @test norm([minimum(ps.values),   maximum(ps.values)] - [-5.9919354276389285, 461.8866079275564]) <= 1.0e-3

     # Write out mesh with principal stresses, elementwise
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-principal-3-ew",
     "quantity"=> :princCauchy, "component"=> 3)
     modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
     ps  = modeldata["postprocessing"]["exported"][1]["field"]
     # println("extremes of third principal stress: $([minimum(ps.values),   maximum(ps.values)])")
     @test norm([minimum(ps.values),   maximum(ps.values)] - [-461.88660792935354, 5.9919354278587855]) <= 1.0e-3

     # Write out mesh with principal stresses, elementwise
     modeldata["postprocessing"] = FDataDict("file"=>"twisted_beam-press-ew",
     "quantity"=> :pressure, "component"=> 1)
     modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)
     ps  = modeldata["postprocessing"]["exported"][1]["field"]
     # println("extremes of pressure: $([minimum(ps.values),   maximum(ps.values)])")
     @test norm([minimum(ps.values),   maximum(ps.values)] - [-167.8998735447784, 167.89987354544647]) <= 1.0e-3
     try
         for f    in  ["twisted_beam-Cauchyxy-region 1.vtk",  "twisted_beam-press-ew-pressure1-region 1.vtk",
             "twisted_beam-Cauchyxz-region 1.vtk",     "twisted_beam-principal-1-ew-princCauchy1-region 1.vtk",
             "twisted_beam-ew-Cauchyxz-region 1.vtk",  "twisted_beam-principal-3-ew-princCauchy3-region 1.vtk",
             "twisted_beam1.vtk",  "twisted_beam-ew-vm1-region 1.vtk",       "twisted_beam-vm1-region 1.vtk"]
             rm(f)
         end
     catch end

     # println("Done")
     true
 end
 end
 using .mtwistedalgofin
 mtwistedalgofin.test()
