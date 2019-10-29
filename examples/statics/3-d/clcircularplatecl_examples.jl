module clcircularplatecl_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Statistics: mean

# Clamped square plate with concentrated force
#  Data listed in the Simo 1990 paper 'A class of... '
#  Analytical solution for the vertical deflection under the load
#    Zienkiewicz, Taylor, The finite element method, fifth edition, volume 2,
#     analyt_sol = 3*(1-nu^2)*Magnitude*R^2/(4*pi*E*thickness^3);
E = 1e7;
nu = 0.3;
Magnitude = 10;
R = 100.0;
thickness = 0.5;
tolerance = 0.0001*thickness;
analyt_sol = 3*(1-nu^2)*Magnitude*R^2/(4*pi*E*thickness^3);

function clcircularplatecl_h8_full()
	nt = 1
	for nperradius in [2, 4, 8]
		nt=nt+1;
    	fens, fes  = Q4circlen(R, nperradius)
    	fens, fes = H8extrudeQ4(fens, fes, nt, (x, k) -> [x[1],x[2],k*thickness/nt])

    	MR = DeforModelRed3D
    	material = MatDeforElastIso(MR, E, nu)
    	femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    	boundaryfes = meshboundary(fes);
    	topl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf thickness thickness], inflate =   tolerance);
    	botl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf 0.0 0.0], inflate =   tolerance);
    	x0l = selectelem(fens,boundaryfes, box =  [0.0 0.0 -Inf Inf 0.0 thickness], inflate =   tolerance);
    	y0l = selectelem(fens,boundaryfes, box =  [-Inf Inf 0.0 0.0 0.0 thickness], inflate =   tolerance);
    	cyll = setdiff(1:count(boundaryfes),topl,botl,x0l,y0l);

    	geom = NodalField(fens.xyz)
    	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field

    	lx1 = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,lx1,true,1,0.0)

    	lx1 = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,lx1,true,2,0.0)

    	lx1 = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,lx1,true,1,0.0)
    	setebc!(u,lx1,true,2,0.0)
    	setebc!(u,lx1,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	enl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(enl, length(enl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; Magnitude/4/length(enl)]), 3)

  		associategeometry!(femm, geom)

  		K = stiffness(femm, geom, u)

  		scattersysvec!(u, K\F)

  		u0z = mean(u.values[enl, 3]);
  		println("Deflection under the load: $(round((u0z / analyt_sol)* 100000)/100000*100) %")
  		   
    	# File =  "clcircularplatecl_1_$(nperradius).vtk"
    	# vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	# @async run(`"paraview.exe" $File`)

    end

    true
end # clcircularplatecl_h8_full

function clcircularplatecl_h8_uri()
	nt = 1
	for nperradius in [2, 4, 8]
		nt=nt+1;
		fens, fes  = Q4circlen(R, nperradius)
		fens, fes = H8extrudeQ4(fens, fes, nt, (x, k) -> [x[1],x[2],k*thickness/nt])

    	MR = DeforModelRed3D
    	material = MatDeforElastIso(MR, E, nu)
    	femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)

    	boundaryfes = meshboundary(fes);
    	topl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf thickness thickness], inflate =   tolerance);
    	botl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf 0.0 0.0], inflate =   tolerance);
    	x0l = selectelem(fens,boundaryfes, box =  [0.0 0.0 -Inf Inf 0.0 thickness], inflate =   tolerance);
    	y0l = selectelem(fens,boundaryfes, box =  [-Inf Inf 0.0 0.0 0.0 thickness], inflate =   tolerance);
    	cyll = setdiff(1:count(boundaryfes),topl,botl,x0l,y0l);

    	geom = NodalField(fens.xyz)
    	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field

    	lx1 = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,lx1,true,1,0.0)

    	lx1 = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,lx1,true,2,0.0)

    	lx1 = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,lx1,true,1,0.0)
    	setebc!(u,lx1,true,2,0.0)
    	setebc!(u,lx1,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	enl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(enl, length(enl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; Magnitude/4/length(enl)]), 3)

  		associategeometry!(femm, geom)

  		K = stiffness(femm, geom, u)

  		scattersysvec!(u, K\F)

  		u0z = mean(u.values[enl, 3]);
  		println("Deflection under the load: $(round((u0z / analyt_sol)* 100000)/100000*100) %")
  		   
    	# File =  "clcircularplatecl_1_$(nperradius).vtk"
    	# vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	# @async run(`"paraview.exe" $File`)

    end

    true
end # clcircularplatecl_h8_uri

function clcircularplatecl_h8_ms()
    nt = 1
    for nperradius in [2, 4, 8]
    	nt=nt+1;
    	fens, fes  = Q4circlen(R, nperradius)
    	fens, fes = H8extrudeQ4(fens, fes, nt, (x, k) -> [x[1],x[2],k*thickness/nt])

    	MR = DeforModelRed3D
    	material = MatDeforElastIso(MR, E, nu)
    	femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    	boundaryfes = meshboundary(fes);
    	topl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf thickness thickness], inflate =   tolerance);
    	botl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf 0.0 0.0], inflate =   tolerance);
    	x0l = selectelem(fens,boundaryfes, box =  [0.0 0.0 -Inf Inf 0.0 thickness], inflate =   tolerance);
    	y0l = selectelem(fens,boundaryfes, box =  [-Inf Inf 0.0 0.0 0.0 thickness], inflate =   tolerance);
    	cyll = setdiff(1:count(boundaryfes),topl,botl,x0l,y0l);

    	geom = NodalField(fens.xyz)
    	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field

    	lx1 = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,lx1,true,1,0.0)

    	lx1 = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,lx1,true,2,0.0)

    	lx1 = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,lx1,true,1,0.0)
    	setebc!(u,lx1,true,2,0.0)
    	setebc!(u,lx1,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	enl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(enl, length(enl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; Magnitude/4/length(enl)]), 3)

  		associategeometry!(femm, geom)

  		K = stiffness(femm, geom, u)

  		scattersysvec!(u, K\F)

  		u0z = mean(u.values[enl, 3]);
  		println("Deflection under the load: $(round((u0z / analyt_sol)* 100000)/100000*100) %")
  		   
    	# File =  "clcircularplatecl_1_$(nperradius).vtk"
    	# vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	# @async run(`"paraview.exe" $File`)

    end

    true
end # clcircularplatecl_h8_uri

# function twisted_beam_export()
#     println("""
#     Refer to twisted_beam.jl.

#     This example EXPORTS the model to Abaqus. Import the .inp file
#     into Abaqus and run the job.
#     """)
#     E = 0.29e8;
#     nu = 0.22;
#     W = 1.1;
#     L = 12.;
#     t =  0.32;
#     nl = 2; nt = 1; nw = 1; ref = 5;
#     p =   1/W/t;
#     #  Loading in the Z direction
#     loadv = [0;0;p]; dir = 3; uex = 0.005424534868469; # Harder: 5.424e-3;
#     #   Loading in the Y direction
#     #loadv = [0;p;0]; dir = 2; uex = 0.001753248285256; # Harder: 1.754e-3;
#     tolerance  = t/1000;

#     fens,fes  = H8block(L,W,t, nl*ref,nw*ref,nt*ref)

#     # Reshape into a twisted beam shape
#     for i = 1:count(fens)
#         a = fens.xyz[i,1]/L*(pi/2); y = fens.xyz[i,2]-(W/2); z = fens.xyz[i,3]-(t/2);
#         fens.xyz[i,:] = [fens.xyz[i,1],y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
#     end

#     # Clamped end of the beam
#     l1  = selectnode(fens; box = [0 0 -100*W 100*W -100*W 100*W], inflate  =  tolerance)
#     e1 = FDataDict("node_list"=>l1, "component"=>1, "displacement"=>0.0)
#     e2 = FDataDict("node_list"=>l1, "component"=>2, "displacement"=>0.0)
#     e3 = FDataDict("node_list"=>l1, "component"=>3, "displacement"=>0.0)

#     # Traction on the opposite edge
#     boundaryfes  =   meshboundary(fes);
#     Toplist   = selectelem(fens,boundaryfes, box =  [L L -100*W 100*W -100*W 100*W], inflate =   tolerance);
#     el1femm  = FEMMBase(IntegDomain(subset(boundaryfes,Toplist), GaussRule(2, 2)))
#     flux1 = FDataDict("femm"=>el1femm, "traction_vector"=>loadv)


#     # Make the region
#     MR = DeforModelRed3D
#     material = MatDeforElastIso(MR, 00.0, E, nu, 0.0)
#     region1 = FDataDict("femm"=>FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3,2)),
#     material))

#     # Make model data
#     modeldata =  FDataDict(
#     "fens"=> fens, "regions"=>  [region1],
#     "essential_bcs"=>[e1, e2, e3], "traction_bcs"=>  [flux1])


#     AE = AbaqusExporter("twisted_beam");
#     HEADING(AE, "Twisted beam example");
#     PART(AE, "part1");
#     END_PART(AE);
#     ASSEMBLY(AE, "ASSEM1");
#     INSTANCE(AE, "INSTNC1", "PART1");
#     NODE(AE, fens.xyz);
#     ELEMENT(AE, "c3d8rh", "AllElements", 1, region1["femm"].integdata.fes.conn)
#     ELEMENT(AE, "SFM3D4", "TractionElements",
#     1+count(region1["femm"].integdata.fes), flux1["femm"].integdata.fes.conn)
#     NSET_NSET(AE, "l1", l1)
#     ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
#     SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
#     SURFACE_SECTION(AE, "TractionElements")
#     END_INSTANCE(AE);
#     END_ASSEMBLY(AE);
#     MATERIAL(AE, "elasticity")
#     ELASTIC(AE, E, nu)
#     SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
#     STEP_PERTURBATION_STATIC(AE)
#     BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 1)
#     BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 2)
#     BOUNDARY(AE, "ASSEM1.INSTNC1.l1", 3)
#     DLOAD(AE, "ASSEM1.INSTNC1.TractionElements", vec(flux1["traction_vector"]))
#     END_STEP(AE)
#     close(AE)

#     true
# end # twisted_beam_export

function allrun()
    println("#####################################################")
    println("# clcircularplatecl_h8_full ")
    clcircularplatecl_h8_full()
    println("#####################################################")
    println("# clcircularplatecl_h8_uri ")
    clcircularplatecl_h8_uri()
    println("#####################################################")
    println("# clcircularplatecl_h8_ms ")
    clcircularplatecl_h8_ms()
    return true
end # function allrun

end # module clcircularplatecl_examples
