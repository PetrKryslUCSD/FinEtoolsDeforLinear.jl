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
thickness = 2.5;
tolerance = 0.0001*thickness;
@show analyt_sol = 3*(1-nu^2)*Magnitude*R^2/(4*pi*E*thickness^3);

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

    	x0nl = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,x0nl,true,1,0.0)

    	y0nl = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,y0nl,true,2,0.0)

    	cylnl = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,cylnl,true,1,0.0)
    	setebc!(u,cylnl,true,2,0.0)
    	setebc!(u,cylnl,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	cnl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(cnl, length(cnl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; Magnitude/4/length(cnl)]), 3)

  		associategeometry!(femm, geom)

  		K = stiffness(femm, geom, u)

  		scattersysvec!(u, K\F)

  		u0z = mean(u.values[cnl, 3]);
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

    	x0nl = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,x0nl,true,1,0.0)

    	y0nl = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,y0nl,true,2,0.0)

    	cylnl = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,cylnl,true,1,0.0)
    	setebc!(u,cylnl,true,2,0.0)
    	setebc!(u,cylnl,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	cnl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(cnl, length(cnl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; Magnitude/4/length(cnl)]), 3)

  		associategeometry!(femm, geom)

  		K = stiffness(femm, geom, u)

  		scattersysvec!(u, K\F)

  		u0z = mean(u.values[cnl, 3]);
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

    	x0nl = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,x0nl,true,1,0.0)

    	y0nl = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,y0nl,true,2,0.0)

    	cylnl = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,cylnl,true,1,0.0)
    	setebc!(u,cylnl,true,2,0.0)
    	setebc!(u,cylnl,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	cnl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(cnl, length(cnl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; Magnitude/4/length(cnl)]), 3)

  		associategeometry!(femm, geom)
  		@show minimum(femm.phis), maximum(femm.phis)

  		K = stiffness(femm, geom, u)

  		scattersysvec!(u, K\F)

  		u0z = mean(u.values[cnl, 3]);
  		println("Deflection under the load: $(round((u0z / analyt_sol)* 100000)/100000*100) %")
  		   
    	# File =  "clcircularplatecl_1_$(nperradius).vtk"
    	# vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	# @async run(`"paraview.exe" $File`)

    end

    true
end # clcircularplatecl_h8_ms

function clcircularplatecl_h8u_ms()
    nt = 1
    for nperradius in [2, 4, 8]
    	nt=nt+1;
    	fens,fes = T4quartercyln(R, thickness, nperradius, nt)
    	fens,fes = T4toH8(fens,fes)
    	
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

    	x0nl = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,x0nl,true,1,0.0)

    	y0nl = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,y0nl,true,2,0.0)

    	cylnl = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,cylnl,true,1,0.0)
    	setebc!(u,cylnl,true,2,0.0)
    	setebc!(u,cylnl,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	cnl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(cnl, length(cnl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; Magnitude/4/length(cnl)]), 3)

  		associategeometry!(femm, geom)
  		@show minimum(femm.phis), maximum(femm.phis)
  		
  		K = stiffness(femm, geom, u)

  		scattersysvec!(u, K\F)

  		u0z = mean(u.values[cnl, 3]);
  		println("Deflection under the load: $(round((u0z / analyt_sol)* 100000)/100000*100) %")
  		   
    	File =  "clcircularplatecl_h8u_ms_$(nperradius).vtk"
    	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	@async run(`"paraview.exe" $File`)

    end

    true
end # clcircularplatecl_h8u_ms

function clcircularplatecl_h8_export()
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

    	x0nl = connectednodes(subset(boundaryfes, x0l))
    	setebc!(u,x0nl,true,1,0.0)

    	y0nl = connectednodes(subset(boundaryfes, y0l))
    	setebc!(u,y0nl,true,2,0.0)

    	cylnl = connectednodes(subset(boundaryfes, cyll))
    	setebc!(u,cylnl,true,1,0.0)
    	setebc!(u,cylnl,true,2,0.0)
    	setebc!(u,cylnl,true,3,0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	cnl = selectnode(fens; box = [0 0 0 0 0 thickness], inflate  =  tolerance)
    	
    	AE = AbaqusExporter("clcircularplatecl_h8_export_$(nperradius)");
    	HEADING(AE, "Clamped square plate with concentrated force");
    	PART(AE, "part1");
    	END_PART(AE);
    	ASSEMBLY(AE, "ASSEM1");
    	INSTANCE(AE, "INSTNC1", "PART1");
    	NODE(AE, fens.xyz);
    	ELEMENT(AE, "c3d8rh", "AllElements", 1, connasarray(femm.integdomain.fes))
    	NSET_NSET(AE, "cnl", cnl)
    	NSET_NSET(AE, "x0nl", x0nl)
    	NSET_NSET(AE, "y0nl", y0nl)
    	NSET_NSET(AE, "cylnl", cylnl)
    	ORIENTATION(AE, "GlobalOrientation", vec([1. 0 0]), vec([0 1. 0]));
    	SOLID_SECTION(AE, "elasticity", "GlobalOrientation", "AllElements", "Hourglassctl");
    	END_INSTANCE(AE);
    	END_ASSEMBLY(AE);
    	MATERIAL(AE, "elasticity")
    	ELASTIC(AE, E, nu)
    	SECTION_CONTROLS(AE, "Hourglassctl", "HOURGLASS=ENHANCED")
    	STEP_PERTURBATION_STATIC(AE)
    	BOUNDARY(AE, "ASSEM1.INSTNC1.x0nl", 1)
    	BOUNDARY(AE, "ASSEM1.INSTNC1.y0nl", 2)
    	BOUNDARY(AE, "ASSEM1.INSTNC1.cylnl", 1)
    	BOUNDARY(AE, "ASSEM1.INSTNC1.cylnl", 2)
    	BOUNDARY(AE, "ASSEM1.INSTNC1.cylnl", 3)
    	CLOAD(AE, "ASSEM1.INSTNC1.cnl", 3, Magnitude/4/length(cnl))
    	END_STEP(AE)
    	close(AE)
  		   
    	# File =  "clcircularplatecl_h8_export_$(nperradius).vtk"
    	# vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	# @async run(`"paraview.exe" $File`)

    end

    true
end # clcircularplatecl_h8_export



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
    println("#####################################################")
    println("# clcircularplatecl_h8_export ")
    clcircularplatecl_h8_export()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing