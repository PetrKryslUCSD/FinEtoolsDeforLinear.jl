module pinchcyl_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Statistics: mean

E = 3e6;
nu = 0.3;
thickness = 3.0;
uzex = -1.82488e-5; # analytical solution for the vertical deflection under the load
R = 300;
L = 600;
ref = 32
tolerance = thickness/1000;
load = [0; 0; 1.0];

function pinchcyl_h8_full()
	let (n, nt) = (ref, 2)
    	fens, fes = H8block(90/360*2*pi, L/2, thickness, n, n, nt)
    	
    	for i in 1:count(fens)
    		a = fens.xyz[i,1]; y = fens.xyz[i,2]; z = fens.xyz[i,3];
    		fens.xyz[i,:] .= ((R-thickness/2+z)*sin(a), y, (R-thickness/2+z)*cos(a));
    	end

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
    	
    	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
    	setebc!(u, y0nl, true, [1; 3], 0.0)
    	yL2nl = selectnode(fens, box = [-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance);
    	setebc!(u, yL2nl, true, [2], 0.0)
    	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
    	setebc!(u, x0nl, true, [1], 0.0)
    	z0nl = selectnode(fens, box = [-Inf Inf -Inf Inf 0 0], inflate = tolerance);
    	setebc!(u, z0nl, true, [3], 0.0)

    	applyebc!(u)
    	numberdofs!(u)

    	loadnl = selectnode(fens; box = [0 0 L/2 L/2 -Inf Inf], inflate  =  tolerance)
    	
    	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(loadnl, length(loadnl), 1)), PointRule()))
  		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; -1.0/4/length(loadnl)]), 3)

  		associategeometry!(femm, geom)

  		K = stiffness(femm, geom, u)

  		u = solve!(u, K, F)

  		u0z = mean(u.values[loadnl, 3]);
  		println("Deflection under the load: $(round((u0z / uzex)* 100000)/100000*100) %")
  		   
    	File =  "pinchcyl_h8_full_$(n)x$(nt).vtk"
    	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	@async run(`"paraview.exe" $File`)

    end

    true
end # pinchcyl_h8_full

function pinchcyl_h8_uri()
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

  		u = solve!(u, K, F)

  		u0z = mean(u.values[cnl, 3]);
  		println("Deflection under the load: $(round((u0z / analyt_sol)* 100000)/100000*100) %")
  		   
    	# File =  "pinchcyl_1_$(nperradius).vtk"
    	# vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	# @async run(`"paraview.exe" $File`)

    end

    true
end # pinchcyl_h8_uri

function pinchcyl_h8_ms()
    let (n, nt) = (ref, 4)
        	fens, fes = H8block(90/360*2*pi, L/2, thickness, n, n, nt)
        	
        	for i in 1:count(fens)
        		a = fens.xyz[i,1]; y = fens.xyz[i,2]; z = fens.xyz[i,3];
        		fens.xyz[i,:] .= ((R-thickness/2+z)*sin(a), y, (R-thickness/2+z)*cos(a));
        	end

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
        	
        	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
        	setebc!(u, y0nl, true, [1; 3], 0.0)
        	yL2nl = selectnode(fens, box = [-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance);
        	setebc!(u, yL2nl, true, [2], 0.0)
        	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
        	setebc!(u, x0nl, true, [1], 0.0)
        	z0nl = selectnode(fens, box = [-Inf Inf -Inf Inf 0 0], inflate = tolerance);
        	setebc!(u, z0nl, true, [3], 0.0)

        	applyebc!(u)
        	numberdofs!(u)

        	loadnl = selectnode(fens; box = [0 0 L/2 L/2 -Inf Inf], inflate  =  tolerance)
        	
        	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(loadnl, length(loadnl), 1)), PointRule()))
      		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; -1.0/4/length(loadnl)]), 3)

      		associategeometry!(femm, geom)

      		K = stiffness(femm, geom, u)

      		u = solve!(u, K, F)

      		u0z = mean(u.values[loadnl, 3]);
      		println("Deflection under the load: $(round((u0z / uzex)* 100000)/100000*100) %")
      		   
        	File =  "pinchcyl_h8_ms_$(n)x$(nt).vtk"
        	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
        	@async run(`"paraview.exe" $File`)

        end

    true
end # pinchcyl_h8_ms

function pinchcyl_h8_export()
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
    	
    	AE = AbaqusExporter("pinchcyl_h8_export_$(nperradius)");
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
  		   
    	# File =  "pinchcyl_h8_export_$(nperradius).vtk"
    	# vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    	# @async run(`"paraview.exe" $File`)

    end

    true
end # pinchcyl_h8_export

function pinchcyl_h20r()
    let (n, nt) = (ref, 4)
        	fens, fes = H20block(90/360*2*pi, L/2, thickness, n, n, nt)
        	
        	for i in 1:count(fens)
        		a = fens.xyz[i,1]; y = fens.xyz[i,2]; z = fens.xyz[i,3];
        		fens.xyz[i,:] .= ((R-thickness/2+z)*sin(a), y, (R-thickness/2+z)*cos(a));
        	end

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
        	
        	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
        	setebc!(u, y0nl, true, [1; 3], 0.0)
        	yL2nl = selectnode(fens, box = [-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance);
        	setebc!(u, yL2nl, true, [2], 0.0)
        	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
        	setebc!(u, x0nl, true, [1], 0.0)
        	z0nl = selectnode(fens, box = [-Inf Inf -Inf Inf 0 0], inflate = tolerance);
        	setebc!(u, z0nl, true, [3], 0.0)

        	applyebc!(u)
        	numberdofs!(u)

        	loadnl = selectnode(fens; box = [0 0 L/2 L/2 -Inf Inf], inflate  =  tolerance)
        	
        	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(loadnl, length(loadnl), 1)), PointRule()))
      		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; -1.0/4/length(loadnl)]), 3)

      		associategeometry!(femm, geom)

      		K = stiffness(femm, geom, u)

      		u = solve!(u, K, F)

      		u0z = mean(u.values[loadnl, 3]);
      		println("Deflection under the load: $(round((u0z / uzex)* 100000)/100000*100) %")
      		   
        	File =  "pinchcyl_h20r_$(n)x$(nt).vtk"
        	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
        	@async run(`"paraview.exe" $File`)

        end

    true
end # pinchcyl_h20r

function pinchcyl_h20()
    let (n, nt) = (ref, 4)
        	fens, fes = H20block(90/360*2*pi, L/2, thickness, n, n, nt)
        	
        	for i in 1:count(fens)
        		a = fens.xyz[i,1]; y = fens.xyz[i,2]; z = fens.xyz[i,3];
        		fens.xyz[i,:] .= ((R-thickness/2+z)*sin(a), y, (R-thickness/2+z)*cos(a));
        	end

        	MR = DeforModelRed3D
        	material = MatDeforElastIso(MR, E, nu)
        	femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 3)), material)

        	boundaryfes = meshboundary(fes);
        	topl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf thickness thickness], inflate =   tolerance);
        	botl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf 0.0 0.0], inflate =   tolerance);
        	x0l = selectelem(fens,boundaryfes, box =  [0.0 0.0 -Inf Inf 0.0 thickness], inflate =   tolerance);
        	y0l = selectelem(fens,boundaryfes, box =  [-Inf Inf 0.0 0.0 0.0 thickness], inflate =   tolerance);
        	cyll = setdiff(1:count(boundaryfes),topl,botl,x0l,y0l);

        	geom = NodalField(fens.xyz)
        	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
        	
        	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
        	setebc!(u, y0nl, true, [1; 3], 0.0)
        	yL2nl = selectnode(fens, box = [-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance);
        	setebc!(u, yL2nl, true, [2], 0.0)
        	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
        	setebc!(u, x0nl, true, [1], 0.0)
        	z0nl = selectnode(fens, box = [-Inf Inf -Inf Inf 0 0], inflate = tolerance);
        	setebc!(u, z0nl, true, [3], 0.0)

        	applyebc!(u)
        	numberdofs!(u)

        	loadnl = selectnode(fens; box = [0 0 L/2 L/2 -Inf Inf], inflate  =  tolerance)
        	
        	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(loadnl, length(loadnl), 1)), PointRule()))
      		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; -1.0/4/length(loadnl)]), 3)

      		associategeometry!(femm, geom)

      		K = stiffness(femm, geom, u)

      		u = solve!(u, K, F)

      		u0z = mean(u.values[loadnl, 3]);
      		println("Deflection under the load: $(round((u0z / uzex)* 100000)/100000*100) %")
      		   
        	File =  "pinchcyl_h20_$(n)x$(nt).vtk"
        	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
        	@async run(`"paraview.exe" $File`)

        end

    true
end # pinchcyl_h20

function pinchcyl_t10_ms()
    let (n, nt) = (ref, 4)
        	fens, fes = T10block(90/360*2*pi, L/2, thickness, n, n, nt)
        	
        	for i in 1:count(fens)
        		a = fens.xyz[i,1]; y = fens.xyz[i,2]; z = fens.xyz[i,3];
        		fens.xyz[i,:] .= ((R-thickness/2+z)*sin(a), y, (R-thickness/2+z)*cos(a));
        	end

        	MR = DeforModelRed3D
        	material = MatDeforElastIso(MR, E, nu)
        	femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)

        	boundaryfes = meshboundary(fes);
        	topl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf thickness thickness], inflate =   tolerance);
        	botl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf 0.0 0.0], inflate =   tolerance);
        	x0l = selectelem(fens,boundaryfes, box =  [0.0 0.0 -Inf Inf 0.0 thickness], inflate =   tolerance);
        	y0l = selectelem(fens,boundaryfes, box =  [-Inf Inf 0.0 0.0 0.0 thickness], inflate =   tolerance);
        	cyll = setdiff(1:count(boundaryfes),topl,botl,x0l,y0l);

        	geom = NodalField(fens.xyz)
        	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
        	
        	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
        	setebc!(u, y0nl, true, [1; 3], 0.0)
        	yL2nl = selectnode(fens, box = [-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance);
        	setebc!(u, yL2nl, true, [2], 0.0)
        	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
        	setebc!(u, x0nl, true, [1], 0.0)
        	z0nl = selectnode(fens, box = [-Inf Inf -Inf Inf 0 0], inflate = tolerance);
        	setebc!(u, z0nl, true, [3], 0.0)

        	applyebc!(u)
        	numberdofs!(u)

        	loadnl = selectnode(fens; box = [0 0 L/2 L/2 -Inf Inf], inflate  =  tolerance)
        	
        	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(loadnl, length(loadnl), 1)), PointRule()))
      		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; -1.0/4/length(loadnl)]), 3)

      		associategeometry!(femm, geom)

      		K = stiffness(femm, geom, u)

      		u = solve!(u, K, F)

      		u0z = mean(u.values[loadnl, 3]);
      		println("Deflection under the load: $(round((u0z / uzex)* 100000)/100000*100) %")
      		   
        	File =  "pinchcyl_t10_ms_$(n)x$(nt).vtk"
        	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
        	@async run(`"paraview.exe" $File`)

        end

    true
end # pinchcyl_t10_ms

function pinchcyl_t10()
    let (n, nt) = (ref, 4)
        	fens, fes = T10block(90/360*2*pi, L/2, thickness, n, n, nt)
        	
        	for i in 1:count(fens)
        		a = fens.xyz[i,1]; y = fens.xyz[i,2]; z = fens.xyz[i,3];
        		fens.xyz[i,:] .= ((R-thickness/2+z)*sin(a), y, (R-thickness/2+z)*cos(a));
        	end

        	MR = DeforModelRed3D
        	material = MatDeforElastIso(MR, E, nu)
        	femm = FEMMDeforLinear(MR, IntegDomain(fes, TetRule(4)), material)

        	boundaryfes = meshboundary(fes);
        	topl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf thickness thickness], inflate =   tolerance);
        	botl = selectelem(fens,boundaryfes, box =  [-Inf Inf -Inf Inf 0.0 0.0], inflate =   tolerance);
        	x0l = selectelem(fens,boundaryfes, box =  [0.0 0.0 -Inf Inf 0.0 thickness], inflate =   tolerance);
        	y0l = selectelem(fens,boundaryfes, box =  [-Inf Inf 0.0 0.0 0.0 thickness], inflate =   tolerance);
        	cyll = setdiff(1:count(boundaryfes),topl,botl,x0l,y0l);

        	geom = NodalField(fens.xyz)
        	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
        	
        	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
        	setebc!(u, y0nl, true, [1; 3], 0.0)
        	yL2nl = selectnode(fens, box = [-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance);
        	setebc!(u, yL2nl, true, [2], 0.0)
        	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
        	setebc!(u, x0nl, true, [1], 0.0)
        	z0nl = selectnode(fens, box = [-Inf Inf -Inf Inf 0 0], inflate = tolerance);
        	setebc!(u, z0nl, true, [3], 0.0)

        	applyebc!(u)
        	numberdofs!(u)

        	loadnl = selectnode(fens; box = [0 0 L/2 L/2 -Inf Inf], inflate  =  tolerance)
        	
        	nfemm = FEMMBase(IntegDomain(FESetP1(reshape(loadnl, length(loadnl), 1)), PointRule()))
      		F = distribloads(nfemm, geom, u, ForceIntensity([0; 0; -1.0/4/length(loadnl)]), 3)

      		associategeometry!(femm, geom)

      		K = stiffness(femm, geom, u)

      		u = solve!(u, K, F)

      		u0z = mean(u.values[loadnl, 3]);
      		println("Deflection under the load: $(round((u0z / uzex)* 100000)/100000*100) %")
      		   
        	File =  "pinchcyl_t10_$(n)x$(nt).vtk"
        	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
        	@async run(`"paraview.exe" $File`)

        end

    true
end # pinchcyl_t10



#     true
# end # twisted_beam_export

function allrun()
    println("#####################################################")
    println("# pinchcyl_h8_full ")
    pinchcyl_h8_full()
    println("#####################################################")
    println("# pinchcyl_h8_ms ")
    pinchcyl_h8_ms()
    println("#####################################################")
    println("# pinchcyl_t10_ms ")
    pinchcyl_t10_ms()
    println("#####################################################")
    println("# pinchcyl_h20r ")
    pinchcyl_h20r()

    return true
end # function allrun


@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing
