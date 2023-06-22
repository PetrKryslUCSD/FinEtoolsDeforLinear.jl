module bushing_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve!, matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using Statistics: mean
 
E = 3.0;
nu = 0.4999999; # This is the correct Poisson ratio
Ri = 0.25;
Re = 1.0;
L = Ri/2;
ang = 180/180*pi;
iuz = Ri/15
p = 0.27;
tolerance = min(L, Re-Ri, ang)/1000;
nR, nc, nt = 7, 14, 1

function bushing_h8_full()

	fens, fes  = H8block(ang, L, Re-Ri, nc, nt, nR)
	internal_fenids = selectnode(fens, box = [0 ang 0 L 0 0], inflate = tolerance);

	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, E, nu)
	femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

	boundaryfes = meshboundary(fes);
	ibcl = selectelem(fens, boundaryfes, box = [0 ang 0 L 0 0], inflate = tolerance);
	ebcl = selectelem(fens, boundaryfes, box = [0 ang 0 L Re-Ri Re-Ri], inflate = tolerance);
	
	for i in 1:count(fens)
		a = fens.xyz[i,1]; y = fens.xyz[i,2]; r = fens.xyz[i,3]+Ri;
		fens.xyz[i,:] .= (r*sin(a), y, r*cos(a));
	end

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field

	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
	setebc!(u, x0nl, true, 1, 0.0)

	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
	setebc!(u, y0nl, true, 2, 0.0)
	yLnl = selectnode(fens, box = [-Inf Inf L L -Inf Inf], inflate = tolerance);
	setebc!(u, yLnl, true, 2, 0.0)

	enl = connectednodes(subset(boundaryfes, ebcl))
	setebc!(u, enl, true, [1; 2; 3], 0.0)
	inl = connectednodes(subset(boundaryfes, ibcl))
	setebc!(u, inl, true, [1; 2], 0.0)
	setebc!(u, inl, true, [3], iuz)

	applyebc!(u)
	numberdofs!(u)

	associategeometry!(femm, geom)

	K = stiffness(femm, geom, u)

    u = solve!(u, K, fill(0.0, nalldofs(u)))

	File =  "bushing_h8_full_u.vtk"
	vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
	@async run(`"paraview.exe" $File`)

    true
end # bushing_h8_full

function bushing_h8_algo_full()

	fens, fes  = H8block(ang, L, Re-Ri, nc, nt, nR)
	internal_fenids = selectnode(fens, box = [0 ang 0 L 0 0], inflate = tolerance);

	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, E, nu)
	femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
	regions = FDataDict[FDataDict("femm"=>femm)]

	boundaryfes = meshboundary(fes);
	ibcl = selectelem(fens, boundaryfes, box = [0 ang 0 L 0 0], inflate = tolerance);
	ebcl = selectelem(fens, boundaryfes, box = [0 ang 0 L Re-Ri Re-Ri], inflate = tolerance);
	
	for i in 1:count(fens)
		a = fens.xyz[i,1]; y = fens.xyz[i,2]; r = fens.xyz[i,3]+Ri;
		fens.xyz[i,:] .= (r*sin(a), y, r*cos(a));
	end

	essential_bcs = FDataDict[]
	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>x0nl))

	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>y0nl))
	yLnl = selectnode(fens, box = [-Inf Inf L L -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>yLnl))

	enl = connectednodes(subset(boundaryfes, ebcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 3, "node_list"=>enl))
	inl = connectednodes(subset(boundaryfes, ibcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> iuz, "component"=> 3, "node_list"=>inl))

	modeldata = FDataDict("fens"=>fens,
	"regions"=>regions,
	"essential_bcs"=>essential_bcs
	)
	modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

	filebase = "bushing_h8_algo_full"
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-s", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstress(modeldata)

    true
end # bushing_h8_algo_full

function bushing_h8_algo_ms()

	fens, fes  = H8block(ang, L, Re-Ri, nc, nt, nR)
	internal_fenids = selectnode(fens, box = [0 ang 0 L 0 0], inflate = tolerance);

	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, E, nu)
	femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)
	regions = FDataDict[FDataDict("femm"=>femm)]

	boundaryfes = meshboundary(fes);
	ibcl = selectelem(fens, boundaryfes, box = [0 ang 0 L 0 0], inflate = tolerance);
	ebcl = selectelem(fens, boundaryfes, box = [0 ang 0 L Re-Ri Re-Ri], inflate = tolerance);
	
	for i in 1:count(fens)
		a = fens.xyz[i,1]; y = fens.xyz[i,2]; r = fens.xyz[i,3]+Ri;
		fens.xyz[i,:] .= (r*sin(a), y, r*cos(a));
	end

	essential_bcs = FDataDict[]
	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>x0nl))

	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>y0nl))
	yLnl = selectnode(fens, box = [-Inf Inf L L -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>yLnl))

	enl = connectednodes(subset(boundaryfes, ebcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 3, "node_list"=>enl))
	inl = connectednodes(subset(boundaryfes, ibcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> iuz, "component"=> 3, "node_list"=>inl))

	modeldata = FDataDict("fens"=>fens,
	"regions"=>regions,
	"essential_bcs"=>essential_bcs
	)
	modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

	filebase = "bushing_h8_algo_ms"
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-nodal-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstress(modeldata)
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-elwise-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

    true
end # bushing_h8_algo_ms

function bushing_h8u_algo_ms()

	fens, fes  = T4block(ang, L, Re-Ri, Int.(round.((nc, nt, nR)./2).+1)...)
	fens, fes  = T4toH8(fens, fes)
	internal_fenids = selectnode(fens, box = [0 ang 0 L 0 0], inflate = tolerance);

	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, E, nu)
	femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)
	regions = FDataDict[FDataDict("femm"=>femm)]

	boundaryfes = meshboundary(fes);
	ibcl = selectelem(fens, boundaryfes, box = [0 ang 0 L 0 0], inflate = tolerance);
	ebcl = selectelem(fens, boundaryfes, box = [0 ang 0 L Re-Ri Re-Ri], inflate = tolerance);
	
	for i in 1:count(fens)
		a = fens.xyz[i,1]; y = fens.xyz[i,2]; r = fens.xyz[i,3]+Ri;
		fens.xyz[i,:] .= (r*sin(a), y, r*cos(a));
	end

	essential_bcs = FDataDict[]
	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>x0nl))

	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>y0nl))
	yLnl = selectnode(fens, box = [-Inf Inf L L -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>yLnl))

	enl = connectednodes(subset(boundaryfes, ebcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 3, "node_list"=>enl))
	inl = connectednodes(subset(boundaryfes, ibcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> iuz, "component"=> 3, "node_list"=>inl))

	modeldata = FDataDict("fens"=>fens,
	"regions"=>regions,
	"essential_bcs"=>essential_bcs
	)
	modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

	filebase = "bushing_h8u_algo_ms"
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-nodal-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstress(modeldata)
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-elwise-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

    true
end # bushing_h8u_algo_ms

function bushing_h8_export()
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
end # bushing_h8_export

function bushing_t10_algo_ms()

	fens, fes  = T10block(ang, L, Re-Ri, nc, nt, nR)
	internal_fenids = selectnode(fens, box = [0 ang 0 L 0 0], inflate = tolerance);

	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, E, nu)
	femm = FEMMDeforLinearMST10(MR, IntegDomain(fes, TetRule(4)), material)
	regions = FDataDict[FDataDict("femm"=>femm)]

	boundaryfes = meshboundary(fes);
	ibcl = selectelem(fens, boundaryfes, box = [0 ang 0 L 0 0], inflate = tolerance);
	ebcl = selectelem(fens, boundaryfes, box = [0 ang 0 L Re-Ri Re-Ri], inflate = tolerance);
	
	for i in 1:count(fens)
		a = fens.xyz[i,1]; y = fens.xyz[i,2]; r = fens.xyz[i,3]+Ri;
		fens.xyz[i,:] .= (r*sin(a), y, r*cos(a));
	end

	essential_bcs = FDataDict[]
	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>x0nl))

	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>y0nl))
	yLnl = selectnode(fens, box = [-Inf Inf L L -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>yLnl))

	enl = connectednodes(subset(boundaryfes, ebcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 3, "node_list"=>enl))
	inl = connectednodes(subset(boundaryfes, ibcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> iuz, "component"=> 3, "node_list"=>inl))

	modeldata = FDataDict("fens"=>fens,
	"regions"=>regions,
	"essential_bcs"=>essential_bcs
	)
	modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

	filebase = "bushing_T10_algo_ms"
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-nodal-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstress(modeldata)
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-elwise-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

    true
end # bushing_t10_algo_ms

function bushing_t10_algo()

	fens, fes  = T10block(ang, L, Re-Ri, nc, nt, nR)
	internal_fenids = selectnode(fens, box = [0 ang 0 L 0 0], inflate = tolerance);

	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, E, nu)
	femm = FEMMDeforLinear(MR, IntegDomain(fes, TetRule(4)), material)
	regions = FDataDict[FDataDict("femm"=>femm)]

	boundaryfes = meshboundary(fes);
	ibcl = selectelem(fens, boundaryfes, box = [0 ang 0 L 0 0], inflate = tolerance);
	ebcl = selectelem(fens, boundaryfes, box = [0 ang 0 L Re-Ri Re-Ri], inflate = tolerance);
	
	for i in 1:count(fens)
		a = fens.xyz[i,1]; y = fens.xyz[i,2]; r = fens.xyz[i,3]+Ri;
		fens.xyz[i,:] .= (r*sin(a), y, r*cos(a));
	end

	essential_bcs = FDataDict[]
	x0nl = selectnode(fens, box = [0 0 -Inf Inf -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>x0nl))

	y0nl = selectnode(fens, box = [-Inf Inf 0 0 -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>y0nl))
	yLnl = selectnode(fens, box = [-Inf Inf L L -Inf Inf], inflate = tolerance);
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>yLnl))

	enl = connectednodes(subset(boundaryfes, ebcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>enl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 3, "node_list"=>enl))
	inl = connectednodes(subset(boundaryfes, ibcl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 1, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> 0.0, "component"=> 2, "node_list"=>inl))
	push!(essential_bcs, FDataDict("displacement"=> iuz, "component"=> 3, "node_list"=>inl))

	modeldata = FDataDict("fens"=>fens,
	"regions"=>regions,
	"essential_bcs"=>essential_bcs
	)
	modeldata = AlgoDeforLinearModule.linearstatics(modeldata)

	filebase = "bushing_T10_algo"
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-nodal-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstress(modeldata)
	modeldata["postprocessing"] = FDataDict("file"=>filebase * "-elwise-", "quantity"=>:pressure, "component"=>collect(1:1), "nodevalmethod"=>:averaging, "reportat"=>:extraptrend)
	modeldata = AlgoDeforLinearModule.exportstresselementwise(modeldata)

    true
end # bushing_t10_algo_ms

#     true
# end # twisted_beam_export

function allrun()
    println("#####################################################")
    println("# bushing_h8_full ")
    bushing_h8_full()
    println("#####################################################")
    println("# bushing_h8_algo_full ")
    bushing_h8_algo_full()
    println("#####################################################")
    println("# bushing_h8_algo_ms ")
    bushing_h8_algo_ms()
    println("#####################################################")
    println("# bushing_h8u_algo_ms ")
    bushing_h8u_algo_ms()
    println("#####################################################")
    println("# bushing_t10_algo_ms")
    bushing_t10_algo_ms()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing
    
