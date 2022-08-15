module hex_ip_constr_examples
using FinEtools
using FinEtools.MeshExportModule.CSV: savecsv
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtoolsDeforLinear.AlgoDeforLinearModule: linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky, eigen, svd
using PGFPlotsX

mesh() = (FinEtools.FENodeSetModule.FENodeSet([
	0.0 0.0 0.0; 
	1.0 0.0 0.0; 
	1.0 1.0 0.0; 	
	0.0 1.0 0.0; 
	0.0 0.0 1.0;
	1.0 0.0 1.0; 
	1.0 1.0 1.0
	0.0 1.0 1.0; 
	]), FinEtools.FESetModule.FESetH8(reshape([1, 2, 3, 4, 5, 6, 7, 8], 1, 8)))

function hex_ip_constr_par()
	fens,fes = mesh()

	# integration_rule = TrapezoidalRule(3)
	integration_rule = GaussRule(3, 2)
	pc = integration_rule.param_coords;
	w  =  integration_rule.weights
	npts = integration_rule.npts;

	Q = fill(0.0, npts, 24)
	for j in 1:npts
		gradNpar = bfundpar(fes, vec(pc[j,:]))
		k = 1
		for i in 1:size(gradNpar, 1)
			Q[j, k:(k+2)] .= gradNpar[i, :]
			k = k + 3
		end
	end

	decomp = svd(Q)
	@show decomp.S
	
		true
end # hex_ip_constr_par

function hex_ip_constr_xyz()
	xyzperturbation = [         
	0.0767656  -0.983206    -0.14343     
	0.45767     0.981479     0.450997    
	-0.295854    0.542922     0.321333    
	-0.85204    -0.97824     -0.772874    
	-0.0539756   0.994907     0.822798    
	0.447173    0.528742     0.0445352   
	-0.468417    0.00673427   0.953151    
	-0.898513   -0.915871     0.933237  ] ./ 5.0
	fens,fes = mesh()
	fens.xyz .+= xyzperturbation

	# integration_rule = TrapezoidalRule(3)
	integration_rule = GaussRule(3, 2)
	pc = integration_rule.param_coords;
	w  =  integration_rule.weights
	npts = integration_rule.npts;

	Q = fill(0.0, npts, 24)
	for j in 1:npts
		gradNpar = bfundpar(fes, vec(pc[j,:]))
		J = (transpose(fens.xyz) * gradNpar)
		gradN = gradNpar / J
		k = 1
		for i in 1:size(gradN, 1)
			Q[j, k:(k+2)] .= gradN[i, :]
			k = k + 3
		end
	end

	decomp = svd(Q)
	@show decomp.S
	
	true
end # hex_ip_constr_xyz

function allrun()
	println("#####################################################")
	println("# hex_ip_constr_par ")
	hex_ip_constr_par()
	println("#####################################################")
	println("# hex_ip_constr_xyz ")
	hex_ip_constr_xyz()
	return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing