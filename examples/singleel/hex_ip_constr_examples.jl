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

function hex_ip_constr_full()
	fens,fes = mesh()

	# MR = DeforModelRed3D
	# material = MatDeforElastIso(MR, E, nu)

	# geom = NodalField(fens.xyz)
	# u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
	# applyebc!(u)
	# numberdofs!(u)

	integration_rule = TrapezoidalRule(3)
	pc = integration_rule.param_coords;
	w  =  integration_rule.weights ;
	npts = integration_rule.npts;

	C = fill(0.0, npts, 24)
	for j in 1:npts
		gradN = bfundpar(fes, vec(pc[j,:]))
		k = 1
		for i in 1:size(gradN, 1)
			C[j, k:(k+2)] .= gradN[i, :]
			k = k + 3
		end
	end

	decomp = svd(C)
	@show decomp.S
	# femm = FEMMDeforLinear(MR, IntegDomain(fes, ), material)

	# vol = integratefunction(femm, geom, x -> 1.0, 3)

	# associategeometry!(femm, geom)
	# K = stiffness(femm, geom, u)

	# D = eigen(Matrix(K))

		# File =  "comp_hex_spectrum_full.vtk"
		# vectors = [("ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)", deepcopy(scattersysvec!(u, D.vectors[:,idx]).values)) for idx in 1:length(D.values)] 
		# vtkexportmesh(File, fens, fes;  vectors=vectors)
		# @async run(`"paraview.exe" $File`)

		# savecsv("comp_hex_spectrum_full-nu=$(nu).csv", eigenvalues = vec(D.values))
		# @pgf _a = SemiLogYAxis({
		#     xlabel = "Mode [ND]",
		#     ylabel = "Generalized stiffness [N/m]",
		#     grid="major",
		#     legend_pos  = "south east",
		#     title = "Hexahedron spectrum, \\nu=$(nu)"
		# },
		# Plot({"red", mark="triangle"}, Table([:x => vec(7:size(K, 1)), :y => vec(D.values[7:end])])), LegendEntry("FEA"))
		# display(_a)

		true
end # hex_ip_constr_full

function allrun()
	println("#####################################################")
	println("# hex_ip_constr_full ")
	hex_ip_constr_full()
	return true
end # function allrun

end # module hex_ip_constr_examples
