module bend_hex_spectrum_examples
using FinEtools
using FinEtools.MeshExportModule.CSV: savecsv
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtoolsDeforLinear.AlgoDeforLinearModule: linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky, eigen
using PGFPlotsX

# Isotropic material
E = 1.0;
nu = 0.3
aspect = 100.0
# julia> rand(8,3) .* 2 .- 1.0           
# xyzperturbation = [         
# 0.0767656  -0.983206    -0.14343     
# 0.45767     0.981479     0.450997    
# -0.295854    0.542922     0.321333    
# -0.85204    -0.97824     -0.772874    
# -0.0539756   0.994907     0.822798    
# 0.447173    0.528742     0.0445352   
# -0.468417    0.00673427   0.953151    
# -0.898513   -0.915871     0.933237  ] ./ 15.0
# Lambda1 = vec([-1.0 1.0 1.0 -1.0 -1.0 1.0 1.0 -1.0])
# Lambda2 = vec([-1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0])
# Gamma1 = vec([1.0 1.0 -1.0 -1.0 -1.0 -1.0 1.0 1.0])
# xyzperturbation = zeros((8,3))  
# xyzperturbation[:, 1] .= Lambda2 ./ 10.0
# xyzperturbation[:, 2] .= Lambda1 ./ 10.0
# xyzperturbation[:, 3] .= Gamma1 ./ 20.0

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

function bend_hex_spectrum_full()
   	fens,fes = mesh()
    fens.xyz[:, 1] .*= aspect
    
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
	applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    vol = integratefunction(femm, geom, x -> 1.0, 3)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))
    
    # File =  "bend_hex_spectrum_full.vtk"
    # vectors = [("ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)", deepcopy(scattersysvec!(u, D.vectors[:,idx]).values)) for idx in 1:length(D.values)] 
    # vtkexportmesh(File, fens, fes;  vectors=vectors)
    # @async run(`"paraview.exe" $File`)

    savecsv("bend_hex_spectrum_full-aspect=$(aspect).csv", eigenvalues = vec(D.values))
    # @pgf _a = SemiLogYAxis({
    #     xlabel = "Mode [ND]",
    #     ylabel = "Generalized stiffness [N/m]",
    #     grid="major",
    #     legend_pos  = "south east",
    #     title = "Hexahedron spectrum, \\aspect=$(aspect)"
    # },
    # Plot({"red", mark="triangle"}, Table([:x => vec(7:size(K, 1)), :y => vec(D.values[7:end])])), LegendEntry("FEA"))
    # display(_a)

    true

end # bend_hex_spectrum_full

function bend_hex_spectrum_underintegrated()
	fens,fes = mesh()
	fens.xyz[:, 1] .*= aspect

	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, E, nu)

	geom = NodalField(fens.xyz)
	u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
	applyebc!(u)
	numberdofs!(u)

	femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)

	vol = integratefunction(femm, geom, x -> 1.0, 3)

	associategeometry!(femm, geom)
	K = stiffness(femm, geom, u)

	D = eigen(Matrix(K))

	# File =  "bend_hex_spectrum_full.vtk"
	# vectors = [("ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)", deepcopy(scattersysvec!(u, D.vectors[:,idx]).values)) for idx in 1:length(D.values)] 
	# vtkexportmesh(File, fens, fes;  vectors=vectors)
	# @async run(`"paraview.exe" $File`)

	savecsv("bend_hex_spectrum_underintegrated-aspect=$(aspect).csv", eigenvalues = vec(D.values))
	
	true

end # bend_hex_spectrum_underintegrated

function bend_hex_spectrum_ms()
   	fens,fes = mesh()
    fens.xyz[:, 1] .*= aspect
    
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
	applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    vol = integratefunction(femm, geom, x -> 1.0, 3)

    associategeometry!(femm, geom)
    
    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))
    
    savecsv("bend_hex_spectrum_ms-aspect=$(aspect).csv", eigenvalues = vec(D.values))
    
    # File =  "bend_hex_spectrum_ms.vtk"
    # vectors = [("ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)", deepcopy(scattersysvec!(u, D.vectors[:,idx]).values)) for idx in 1:length(D.values)] 
    # vtkexportmesh(File, fens, fes;  vectors=vectors)
    # @async run(`"paraview.exe" $File`)

    true

end # bend_hex_spectrum_ms

function allrun()
	println("#####################################################")
    println("# bend_hex_spectrum_full ")
    bend_hex_spectrum_full()
    println("#####################################################")
    println("# bend_hex_spectrum_underintegrated ")
    bend_hex_spectrum_underintegrated()
    println("#####################################################")
    println("# bend_hex_spectrum_ms ")
    bend_hex_spectrum_ms()
    return true
end # function allrun

end # module bend_hex_spectrum_examples
