module single_hex_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtoolsDeforLinear.AlgoDeforLinearModule: linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky, eigen

# Isotropic material
E = 1.0;
nu = 0.3
a, b, c = 1.0, 1.0, 1.0
# julia> rand(8,3) .* 2 .- 1.0           
perturbation = [         
  0.0767656  -0.983206    -0.14343     
  0.45767     0.981479     0.450997    
 -0.295854    0.542922     0.321333    
 -0.85204    -0.97824     -0.772874    
 -0.0539756   0.994907     0.822798    
  0.447173    0.528742     0.0445352   
 -0.468417    0.00673427   0.953151    
 -0.898513   -0.915871     0.933237  ] ./ 10.0

function single_hex_full()
   	fens,fes = H8block(a, b, c, 1, 1, 1)
    fens.xyz +=  perturbation
    
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
	applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    @show vol = integratefunction(femm, geom, x -> 1.0, 3)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))
    
    File =  "single_hex_full.vtk"
    vectors = [("ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)", deepcopy(scattersysvec!(u, D.vectors[:,idx]).values)) for idx in 1:length(D.values)] 
    vtkexportmesh(File, fens, fes;  vectors=vectors)
    @async run(`"paraview.exe" $File`)

    true

end # single_hex_full

function single_hex_underintegrated()
   	fens,fes = H8block(a, b, c, 1, 1, 1)
    fens.xyz +=  perturbation
    
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1), 3)) # displacement field
	applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)

    @show vol = integratefunction(femm, geom, x -> 1.0, 3)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))
    
    File =  "single_hex_underintegrated.vtk"
    vectors = [("ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)", deepcopy(scattersysvec!(u, D.vectors[:,idx]).values)) for idx in 1:length(D.values)] 
    vtkexportmesh(File, fens, fes;  vectors=vectors)
    @async run(`"paraview.exe" $File`)

    true

end # single_hex_underintegrated

function allrun()
    println("#####################################################")
    println("# single_hex_full ")
    single_hex_full()
    println("#####################################################")
    println("# single_hex_underintegrated ")
    single_hex_underintegrated()
    return true
end # function allrun

end # module single_hex_examples
