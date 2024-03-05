module single_hex_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.AlgoBaseModule: evalconvergencestudy
using FinEtoolsDeforLinear.AlgoDeforLinearModule:
    linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra: Symmetric, cholesky, eigen, norm

# Isotropic material
E = 1.0
nu = 0.3
# julia> rand(8,3) .* 2 .- 1.0           
xyzperturbation =
    [
        0.0767656 -0.983206 -0.14343
        0.45767 0.981479 0.450997
        -0.295854 0.542922 0.321333
        -0.85204 -0.97824 -0.772874
        -0.0539756 0.994907 0.822798
        0.447173 0.528742 0.0445352
        -0.468417 0.00673427 0.953151
        -0.898513 -0.915871 0.933237
    ] ./ 5.0
# Lambda1 = vec([-1.0 1.0 1.0 -1.0 -1.0 1.0 1.0 -1.0])
# Lambda2 = vec([-1.0 -1.0 1.0 1.0 -1.0 -1.0 1.0 1.0])
# Gamma1 = vec([1.0 1.0 -1.0 -1.0 -1.0 -1.0 1.0 1.0])
# xyzperturbation = zeros((8,3))  
# xyzperturbation[:, 1] .= Lambda2 ./ 10.0
# xyzperturbation[:, 2] .= Lambda1 ./ 10.0
# xyzperturbation[:, 3] .= Gamma1 ./ 20.0

function mesh(alpha = 0)
    (
        FinEtools.FENodeSetModule.FENodeSet(
            [
                0.0 0.0 0.0
                1.0 0.0 0.0
                1.0 1.0 0.0
                0.0 1.0 0.0
                0.0 0.0 1.0
                1.0 0.0 1.0
                1.0 1.0 1.0
                0.0 1.0 1.0
            ] + alpha * xyzperturbation,
        ),
        FinEtools.FESetModule.FESetH8(reshape([1, 2, 3, 4, 5, 6, 7, 8], 1, 8)),
    )
end

function single_hex_perfect_cube()
    fens, fes = mesh()

    File = "single_hex_perfect_cube.vtk"
    vtkexportmesh(File, fens, fes)
    @async run(`"paraview.exe" $File`)

    true
end # single_hex_perfect_cube

function single_hex_full()
    fens, fes = mesh()
    fens.xyz += xyzperturbation

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    @show vol = integratefunction(femm, geom, x -> 1.0; m = 3)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))

    File = "single_hex_full.vtk"
    vectors = [
        (
            "ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)",
            deepcopy(scattersysvec!(u, D.vectors[:, idx]).values),
        ) for idx = 1:length(D.values)
    ]
    vtkexportmesh(File, fens, fes; vectors = vectors)
    @async run(`"paraview.exe" $File`)

    true
end # single_hex_full

function single_hex_underintegrated()
    fens, fes = mesh()
    fens.xyz += xyzperturbation

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 1)), material)

    @show vol = integratefunction(femm, geom, x -> 1.0; m = 3)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))

    File = "single_hex_underintegrated.vtk"
    vectors = [
        (
            "ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)",
            deepcopy(scattersysvec!(u, D.vectors[:, idx]).values),
        ) for idx = 1:length(D.values)
    ]
    vtkexportmesh(File, fens, fes; vectors = vectors)
    @async run(`"paraview.exe" $File`)

    true
end # single_hex_underintegrated

function single_hex_ms()
    fens, fes = mesh()
    fens.xyz += xyzperturbation

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    @show vol = integratefunction(femm, geom, x -> 1.0; m = 3)

    associategeometry!(femm, geom)
    fill!(femm.phis, 0.0)
    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))

    File = "single_hex_ms.vtk"
    vectors = [
        (
            "ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)",
            deepcopy(scattersysvec!(u, D.vectors[:, idx]).values),
        ) for idx = 1:length(D.values)
    ]
    vtkexportmesh(File, fens, fes; vectors = vectors)
    @async run(`"paraview.exe" $File`)

    true
end # single_hex_ms

function single_hex_full_scaling()
    fens, fes = mesh(1.0)
    @show fens.xyz

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    # @show vol = integratefunction(femm, geom, x -> 1.0; m=3)

    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K1 = Matrix(K)

    fens, fes = mesh(1.0)
    fens.xyz .*= 10
    geom = NodalField(fens.xyz)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K10 = Matrix(K)

    @show norm(K1 - K10) / norm(K1)

    fens, fes = mesh(1.0)
    fens.xyz .*= 1
    geom = NodalField(fens.xyz)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K1b = Matrix(K)

    @show norm(K1 - K1b) / norm(K1)

    fens, fes = mesh(1.0)
    fens.xyz .*= 2
    geom = NodalField(fens.xyz)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K2 = Matrix(K)

    @show norm(K1 - K2) / norm(K1)

    fens, fes = mesh(1.0)
    fens.xyz .*= 200
    geom = NodalField(fens.xyz)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K200 = Matrix(K)

    @show norm(K1 - K200) / norm(K1)

    true
end # single_hex_full

function allrun()
    println("#####################################################")
    println("# single_hex_full_scaling ")
    single_hex_full_scaling()
    println("#####################################################")
    println("# single_hex_perfect_cube ")
    single_hex_perfect_cube()
    println("#####################################################")
    println("# single_hex_full ")
    single_hex_full()
    println("#####################################################")
    println("# single_hex_underintegrated ")
    single_hex_underintegrated()
    println("#####################################################")
    println("# single_hex_ms ")
    single_hex_ms()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
