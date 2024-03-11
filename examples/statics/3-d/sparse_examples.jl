module sparse_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy, solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule:
    linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra
using SparseArrays
using LinearSolve
using SciMLOperators
using IncompleteLU
using Printf
using SymRCM
using SimplySparse
using SparseMatricesCSR
using ThreadedSparseCSR
using UnicodePlots
using PlotlyJS
using Random
using DataDrop

# Isotropic material
E = 1000.0
nu = 0.4999 # Taylor data: nearly incompressible material
nu = 0.3 # Compressible material
W = 25.0
H = 50.0
L = 50.0
htol = minimum([L, H, W]) / 1000
uzex = -0.16
magn = 0.2 * (-12.6) / 4
Force = magn * W * H * 2
CTE = 0.0
n = 5 #

function getfrcL!(forceout, XYZ, tangents, feid, qpid)
    copyto!(forceout, [0.0; 0.0; magn])
end

function example(n = 10; precond = :ilu, alg = :cg, other...)
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag), n=$(n)
    """)

    fens, fes = H8block(W, L, H, n, 2 * n, 2 * n)
    println("Number of elements: $(count(fes))")
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid, qpid)
        copyto!(csmatout, csmat)
    end

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    # Renumber the nodes
    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    C = connectionmatrix(femm, count(fens))
    perm = symrcm(C)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u, perm)
    # numberdofs!(u)
    println("nalldofs(u) = $(nalldofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    F_f = vector_blocked_f(F, nfreedofs(u))
    associategeometry!(femm, geom)
    ass = SysmatAssemblerSparse(0.0)
    setnomatrixresult(ass, false)
    K = stiffness(femm, ass, geom, u)
    I, J, V = deepcopy(ass._rowbuffer), deepcopy(ass._colbuffer), deepcopy(ass._matbuffer)
    println("Stiffness: number of non zeros = $(length(I)) [ND]")
    
    @time S = sparse(I, J, V, nalldofs(u), nalldofs(u))
    # @time S = SimplySparse.sparse(I, J, V, nalldofs(u), nalldofs(u))
    I, J, V = deepcopy(ass._rowbuffer), deepcopy(ass._colbuffer), deepcopy(ass._matbuffer)
    @time S = SimplySparse.par_sparse(I, J, V, nalldofs(u), nalldofs(u))
    true
end # example

function allrun(n = 30; args...)
    println("#####################################################")
    println("# example ")
    example(n; args...)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
