module solver_examples
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy, solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule:
    linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra
using SparseArrays
using LinearSolve
using IncompleteLU
using Printf
using SymRCM
using UnicodePlots
using PlotlyJS
# using Infiltrator
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

function example(n = 10, precond = :ilu, alg = :cg)
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag)
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
    println("nfreedofs(u) = $(nfreedofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    F_f = vector_blocked_f(F, nfreedofs(u))
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    K_ff = matrix_blocked_ff(K, nfreedofs(u))
    K = nothing
    println("Stiffness: number of non zeros = $(nnz(K_ff)) [ND]")
    println("Sparsity = $(nnz(K_ff)/size(K_ff, 1)/size(K_ff, 2)) [ND]")
    # display(spy(K_ff, canvas = DotCanvas))
    
    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)

    if precond == :ilu
        mK_ffd = mean(diag(K_ff))
        PRECND = ilu(K_ff, τ = mK_ffd / 100.0)
    elseif precond == :kdiag
        PRECND = Diagonal(diag(K_ff))
    end

    if alg == :cg
        ALG = KrylovJL_CG
    elseif alg == :gmres
        ALG = KrylovJL_GMRES
    end

    @time prob = LinearProblem(K_ff, F_f)
    @time sol = solve(prob, ALG(), Pl=PRECND)
    scattersysvec!(u, sol.u[:])

    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "example-n=$(n).vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    @async run(`"paraview.exe" $File`)
 
    true
end # example

function allrun(args...)
    println("#####################################################")
    println("# example ")
    example(args...)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
