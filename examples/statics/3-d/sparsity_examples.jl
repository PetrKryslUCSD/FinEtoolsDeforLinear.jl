module sparsity_examples
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
using DataStructures
using AddToSparse

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

function __collect_unique_node_neighbors(ellist, conn, npe)
    totn = length(ellist) * npe
    nodes = fill(zero(eltype(conn[1])), totn)
    p = 1
    for i in ellist
        for k in conn[i]
            nodes[p] = k 
            p += 1
        end
    end
    sort!(nodes)
    unique!(nodes)
    return nodes
end

function _unique_nodes(n2e, conn)
    npe = length(conn[1])
    unique_nodes = fill(Vector{eltype(n2e.map[1])}(), length(n2e.map))
    Base.Threads.@threads for i in 1:length(n2e.map) # run this in PARALLEL
        unique_nodes[i] = __collect_unique_node_neighbors(n2e.map[i], conn, npe)
    end
    return unique_nodes
end

function _populate_dofs(n, n2n, dofnums, start, dofs)
    nd = size(dofnums, 2)
    totd = length(n2n[n]) * nd
    _dofs = fill(zero(eltype(dofs)), totd)
    p = 1
    for k in n2n[n]
        for d in axes(dofnums, 2)
            _dofs[p] = dofnums[k, d]
            p += 1
        end
    end
    sort!(_dofs)
    for d in axes(dofnums, 2)
        j = dofnums[n, d]
        s = start[j]
        p = 0
        for m in eachindex(_dofs)
            dofs[s+p] = _dofs[m]
            p += 1
        end
    end
    return nothing
end

function _prepare_start_dofs(IT, n2n, dofnums)
    nd = size(dofnums, 2)
    total_dofs = length(n2n) * nd
    lengths = Vector{IT}(undef, total_dofs+1)
    for k in eachindex(n2n)
        kl = length(n2n[k]) * nd
        for d in axes(dofnums, 2)
            j = dofnums[k, d]
            lengths[j] = kl
        end
    end
    lengths[end] = 0
    # Now we start overwriting the lengths array with the starts
    start = lengths
    sumlen = 0
    len = start[1]
    sumlen += len
    start[1] = 1
    plen = len
    for k in 2:total_dofs
        len = start[k]
        sumlen += len
        start[k] = start[k-1] + plen 
        plen = len
    end
    start[end] = sumlen+1
    dofs = Vector{IT}(undef, sumlen)
    return start, dofs
end

function example(n=10; precond=:ilu, alg=:cg, other...)
    elementtag = "H8"
    println("""
    Stubby corbel example. Element: $(elementtag), n=$(n)
    """)

    fens, fes = H8block(W, L, H, n, 2 * n, 2 * n)
    println("Number of elements: $(count(fes))")
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing=true, direction=[0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing=true, direction=[0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing=true, direction=[1.0 0.0 0.0])

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

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 3)), material)

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
    associategeometry!(femm, geom)
    ass = SysmatAssemblerSparse(0.0)
    setnomatrixresult(ass, false)
    K = stiffness(femm, ass, geom, u)
    I, J, V = deepcopy(ass._rowbuffer), deepcopy(ass._colbuffer), deepcopy(ass._matbuffer)
    @time K = sparse(I, J, V, nalldofs(u), nalldofs(u))
    # @show K.colptr
    # @show K.rowval

    IT = eltype(u.dofnums)
    @time n2e = FENodeToFEMap(fes.conn, count(fens))
    @time n2n = _unique_nodes(n2e, fes.conn)
    @time start, dofs = _prepare_start_dofs(IT, n2n, u.dofnums)
    
    @time Base.Threads.@threads for n in 1:count(fens)
        _populate_dofs(n, n2n, u.dofnums, start, dofs)
    end
    # @show start
    # @show dofs
    S = SparseMatrixCSC(size(K)..., start, dofs, fill(0.0, length(dofs)))
    AddToSparse.addtosparse(S, I, J, V)
    @show norm(K - S) / norm(K)
    
    true
end # example

function allrun(n=3; args...)
    println("#####################################################")
    println("# example ")
    example(n; args...)
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
