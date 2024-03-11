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

function _traverse(node, fun)
    if node.leftChild !== nothing
        _traverse(node.leftChild, fun)
    end
    fun(node.data)
    if node.rightChild !== nothing
        _traverse(node.rightChild, fun)
    end
    nothing
end

function traverse(tree, fun)
    _traverse(tree.root, fun)
    nothing
end

mutable struct Collector{IT}
    start::IT
    counter::IT
    v::Vector{IT}
end

function (c::Collector{T})(d) where {T}
    c.v[c.start+c.counter] = d
    c.counter += 1
end

function _col_collector(ellist, conn, dofnums, start, v)
    tree = AVLTree{Int}()
    for i in ellist
        for k in conn[i]
            for d in axes(dofnums, 2)
                push!(tree, dofnums[k, d])
            end
        end
    end
    c = Collector(start, 0, v)
    traverse(tree, c)
    return c.counter
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

    @time femap = FENodeToFEMap(fes.conn, count(fens))
    @time begin
        npe = nodesperelem(fes)
        lens = [(npe*length(femap.map[k])+1)*ndofs(u) for k in 1:length(femap.map)]
        totallen = sum(lens)
        c = fill(zero(Int), totallen)   
        counter = fill(zero(Int), count(fens)) 
        start = similar(counter)
        start[1] = 1
        for k in 2:length(lens)
            start[k] = start[k-1] + lens[k-1]
        end
    end
      
    @time Base.Threads.@threads for n in 1:count(fens)
        counter[n] = _col_collector(femap.map[n], fes.conn, u.dofnums, start[n], c)
    end
    true
end # example

function example0(n=10; precond=:ilu, alg=:cg, other...)
    elementtag = "H8"
    println("""
    Stubby corbel example0. Element: $(elementtag), n=$(n)
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

    @time femap = FENodeToFEMap(fes.conn, count(fens))

    for n in 1:count(fens)
        tree = AVLTree{Int}()
        for i in femap.map[n]
            for k in fes.conn[i]
                for d in 1:ndofs(u)
                    push!(tree, u.dofnums[k, d])
                end
            end
        end
        c = Collector(sizehint!(Vector{Int}(), length(tree)))
        traverse(tree, c)
        @show u.dofnums[n, :], c.v
    end
    true
end # example0

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
