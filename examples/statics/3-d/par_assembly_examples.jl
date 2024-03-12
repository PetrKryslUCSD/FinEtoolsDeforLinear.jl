"""
This is an example that shows scaling with threading assembly that can be
compared with Ferrite.

"""

println("Current folder: $(pwd())")
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.add(url="https://github.com/PetrKryslUCSD/FinEtools.jl.git#main")
Pkg.add(url="https://github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git")
Pkg.add(url="https://github.com/PetrKryslUCSD/ParFEM.jl.git")
Pkg.add("ChunkSplitters")
Pkg.add("SymRCM")
Pkg.add("ThreadPinning")
Pkg.update()
Pkg.status()


module par_assembly_examples
using FinEtools
using FinEtoolsDeforLinear
using ChunkSplitters
using SymRCM
using ParFEM: parallel_make_csc_matrix

function run_example(N = 10, ntasks = 2, do_serial = false)
    E = 1000.0
    nu = 0.4999 #Taylor data
    W = 25.0
    H = 50.0
    L = 50.0
    CTE = 0.0

    fens, fes = H8block(W, L, H, N, N, 10*N)
    #C = connectionmatrix(FEMMBase(IntegDomain(fes, GaussRule(3, 1))), count(fens))
    #ordering = symrcm(C)
    #fens, fes = reordermesh(fens, fes, ordering)

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)
    ir = GaussRule(3, 2)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    if do_serial
        @info "Serial Assembly"
        @time femm = FEMMDeforLinear(MR, IntegDomain(fes, ir), material)
        @time associategeometry!(femm, geom)
        @time K = stiffness(femm, geom, u)
    end

    @info "Parallel Assembly with $(ntasks) Tasks"

    function matrixcomputation!(femm, assembler)
        stiffness(femm, assembler, geom, u)
    end

    femms = FEMMDeforLinear[]
    @time for (ch, j) in chunks(1:count(fes), ntasks)
        femm = FEMMDeforLinear(MR, IntegDomain(subset(fes, ch), ir), material)
        associategeometry!(femm, geom)
        push!(femms, femm)
    end

    @time assembler = make_assembler(femms, SysmatAssemblerSparseSymm, u)
    @time start_assembler!(assembler)
    @time assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparseSymm, u)
    @time parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    @time K = make_matrix!(assembler)

    true
end # run_example

function run_example2(N = 10, ntasks = 2, do_serial = false)
    E = 1000.0
    nu = 0.4999 #Taylor data
    W = 25.0
    H = 50.0
    L = 50.0
    CTE = 0.0

    fens, fes = H8block(W, L, H, N, N, 10*N)
    #C = connectionmatrix(FEMMBase(IntegDomain(fes, GaussRule(3, 1))), count(fens))
    #ordering = symrcm(C)
    #fens, fes = reordermesh(fens, fes, ordering)

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)
    ir = GaussRule(3, 2)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    if do_serial
        @info "Serial Assembly"
        @time femm = FEMMDeforLinear(MR, IntegDomain(fes, ir), material)
        @time associategeometry!(femm, geom)
        @time K = stiffness(femm, geom, u)
    end

    @info "Parallel Assembly with $(ntasks) Tasks"

    function matrixcomputation!(femm, assembler)
        associategeometry!(femm, geom)
        stiffness(femm, assembler, geom, u)
    end

    function createsubdomain(fessubset)
        FEMMDeforLinear(MR, IntegDomain(fessubset, GaussRule(3, 2)), material)
    end

    @time K = parallel_make_csc_matrix(fes, u, createsubdomain, matrixcomputation!, ntasks)
    
    # femms = FEMMDeforLinear[]
    # @time for (ch, j) in chunks(1:count(fes), ntasks)
    #     femm = FEMMDeforLinear(MR, IntegDomain(subset(fes, ch), ir), material)
    #     associategeometry!(femm, geom)
    #     push!(femms, femm)
    # end

    # @time assembler = make_assembler(femms, SysmatAssemblerSparseSymm, u)
    # @time start_assembler!(assembler)
    # @time assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparseSymm, u)
    # @time parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    # @time K = make_matrix!(assembler)

    true
end # run_example2

end # module


@show Threads.nthreads()
using ThreadPinning
pinthreads(:cores)

@show N = parse(Int, ARGS[1])
@show ntasks = parse(Int, ARGS[2])

using .Main.par_assembly_examples;

ex = Main.par_assembly_examples.run_example2
ex(N, 1, true)
ex(N, 1, true)

ex(N, ntasks, false)
ex(N, ntasks, false)

nothing