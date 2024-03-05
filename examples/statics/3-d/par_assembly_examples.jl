
println("Current folder: $(pwd())")
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.add("FinEtools")
Pkg.add("FinEtoolsDeforLinear")
Pkg.add("ThreadPinning")
Pkg.update()


module par_assembly_examples
using FinEtools

using FinEtoolsDeforLinear

# Isotropic material
E = 1000.0
nu = 0.4999 #Taylor data
# nu=0.3; #Taylor data#.
W = 25.0
H = 50.0
L = 50.0
htol = minimum([L, H, W]) / 1000
CTE = 0.0

function run_example(N = 10, ntasks = 2, do_serial = false)

    fens, fes = H8block(W, L, H, N, N, N)

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(FFlt) : zero(FFlt) for i = 1:3, j = 1:3]

    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    associategeometry!(femm, geom)

    if do_serial
        @time K = stiffness(femm, geom, u)
    end

    function matrixcomputation!(femm, assembler)
        stiffness(femm, assembler, geom, u)
    end

    @time assembler = make_assembler(femms, SysmatAssemblerSparseSymm, u)
    @time start_assembler!(assembler)
    @time assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparseSymm, u)
    @time parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    @time K = make_matrix!(assembler)


    true
end # stubby_corbel_H8_by_hand

end # module 


@show Threads.nthreads()
using ThreadPinning
pinthreads(:cores)

@show N = parse(Int, ARGS[1])
@show ntasks = parse(Int, ARGS[2])

using .Main.par_assembly_examples;
Main.par_assembly_examples.run_example(N, ntasks, true);
using .Main.par_assembly_examples;
Main.par_assembly_examples.run_example(N, ntasks, true);

nothing
