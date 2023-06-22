module multimaterial_nas_examples
using FinEtools
using FinEtools.MeshExportModule: MESH
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtoolsDeforLinear
using LinearAlgebra
using SparseArrays
using Arpack

function multimaterial_nas()
    println("""
    Vibration modes of block composed of two materials. NASTRAN input file.
    """)

    materials = Dict(
        1=>(name = "steel", E = 205000*phun("MPa"), nu = 0.3, rho = 7850*phun("KG*M^-3")), 
        2=>(name = "aluminum", E = 70000*phun("MPa"), nu = 0.34, rho = 2700*phun("KG*M^-3"))
        )
    OmegaShift = (2*pi*100) ^ 2; # to resolve rigid body modes
    neigvs = 20;

    # No need to change anything below this line ##########
    MR = DeforModelRed3D
    output = import_NASTRAN(joinpath(@__DIR__, "twoblocks.nas"))
    fens, fesets, pids = output["fens"], output["fesets"], output["property_ids"] 

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

    numberdofs!(u)

    K = spzeros(nalldofs(u), nalldofs(u))
    M = spzeros(nalldofs(u), nalldofs(u))
    allfes = nothing
    for i in 1:length(fesets)
        pid = pids[i]
        @show E, nu, rho = materials[pid].E, materials[pid].nu, materials[pid].rho
        material = MatDeforElastIso(MR, rho, E, nu, 0.0)
        femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fesets[i], NodalSimplexRule(3)), material)
        femm = associategeometry!(femm, geom)
        K += stiffness(femm, geom, u)
        M += mass(femm, geom, u)
        if allfes == nothing
            allfes = fesets[i]
        else
            allfes = cat(allfes, fesets[i])
        end
    end
    
    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.

    if true
        d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
        d = d .- OmegaShift;
        fs = real(sqrt.(complex(d)))/(2*pi)
        println("Eigenvalues: $fs [Hz]")
        mode = 7
        scattersysvec!(u, v[:,mode])
        File =  "multimaterial_nas.vtk"
        vtkexportmesh(File, fens, allfes; vectors=[("mode$mode", u.values)])
        @async run(`"paraview.exe" $File`)
    end

    # Extract the boundary
    bfes = meshboundary(allfes)
    bconn = connasarray(bfes)

    MESH.write_MESH("multimaterial_nas.mesh", fens, bfes) 

    true

end # multimaterial

function allrun()
    println("#####################################################")
    println("# multimaterial_nas ")
    multimaterial_nas()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")


end # module 
nothing
