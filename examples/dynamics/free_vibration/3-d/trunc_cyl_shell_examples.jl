module trunc_cyl_shell_examples
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using LinearAlgebra
using Arpack
using SubSIt

function trunc_cyl_shell()
    println("""
    Vibration modes of truncated cylindrical shell.
    """)

    # t0 = time()

    E = 205000 * phun("MPa")# Young's modulus
    nu = 0.3# Poisson ratio
    rho = 7850 * phun("KG*M^-3")# mass density
    OmegaShift = (2 * pi * 100)^2 # to resolve rigid body modes
    h = 0.05 * phun("M")
    l = 10 * h
    Rmed = h / 0.2
    psi = 0    # Cylinder
    nh = 5
    nl = 12
    nc = 40
    tolerance = h / nh / 100
    neigvs = 20

    MR = DeforModelRed3D
    fens, fes = H8block(h, l, 2 * pi, nh, nl, nc)
    # Shape into a cylinder
    R = zeros(3, 3)
    for i = 1:count(fens)
        x, y, z = fens.xyz[i, :]
        rotmat3!(R, [0, z, 0])
        Q = [
            cos(psi * pi / 180) sin(psi * pi / 180) 0
            -sin(psi * pi / 180) cos(psi * pi / 180) 0
            0 0 1
        ]
        fens.xyz[i, :] = reshape([x + Rmed - h / 2, y - l / 2, 0], 1, 3) * Q * R
    end
    candidates = selectnode(fens, plane = [0.0 0.0 1.0 0.0], thickness = h / 1000)
    fens, fes = mergenodes(fens, fes, tolerance, candidates)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 3)), material)
    M = mass(femm, geom, u)
    K_ff = matrix_blocked(K, nfreedofs(u), nfreedofs(u))[:ff]
    M_ff = matrix_blocked(M, nfreedofs(u), nfreedofs(u))[:ff]

    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.

    if true
        d, v, nev, nconv = eigs(
            Symmetric(K_ff + OmegaShift * M_ff),
            Symmetric(M_ff);
            nev = neigvs,
            which = :SM,
            explicittransform = :none,
        )
        d = d .- OmegaShift
        fs = real(sqrt.(complex(d))) / (2 * pi)
        println("Eigenvalues: $fs [Hz]")
        mode = 7
        scattersysvec!(u, v[:, mode])
        File = "trunc_cyl_shell.vtk"
        vtkexportmesh(File, fens, fes; vectors = [("mode$mode", u.values)])
        @async run(`"paraview.exe" $File`)
    end

    if true
        solver = SubSIt.ssit
        v0 = rand(size(K_ff, 1), 2 * neigvs)
        tol = 1.0e-2
        maxiter = 20
        lamb, v, nconv, niter, lamberr = solver(
            K_ff + OmegaShift .* M_ff,
            M_ff;
            nev = neigvs,
            X = v0,
            tol = tol,
            maxiter = maxiter,
        )
        if nconv < neigvs
            println("NOT converged")
        end
        lamb = lamb .- OmegaShift
        fs = real(sqrt.(complex(lamb))) / (2 * pi)
        println("Eigenvalues: $fs [Hz]")
        println("Eigenvalue errors: $lamberr [ND]")
        mode = 7
        scattersysvec!(u, v[:, mode])
        File = "trunc_cyl_shell.vtk"
        vtkexportmesh(File, fens, fes; vectors = [("mode$mode", u.values)])
        @async run(`"paraview.exe" $File`)
    end

    true
end # trunc_cyl_shell

function trunc_cyl_shell_nas()
    println("""
    Vibration modes of truncated cylindrical shell. NASTRAN input file.
    """)

    # t0 = time()

    E = 205000 * phun("MPa")# Young's modulus
    nu = 0.3# Poisson ratio
    rho = 7850 * phun("KG*M^-3")# mass density
    OmegaShift = (2 * pi * 100)^2 # to resolve rigid body modes
    h = 0.05 * phun("M")
    l = 10 * h
    Rmed = h / 0.2
    psi = 0    # Cylinder
    nh = 5
    nl = 12
    nc = 40
    tolerance = h / nh / 100
    neigvs = 20

    MR = DeforModelRed3D
    output = import_NASTRAN(joinpath(@__DIR__, "trunc_cyl_shell_2.nas"))
    fens, fes = output["fens"], output["fesets"][1]

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    femm = associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    M = mass(femm, geom, u)

    # eigs returns the nev requested eigenvalues in d, the corresponding Ritz vectors
    # v (only if ritzvec=true), the number of converged eigenvalues nconv, the number
    # of iterations niter and the number of matrix vector multiplications nmult, as
    # well as the final residual vector resid.

    if true
        d, v, nev, nconv = eigs(K + OmegaShift * M, M; nev = neigvs, which = :SM)
        d = d .- OmegaShift
        fs = real(sqrt.(complex(d))) / (2 * pi)
        println("Eigenvalues: $fs [Hz]")
        mode = 7
        scattersysvec!(u, v[:, mode])
        File = "trunc_cyl_shell_nas.vtk"
        vtkexportmesh(File, fens, fes; vectors = [("mode$mode", u.values)])
        @async run(`"paraview.exe" $File`)
    end

    true
end # trunc_cyl_shell

function allrun()
    println("#####################################################")
    println("# trunc_cyl_shell ")
    trunc_cyl_shell()
    println("#####################################################")
    println("# trunc_cyl_shell_nas ")
    trunc_cyl_shell_nas()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
