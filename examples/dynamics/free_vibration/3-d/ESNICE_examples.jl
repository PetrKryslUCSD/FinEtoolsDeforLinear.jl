module ESNICE_examples
using Statistics
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using PGFPlotsX
using Test
using StatsBase

function ESNICE_energies()
    E = 1e6 * phun("PA")
    nu = 0.0
    L = 2 * phun("M")
    hs = L * collect(10 .^ range(-4.0, stop = 0.0, length = 10))
    mag = 0.001

    rs = Float64[]
    PEs = Float64[]
    APEs = Float64[]
    for h in hs
        xs = collect(linearspace(0.0, L, 2))
        ys = collect(linearspace(0.0, h, 2))
        zs = collect(linearspace(0.0, h, 2))
        global fens
        global fes
        fens, fes = T4blockx(xs, ys, zs, :a)
        fens.xyz[:, 3] .-= h / 2
        MR = DeforModelRed3D
        global geom = NodalField(fens.xyz)
        global u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
        numberdofs!(u)

        l1 = selectnode(fens; plane = [1.0 0.0 0.0 0.0], thickness = h / 1000)
        for i in l1
            x, y, z = geom.values[i, :]
            u.values[i, 1] = z * mag
        end
        l1 = selectnode(fens; plane = [1.0 0.0 0.0 L], thickness = h / 1000)
        for i in l1
            x, y, z = geom.values[i, :]
            u.values[i, 1] = -z * mag
        end
        @show h, u.values
        material = MatDeforElastIso(MR, E, nu)
        global femm = FEMMDeforLinear(MR, IntegDomain(fes, TetRule(1)), material)
        associategeometry!(femm, geom)

        K = stiffness(femm, geom, u)
        U = gathersysvec(u)
        PE = dot(U, K * U) / 2.0
        push!(PEs, PE)
        I = h * h^3 / 12
        APE = 2 * E * I * mag^2 / L
        push!(APEs, APE)

        # ars = []
        # for i = 1:count(fes)
        #     res  = FinEtoolsDeforLinear.FEMMDeforLinearESNICEModule.aspectratio(geom.values[collect(fes.conn[i]), :])
        #     ar = sort([res[1], res[2], res[3], res[4]])
        #     push!(ars, mean([ar[2:3]...]))
        # end
        # @show h/L, ars
        # push!(rs, minimum(ars))
        push!(rs, h / L)
    end

    @show rPE = PEs ./ APEs

    # Least-squares fit
    A = hcat([-log10(r) for r in rs], [-1 for r in rPE])
    b = [log10(r - 1) for r in rPE]
    p = A \ b
    @show p
    a = p[1]
    b = 10^p[2]
    @show a, b

    @pgf a = Axis({
            xlabel = "Aspect ratio",
            ylabel = "Relative Potential Energy",
            grid = "major",
            legend_pos = "north east",
        },
        Plot({color = "red"}, Table([:x => log10.(vec(rs)), :y => log10.(vec(rPE))])),
        Plot({"only marks", mark = "x"},
            Table([
                :x => log10.(vec(rs)),
                :y => log10.(vec([1 / (b * r^a) + 1 for r in rs])),
            ])))
    display(a)

    # fld = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)
    # File =  "mt4energy2.vtk"
    # vtkexportmesh(File, fens, fes; scalars=[("sigmax", fld.values)], vectors=[("u", u.values)])
    # @async run(`"paraview.exe" $File`)
    true
end

function ESNICE_vibration()
    E = 70000 * phun("MPa")
    nu = 0.33
    rho = 2700 * phun("KG/M^3")
    radius = 0.5 * phun("ft")
    neigvs = 20                   # how many eigenvalues
    OmegaShift = (10.0 * 2 * pi)^2

    MR = DeforModelRed3D
    output = import_ABAQUS(joinpath(@__DIR__, "alum_cyl.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm, geom)
    @show minimum(vec(femm.nphis)), maximum(vec(femm.nphis))
    @pgf a = Axis({
            xlabel = "Entity",
            ylabel = "Stabilization factor",
            grid = "major",
            legend_pos = "north east",
        },
        Plot({mark = "circle"}, Table([:x => vec(1:count(fes)), :y => vec(femm.ephis)])))
    display(a)
    K = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d, v, nconv = eigs(K + OmegaShift * M,
        M;
        nev = neigvs,
        which = :SM,
        explicittransform = :none)
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    println("Eigenvalues: $fs [Hz]")

    vectors = []
    for i in 7:length(fs)
        scattersysvec!(u, v[:, i])
        push!(vectors, ("Mode_$i", deepcopy(u.values)))
    end
    File = "alum_cyl_mode_shapes.vtk"
    vtkexportmesh(File,
        connasarray(fes),
        fens.xyz,
        FinEtools.MeshExportModule.VTK.T4;
        vectors = vectors)
    @async run(`"paraview.exe" $File`)

    true
end # function

function allrun()
    println("#####################################################")
    println("# ESNICE_energies ")
    ESNICE_energies()
    println("#####################################################")
    println("# ESNICE_vibration ")
    ESNICE_vibration()
    return true
end # function allrun

@info "All examples may be executed with "
println("using .$(@__MODULE__); $(@__MODULE__).allrun()")

end # module 
nothing
