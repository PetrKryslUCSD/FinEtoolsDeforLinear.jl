module test_alum_cyl_mode_examples

using FinEtools
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using LinearAlgebra
using Arpack
using Test

# Mesh alum_cyl.inp
# Abaqus with the Standard isoparametric C3D4 tetrahedron:
C3D4 = [
    1 0
    2 0
    3 6.63586E-005
    4 0.000171053
    5 0.000211299
    6 0.000244378
    7 2564.63
    8 2568.09
    9 2597.26
    10 4094.38
    11 4714.36
    12 4717.19
    13 5181.98
    14 6865.13
    15 6868.17
    16 6962.86
    17 6965.67
    18 7024.97
    19 7029.44
    20 7108.54
]
# Abaqus with the standard quadratic tetrahedron:
C3D10 = [
    1 0
    2 0
    3 0
    4 0.000139365
    5 0.000221551
    6 0.000291805
    7 2546.81
    8 2546.81
    9 2560.69
    10 4100
    11 4693.55
    12 4693.56
    13 5121.57
    14 6841.21
    15 6841.24
    16 6914.22
    17 6914.23
    18 6950.64
    19 6950.66
    20 7000.64
]

neigvs = 24

function alum_cyl_mode_esnice_t4()

    E = 70000 * phun("MPa")
    nu = 0.33
    rho = 2700 * phun("KG/M^3")
    radius = 0.5 * phun("ft")
    OmegaShift = (10.0 * 2 * pi)^2

    MR = DeforModelRed3D
    output = import_ABAQUS(joinpath(@__DIR__, "alum_cyl.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    fens.xyz .*= phun("mm") # The input is provided in SI(mm) units
    fens, fes = T10toT4(fens, fes)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)
    # println("Number of degrees of freedom: $(nfreedofs(u))")

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d, v, nev, nconv = eigs(
        Symmetric(K + OmegaShift * M),
        Symmetric(M);
        nev = neigvs,
        which = :SM,
        explicittransform = :none,
    )
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # println("Eigenvalues: $fs [Hz]")
    # @show     v' * M * v
    @test norm(
        fs - [
            0.00000e+00,
            0.00000e+00,
            0.00000e+00,
            5.54160e-06,
            8.64750e-05,
            1.18749e-04,
            2.49815e+03,
            2.49888e+03,
            2.51331e+03,
            4.08265e+03,
            4.58599e+03,
            4.58642e+03,
            4.98701e+03,
            6.64802e+03,
            6.64848e+03,
            6.67904e+03,
            6.68216e+03,
            6.77789e+03,
            6.78059e+03,
            6.79936e+03,
            6.80400e+03,
            7.38167e+03,
            7.45600e+03,
            7.47771e+03,
        ],
    ) < 0.01


    true

end # alum_cyl_modes

function alum_cyl_mode_esnice_h8()

    E = 70000 * phun("MPa")
    nu = 0.33
    rho = 2700 * phun("KG/M^3")
    radius = 0.5 * phun("ft")
    OmegaShift = (10.0 * 2 * pi)^2
    neigvs = 15

    MR = DeforModelRed3D
    fens, fes = H8cylindern(radius, 4 * radius, 7, 28)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    numberdofs!(u)
    # println("Number of degrees of freedom: $(nfreedofs(u))")

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm =
        FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)
    associategeometry!(femm, geom)
    K = stiffness(femm, geom, u)
    M = mass(femm, geom, u)
    d, v, nev, nconv = eigs(
        Symmetric(K + OmegaShift * M),
        Symmetric(M);
        nev = neigvs,
        which = :SM,
        explicittransform = :none,
    )
    d = d .- OmegaShift
    fs = real(sqrt.(complex(d))) / (2 * pi)
    # println("Eigenvalues: $fs [Hz]")
    # @show     v' * M * v
    @test norm(
        fs - [
            0.00000e+00,
            0.00000e+00,
            0.00000e+00,
            0.00000e+00,
            0.00000e+00,
            0.00000e+00,
            2.55691e+03,
            2.55691e+03,
            2.56007e+03,
            4.09870e+03,
            4.68757e+03,
            4.68757e+03,
            5.10331e+03,
            6.81594e+03,
            6.81594e+03,
        ],
    ) < 0.001 * norm(fs)


    true

end # alum_cyl_modes

alum_cyl_mode_esnice_t4()
alum_cyl_mode_esnice_h8()
# alum_cyl_mode_esnice_t4_ssit()
# alum_cyl_mode_esnice_t4_ssit2()


end # module 
nothing

module imspectrum01
using FinEtools
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using LinearAlgebra
using Arpack
using Test

mesh() = (
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
        ],
    ),
    FinEtools.FESetModule.FESetH8(reshape([1, 2, 3, 4, 5, 6, 7, 8], 1, 8)),
)

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
    ] ./ 15.0

function bend_hex_spectrum_im()
    aspect = 1.0
    fens, fes = mesh()
    fens.xyz += xyzperturbation
    fens.xyz[:, 1] .*= aspect
    E = 1.0
    nu = 0.3
    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, E, nu)

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field
    applyebc!(u)
    numberdofs!(u)

    femm = FEMMDeforLinearIMH8(MR, IntegDomain(fes, GaussRule(3, 2)), material, 9)

    vol = integratefunction(femm, geom, x -> 1.0; m = 3)

    associategeometry!(femm, geom)

    K = stiffness(femm, geom, u)

    D = eigen(Matrix(K))
    @test norm(
        vec(D.values) - [
            -4.58455e-17
            -4.30053e-17
            -2.66920e-17
            -5.87832e-18
            5.79305e-18
            1.75741e-16
            5.75036e-02
            7.27270e-02
            1.13685e-01
            1.17365e-01
            1.23125e-01
            1.24759e-01
            1.31453e-01
            1.34094e-01
            2.26511e-01
            2.48770e-01
            2.52006e-01
            2.61382e-01
            3.66797e-01
            3.67082e-01
            3.99538e-01
            4.06751e-01
            4.11962e-01
            1.27359e+00
        ],
    ) < 0.001
    # savecsv("bend_hex_spectrum_im-aspect=$(aspect).csv", eigenvalues = vec(D.values))
    # @pgf _a = SemiLogYAxis({
    #     xlabel = "Mode [ND]",
    #     ylabel = "Generalized stiffness [N/m]",
    #     grid="major",
    #     legend_pos  = "south east",
    #     title = "Hexahedron spectrum, aspect=$(aspect)"
    # },
    # Plot({"red", mark="triangle"}, Table([:x => vec(7:size(K, 1)), :y => vec(D.values[7:end])])), LegendEntry("IM"))
    # display(_a)
    # File =  "bend_hex_spectrum_im.vtk"
    # vectors = [("ev_$(idx)_$(round(D.values[idx] * 10000) / 10000)", deepcopy(scattersysvec!(u, D.vectors[:,idx]).values)) for idx in 1:length(D.values)]
    # vtkexportmesh(File, fens, fes;  vectors=vectors)
    # @async run(`"paraview.exe" $File`)

    true
end

bend_hex_spectrum_im()
end

nothing
