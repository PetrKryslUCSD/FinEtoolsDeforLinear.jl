module mmmtdt1
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
function test()
    MR = DeforModelRed3D
    symmet(a) = a + transpose(a)

    a = rand(3, 3)
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
    a = rand(3, 3)
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
    a = rand(3, 3)
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
    a = rand(3, 3)
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
    a = rand(3, 3)
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
    a = rand(3, 3)
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    @test abs(det(a) - dett(MR, a)) / abs(det(a)) <= 1.0e-6

    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(det(a) - strainvdet(MR, av)) / abs(det(a)) <= 1.0e-6

    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
    a = symmet(rand(3, 3))
    av = fill(zero(eltype(a)), 6)
    strainttov!(MR, av, a)
    @test abs(tr(a) - strainvtr(MR, av)) / abs(tr(a)) <= 1.0e-6
end
end
using .mmmtdt1
mmmtdt1.test()


module mmmtdtens4A1
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156
    F = [1.02 0.03 -0.04; 0.01 0.99 -0.03; -0.01 0.02 0.95]
    C = fill(0.0, 3, 3, 3, 3)
    for I = 1:3, J = 1:3, K = 1:3, L = 1:3
        C[I, J, K, L] =
            lambda * delta(I, J) * delta(K, L) +
            mu * (delta(I, K) * delta(J, L) + delta(I, L) * delta(J, K))
    end
    Cm = fill(0.0, 6, 6)
    tens4symmtto6x6t!(Cm, C)

    mI = Matrix(Diagonal([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]))
    m1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    D = lambda * m1 * m1' + 2.0 * mu * mI

    @test norm(Cm - D) <= 1.0e-6
    return true
end
end
using .mmmtdtens4A1
mmmtdtens4A1.test()


module mmmtdtens4B1
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156

    t = fill(0.0, 3, 3, 3, 3)
    delta = (I, J) -> I == J ? 1.0 : 0.0
    tens4ijkl!(t, delta, delta)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    @test norm(tr(S) * I - tens4dot2!(tS, t, S)) <= 1.0e-12

    # @btime tens4dot2!($tS, $t, $S)
    return true
end
end
using .mmmtdtens4B1
mmmtdtens4B1.test()

module mmmtdtens4B2
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156

    t = fill(0.0, 3, 3, 3, 3)
    delta = (I, J) -> I == J ? 1.0 : 0.0
    tens4iljk!(t, delta, delta)
    S = rand(3, 3)
    # @show transpose(S) 
    tS = fill(0.0, 3, 3)
    # @show transpose(S) - tens4dot2!(tS, t, S)
    @test norm(transpose(S) - tens4dot2!(tS, t, S)) <= 1.0e-12

    # @btime tens4dot2!($tS, $t, $S)
    return true
end
end
using .mmmtdtens4B2
mmmtdtens4B2.test()

module mmmtdtens4B3
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156

    t = fill(0.0, 3, 3, 3, 3)
    delta = (I, J) -> I == J ? 1.0 : 0.0
    tens4ikjl!(t, delta, delta)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    # @show S - tens4dot2!(tS, t, S)
    @test norm(S - tens4dot2!(tS, t, S)) <= 1.0e-12

    # @btime tens4dot2!($tS, $t, $S)
    return true
end
end
using .mmmtdtens4B3
mmmtdtens4B3.test()

module mmmtdtens4B4
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156

    t = fill(0.0, 3, 3, 3, 3)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4identity!(t)
    @test norm(S - tens4dot2!(tS, t, S)) <= 1.0e-12
    # @btime tens4dot2!($tS, $t, $S)
    return true
end
end
using .mmmtdtens4B4
mmmtdtens4B4.test()

module mmmtdtens4B5
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156

    t = fill(0.0, 3, 3, 3, 3)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4transposor!(t)
    @test norm(S' - tens4dot2!(tS, t, S)) <= 1.0e-12
    # @btime tens4dot2!($tS, $t, $S)
    return true
end
end
using .mmmtdtens4B5
mmmtdtens4B5.test()

module mmmtdtens4B6
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156

    t = fill(0.0, 3, 3, 3, 3)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4tracor!(t)
    @test norm(tr(S) * 1I - tens4dot2!(tS, t, S)) <= 1.0e-12
    # @btime tens4dot2!($tS, $t, $S)
    return true
end
end
using .mmmtdtens4B6
mmmtdtens4B6.test()

module mmmtdtens4B8
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    delta = (I, J) -> I == J ? 1.0 : 0.0
    lambda = 3.3
    mu = 0.156

    t = fill(0.0, 3, 3, 3, 3)
    S = rand(3, 3)
    tens4skewor!(t)
    tS = rand(3, 3)
    @test norm((S - S') / 2 - tens4dot2!(tS, t, S)) <= 1.0e-12
    # @btime tens4dot2!($tS, $t, $S)
    return true
end
end
using .mmmtdtens4B8
mmmtdtens4B8.test()

module mmmtdtens4C1
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    t = fill(0.0, 3, 3, 3, 3)
    tens4identity!(t)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4dot2!(tS, t, S)
    # @show S - tS
    @test norm(S - tS) <= 1.0e-12
    return true
end
end
using .mmmtdtens4C1
mmmtdtens4C1.test()

module mmmtdtens4C2
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    t = fill(0.0, 3, 3, 3, 3)
    tens4transposor!(t)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4dot2!(tS, t, S)
    # @show S - tS
    @test norm(S' - tS) <= 1.0e-12
    return true
end
end
using .mmmtdtens4C2
mmmtdtens4C2.test()

module mmmtdtens4C3
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    t = fill(0.0, 3, 3, 3, 3)
    tens4tracor!(t)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4dot2!(tS, t, S)
    # @show tr(S) * I - tS
    @test norm(tr(S) * I - tS) <= 1.0e-12
    return true
end
end
using .mmmtdtens4C3
mmmtdtens4C3.test()

module mmmtdtens4C4
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    t = fill(0.0, 3, 3, 3, 3)
    tens4symmetrizor!(t)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4dot2!(tS, t, S)
    # @show (S + S')/2 * I - tS
    @test norm((S + S') / 2 - tS) <= 1.0e-12
    return true
end
end
using .mmmtdtens4C4
mmmtdtens4C4.test()

module mmmtdtens4C5
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    t = fill(0.0, 3, 3, 3, 3)
    tens4skewor!(t)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4dot2!(tS, t, S)
    # @show (S - S')/2 * I - tS
    @test norm((S - S') / 2 - tS) <= 1.0e-12
    return true
end
end
using .mmmtdtens4C5
mmmtdtens4C5.test()

module mmmtdtens4C6
using FinEtools
using FinEtoolsDeforLinear
using LinearAlgebra
using Test
# using BenchmarkTools

function test()
    t = fill(0.0, 3, 3, 3, 3)
    tens4deviator!(t)
    S = rand(3, 3)
    tS = fill(0.0, 3, 3)
    tens4dot2!(tS, t, S)
    # @show tr((S - tr(S)/3*I) ), tr(tS)
    @test norm(tr((S - tr(S) / 3 * I)) - tr(tS)) <= 1.0e-12
    return true
end
end
using .mmmtdtens4C6
mmmtdtens4C6.test()
