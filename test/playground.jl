
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
	@show S - tS
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
	@show S - tS
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
	@show tr(S) * I - tS
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
	@show (S + S')/2 * I - tS
	@test norm((S + S')/2 - tS) <= 1.0e-12
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
	@show (S - S')/2 * I - tS
	@test norm((S - S')/2 - tS) <= 1.0e-12
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
	@show tr((S - tr(S)/3*I) ), tr(tS)
	@test norm(tr((S - tr(S)/3*I) ) - tr(tS)) <= 1.0e-12
	return true
end
end
using .mmmtdtens4C6
mmmtdtens4C6.test()
