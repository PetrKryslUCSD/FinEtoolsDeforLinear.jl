
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
	tens4ikjl!(t, delta, delta)
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
	tens4iljk!(t, delta, delta)
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
