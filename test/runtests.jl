using PauliStrings
using Test

@testset "PauliString - Structure & Basic Arithmetic" begin


	p_string = PauliString(4,[1,3],[2,3],3.3)

	#
	#	check construction 
	#
	p_string.N == 4
	p_string.sites == [1,3]
	p_string.baseIdx == [2,3]
	@test coef(p_string) == 3.3

	#
	#	check set_coef!
	#
	set_coef!(p_string, 2.0)
	@test coef(p_string) == 2.0


	#
	#	check multiplication
	#
	p_start = deepcopy(p_string)
	p_string = 3.0*p_string
	@test coef(p_string) == 6.0

	p_string = p_string*2.0
	@test coef(p_string) == 12.0

	#
	#	check equal and parallel:
	#	parallel means same site same baseIdx different coef
	#
	p_string∥p_start == true
	@test iszero(p_string == p_start)

	# check addition
	p_start += p_start # 4
	p_start *= 3.0 
	@test isone(p_string == p_start)


	#
	#	check conj
	#
	p_string *= 1.0im
	p_string_conj = conj(p_string)

	@test coef(p_string_conj) == -12.0im

	conj!(p_string)
	@test isone(p_string_conj == p_string)

end



@testset "PauliStrings - Convertion to SparseArray" begin

	@show nothing


end



@testset "PauliStrings - Commutation & Multiplication" begin

	# test multiplication
	p1 = PauliString(1,[1],[2],1.0)
	p2 = PauliString(1,[1],[3],1.0)
	p3 = PauliString(1,[1],[1],1.0im)
	@test p3 ==  p1*p2
	@test p3 == conj(p2*p1)


	p1 = PauliString(2,[1,2],[2,1],1.0)
	p2 = PauliString(2,[1],[3],1.0)
	p3 = PauliString(2,[1,2],[1,1],1.0im)

	@test p3 ==  p1*p2
	@test p3 == conj(p2*p1)

	#commutation
	@test commutePauli(p1,p2) == 2*p3

	p1 = PauliString(2,[1,2],[2,1],1.0)
	p2 = PauliString(2,[1,2],[3,2],1.0)
	p3 = PauliString(2,[1,2],[1,3],-1.0)

	@test p3 ==  p1*p2
	@test p3 == p2*p1


	# commutation
	@test commutePauli(p1,p2) == nothing

end
