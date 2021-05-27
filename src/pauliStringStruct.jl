#
# PauliString struct
#

struct PauliString
	N::Int  # length 
	sites::Vector{Int} # To do: change this to bit number for "used sites" 000101 = operator at site 4 and 6
	# vector which contains site indicator where the string has non-trivial operators 
	baseIdx::Vector{Int} # indicates the operator (σx, σy, σz)
	coef::Base.RefValue{ComplexF64}


	function PauliString(N::Int, site::Vector{Int}, baseIdx::Vector{Int}, coef::Number)
		length(baseIdx[baseIdx.>3]) == 0 || throw(ArgumentError("baseIdx $(baseIdx) containes unvalid arguement >3.
									Note the convention is: (1,2,3) = (σx, σy, σz)"))
		return new(N, site, baseIdx, Ref(convert(ComplexF64, coef)))
	end
end


coef(p::PauliString) = p.coef[]
set_coef!(p::PauliString, a::Number) = p.coef[] = a

import Base.==
import Base.+
import Base.*
import Base.conj
(==)(a::PauliString, b::PauliString) = a.sites == b.sites && a.baseIdx == b.baseIdx && coef(a) == coef(b) ? true : false
(∥)(a::PauliString, b::PauliString) = a.sites == b.sites && a.baseIdx == b.baseIdx ? true : false
(+)(a::PauliString, b::PauliString) = a∥b ? PauliString(a.N, a.sites, a.baseIdx, coef(a)+coef(b)) : [a,b]
(*)(a::PauliString, c::Number) = PauliString(a.N, a.sites, a.baseIdx, coef(a)*c)
(*)(c::Number, a::PauliString) = PauliString(a.N, a.sites, a.baseIdx, coef(a)*c)

function (*)(a::PauliString, b::PauliString)

	p1_sites = a.sites
	p2_sites = b.sites
	p1_idx = a.baseIdx.+1 #  in order that 0 (id) are not get cut out later on 
	p2_idx = b.baseIdx.+1

	pArray = sparse([[1 for i in p1_sites]...,
			[2 for i in p2_sites]...],
			[p1_sites...,p2_sites...],[p1_idx..., p2_idx...])
	new_σ = []
	sites = [Set([x[2] for x in  findall(x->x!=0,pArray)])...][end:-1:1]
	#@show sites	
	for col in sites
		#@show pArray[1,col], pArray[2,col]
		i = pArray[1,col] == 0 ? 1 : pArray[1,col]
		j = pArray[2,col] == 0 ? 1 : pArray[2,col]
		push!(new_σ, σσ(i, j))  # +1 to match the convention 
     					   # 0 = σ⁰, 1 = σˣ , ... 
	end
	
	new_idx = [Int(new_σ[i][1])-1 for i in 1:size(sites)[1]]
	c = *([new_σ[i][2] for i in 1:size(sites)[1]]...)
	
	idx = sortperm(sites)
	sites  = sites[idx]
	new_idx = new_idx[idx]

	return PauliString(a.N, sites, new_idx, c*coef(a)*coef(b))
end

function Base.conj!(a::PauliString)	
 	set_coef!(a, (coef(a))')
	return a
end

function Base.conj(a::PauliString)	
	b = deepcopy(a)
	set_coef!(b, (coef(b))')
	return b
end

function Base.show(io::IO, p::PauliString)
	print_string = fill("𝟙-",p.N)
	op_symbols = ["σˣ-","σʸ-","σᶻ-"]
	print_string[p.sites] .= op_symbols[p.baseIdx]

	output_string = chop(string("($(coef(p)))⋅",print_string...))
	print(io, output_string)
end

function rev_sites(a::PauliString)
	return PauliString(a.N, [a.N.-(a.sites.-1)...], a.baseIdx, coef(a))
end

# get cycles of a permutation 
function getcycles(p::Vector{T}) where {T<:Any}
	isperm(p) || throw(DomainError("p is not a permutation"))

	n = length(p)
	used = falses(n)
	cycles = []
	
	for k = 1:n
		if used[k]; continue; end
		used[k] = true
		j = p[k]
			
		cycle = [j]
		# go through one cycle
		while !used[j]
			used[j] = true
			j = p[j]
			push!(cycle, j) 
		end

		if length(cycle)>1 
			cycle[1:end] = [cycle[end], cycle[1:end-1]...]
		end
		push!(cycles, cycle)

	end
	return cycles
end


# see 
# https://github.com/toivoh/ConicHulls.jl/blob/888109ab6fe7d334ffa48e649947f0c293a967ae/src/Common.jl#L37
function isevenperm(p::Vector{T}) where {T<:Any}
	isperm(p) || throw(DomainError("p is not a permutation"))

	n = length(p)
	used = falses(n)
	even = true

	for k = 1:n
		if used[k]; continue; end

		# Each even cycle flips even (an odd number of times)
		used[k] = true
		j = p[k]
	
		# go through one cycle
		while !used[j]
			used[j] = true
			j = p[j]
			even = !even
		end
	end
	
	return even 
end
isoddperm(p) = !isevenperm(p)

"""
    σσ(i,j)

calculates σᵢσⱼ = δᵢⱼ I   + im*ϵᵢⱼₗσₗ


# Arguments:
    - i:
    - j:

return:
    
"""
function σσ(i::Int, j::Int)
	coef = 1.0
	#@show i, j
	# input is σ⁰ = 1, σˣ = 2, ...
	i-=1; j-=1 # change to σ⁰ = 0, σˣ = 1, ...

	if i == 0 || j == 0 #either i or j is σ⁰
		return i+j+1, coef
	elseif i == j 
		return 1, coef
	else

		# determin σ	
		σ_res = filter(x->x ∉ [i,j],[1,2,3])[1]
		
		# determin coeff
		coef *= isoddperm([i,j,σ_res]) ? -1.0im : 1.0im 
	
		return σ_res+1, coef
	end
end


"""
	[σᵢ , σⱼ] = 2*im*ϵᵢⱼₖ σₖ


"""
function com_σσ(i::Int, j::Int)
	i-=1; j-=1
	coef = 1
	if i == 0 || j == 0
		return 0,0
	elseif i == j
		return 0,0
	else
		# determin σ	
		σ_res = filter(x->x ∉ [i,j],[1,2,3])[1]
	
		# determin coeff
		coef *= isoddperm([i,j,σ_res]) ? -1.0im : 1.0im 
	
		return σ_res+1, 2*coef
	end
end

function commutePauli(p1::PauliString, p2::PauliString)	
	p1_sites = p1.sites
	p2_sites = p2.sites
	p1_idx = p1.baseIdx
	p2_idx = p2.baseIdx

	pArray = sparse([[1 for i in p1_sites]...,
			[2 for i in p2_sites]...],
			[p1_sites...,p2_sites...],[p1_idx..., p2_idx...])

	# get unique sites to avoide double counting
	sites = [Set([x[2] for x in  findall(x->x!=0,pArray)])...][end:-1:1]
	new_σσ = []	
	for x in sites

		# +1 for convention: label 1 is used for σ0 in σσ routine
		o1 = pArray[1,x]+1
		o2 = pArray[2,x]+1
		
		o1_o2 = σσ(o1, o2)
		o2_o1 = σσ(o2, o1)
		push!(new_σσ,(o1_o2, o2_o1))
	end


	coef_op = *([x[1][2] for x in new_σσ]...)-*([x[2][2] for x in new_σσ]...)
	
	if !iszero(coef_op)
		new_p = []
		new_site = []
		for (i,x) in enumerate(new_σσ)
			new_o = x[1][1]-1
			if !iszero(new_o)
				push!(new_p, new_o)
				push!(new_site, [sites...][i])
			end
		end
		res_p = PauliString(p1.N, Vector{Int}(new_site), Vector{Int}(new_p), coef_op*coef(p1)*coef(p2))
		
		return res_p
	
	end
end




"""
    adH(H, p, order)
calculates all elements of adH up to order n

# Arguments:


return: 
!!! To do put output in a dict specifying each order.
"""
function adH(H, p, order)
	pauliList = Vector{PauliString}([deepcopy(p)])
	new_pauliList = []
	for n in 1:order
		for pauli_string in pauliList			
			
			# get operator strings in H with support overlaping with pauli_string
			sup_string = getSupport(pauli_string)
			H_accesible = filter(x-> getSupport(x)[1]<=sup_string[2] && getSupport(x)[2]>=sup_string[1], H)
			res = filter(x->x!=nothing, [commutePauli(x,pauli_string)  for x in H_accesible])
			
			#res = filter(x->x!=nothing, [commutePauli(x,y)  for x in H])
			
			push!(new_pauliList, res...)

		end
		pauliList = Vector{PauliString}(new_pauliList)
	end

	if new_pauliList != []
		return pauliList
	else
		return nothing
	end
end



"""
    simplify_stringList(L; tol = nothing)



# Arguments:


return: 
"""
function simplify_stringList(L::Vector{PauliString}; tol = nothing)
	
	# remove operator strings with coef = 0 
	filter!(x->coef(x) != 0, L)
	# sort by number of sites wiht non-trivial operator
	sort!(L, by = x -> size(x.sites)[1]) 

	number_op = 1
	new_L = Vector{PauliString}([])
	while L != []

		# retrieve all operator with number_op non-trivial operator
		# and store them temporary in tmp
		idx = findlast(x->size(x.sites)[1] == number_op, L)
		tmp = idx != nothing ? [popfirst!(L) for i in 1:idx] : []

		# while tmp is not empty check for equivalent operator strings
		while tmp != []

			idx = findall(x->x∥tmp[1], tmp) # find all similar PauliStrings
			# if they sum to coef != 0  add them to the new_L
			if sum([coef.(tmp[idx])...]) != 0  

				# erase zeros from string
				non_zero_idx = findall(x->x!=0,tmp[idx][1].baseIdx)
				if non_zero_idx != []
					new_base_idx = tmp[idx][1].baseIdx[non_zero_idx]
					new_sites = tmp[idx][1].sites[non_zero_idx]
			
					sort_sites = sortperm(new_sites)
					new_sites = new_sites[sort_sites]
					new_base_idx = new_base_idx[sort_sites]
				
				else
				
					new_base_idx = [0]
					new_sites = [1]
				end

				push!(new_L, PauliString(tmp[idx][1].N, new_sites, new_base_idx, sum([coef.(tmp[idx])...])))
			end
			deleteat!(tmp, idx) # remove visited strings from tmp

		end
		
		number_op += 1 # go to next number_op sector

	end
	
	# filter strings which coef are below tol if tol!=nothing
	if tol !=nothing new_L = filter(x->abs(coef(x))>=tol,new_L) end

	return new_L
end

getSupport(p::PauliString)  = (minimum(p.sites),maximum(p.sites))

function estimateLoc(pvec::Vector{PauliString})
	norm = 1/sum(abs.(coef.(pvec)))

	loc = norm*sum(abs.(coef.(pvec)).*[x[2]-x[1] for x in getSupport.(pvec)])

	return loc
end
