import Base.convert

function convert(::Type{Array{ComplexF64,2}}, p::PauliString)
   

end



#
#	many-body pauli-matrices for operator decomposition 
#
function get_Sz(N::Int, site::Int)
	# Important convention we start chain with site : 1 ..... and end with site N
	site = site - 1 # so site 1 corresponds to a bitshift of 0  

	
	vec_nidx = [0:1:(2^N-1)...] # vector with all nidx labels 0,1,2,3 for 0..0, 10..0, 010..0, 110..0,...

	vec = vec_nidx.>>site  # shift bit by site to the left so if there was a 1 at site: 
			 	# site in orginal it would be at leading bit now.
	
	vec = Vector{Int}(-1 .*(2 .*(vec.%2).-1)) # check if there is a 1 at bit at site: site it is now at the leading bit.
	 				    # we can check if vec entries are now odd or not .. if odd the bit at site :site was 
					    # previously 1 otherwise it was 0 ---> for 1 assigne +1 for 0 assigne -1.


	#vec = Vector{Int}(1 .*(2 .*(vec.%2).-1)) # check if there is a 1 at bit at site: site it is now at the leading bit
	 				    # we can check if vec entries are now odd or not .. if odd the bit at site :site was 
					    # previously 1 otherwise it was 0 ---> for 1 assigne +1 for 0 assigne -1



	return sparse(vec_nidx.+1, vec_nidx.+1, vec) # .+1 since julia starts array count with 1
end

function get_Sy(N::Int, site::Int)
	# Important convention we start chain with site : 1 ..... and end with site N
	site = site - 1 # so site 1 corresponds to a bitshift of 0  
	
	vec_nidx = [0:1:(2^N-1)...] # vector with all nidx labels 0,1,2,3 for 0..0, 10..0, 010..0, 110..0,...

	# well Sy ≈ S+ - S-
	# so we can incorporate this with an xor operation:
	
	vec_applied = vec_nidx .⊻ Int(2^site)
	# if there is a bit 1 at site in vec_nidx it will remove that bit 1 but if there isnt it will creat one! 
	
	# according to if their was a spin or not we have -1 or +1 
	vec_data = vec_nidx.>>site # see comment at the beginning
	vec_data = 2 .*(vec_data.%2).-1
	vec_data = 1im.*vec_data
	#vec_data = -1im.*vec_data

	 return sparse(vec_nidx.+1, vec_applied.+1, vec_data)
end

function get_Sx(N::Int, site::Int)
	# Important convention we start chain with site : 1 ..... and end with site N
	site = site - 1 # so site 1 corresponds to a bitshift of 0  
	
	vec_nidx = [0:1:(2^N-1)...] # vector with all nidx labels 0,1,2,3 for 0..0, 10..0, 010..0, 110..0,...
	
	# well Sx ≈ S+ + S-
	# so we can incorporate this with an xor operation:
	vec_applied = vec_nidx .⊻ Int(2^site)
	
	vec_data = ones(2^N)

	return sparse(vec_nidx.+1, vec_applied.+1, vec_data)
end

#
#	convert to PauliString
#
