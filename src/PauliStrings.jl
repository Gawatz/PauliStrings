module PauliStrings

using Combinatorics
using SparseArrays
using LinearAlgebra





export PauliString
export coef, set_coef!, conj, conj!, commutePauli
export adH, simplify_stringList, ∥, show




include("pauliStringStruct.jl")
include("convertPauliString.jl")

# to do: some routines for operator evolution
#include("operatorEvo.jl")






end # module
