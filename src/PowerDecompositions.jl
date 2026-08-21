module PowerDecompositions

using AlgebraicSolvers
using Combinatorics
using DynamicPolynomials
using LinearAlgebra
using MultivariatePolynomials
using Random
using TensorDec

export decompose_powers
export print_decomposition

include("preimages.jl")
include("catalecticants.jl")
include("rank_search.jl")
include("decomposition.jl")

end