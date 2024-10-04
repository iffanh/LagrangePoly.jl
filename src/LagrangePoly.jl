module LagrangePoly

using Combinatorics
using CSDP
using DynamicPolynomials
using JuMP
using LinearAlgebra
using SumOfSquares

include("GenerateLagrangeBases.jl")
include("GenerateLagrangePoly.jl")

# GenerateLagrangePoly.jl
export model_improvement
export compute_poisedness
export _generate_lagrange_poly
export generate_lagrange_poly

# GenerateLagrangeBases.jl
export generate_lagrange_bases
export update_lagrange_bases


#export LagrangePolynomial ##TODO

end
