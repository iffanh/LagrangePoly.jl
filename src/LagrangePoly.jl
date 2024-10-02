module LagrangePoly

using Combinatorics
using CSDP
using DynamicPolynomials
using JuMP
using LinearAlgebra
using SumOfSquares

include("GenerateLagrangeBases.jl")
include("GenerateLagrangePoly.jl")

export generate_lagrange_bases
#export LagrangePolynomial ##TODO

end
