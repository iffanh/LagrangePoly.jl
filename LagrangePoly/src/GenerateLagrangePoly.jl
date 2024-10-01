module GenerateLagrangePoly

# import Pkg
# Pkg.add("Ipopt")
# Pkg.add("JuMP")

include("GenerateLagrangeBases.jl")
using DynamicPolynomials
using Ipopt
using JuMP

function generate_lagrange_poly(Y::Vector{Float64}, LBases)
    poly = 0
    for (y, b) in zip(Y, LBases)
        poly += y*b
    end
    return poly
end

function compute_poisedness(center::Vector{Float64}, radius::Float64, LBases)
    p = Vector{Float64}()
    for LBase in LBases
        N = length(center)
        model = Model(Ipopt.Optimizer)
        @variable(model, x[1:N]) 
        
        the_sum = sum((x[j] - center[j])^2 for j in 1:N)
        @NLconstraint(model, c1, sqrt(the_sum) <= radius)
        
         # Define a function to evaluate the polynomial
        function evaluate_polynomial(LBase, x_vals::Vector{Float64}...)
            values = Dict(Symbol("x$i") => x_vals[i] for i in eachindex(x_vals))
            return subs(LBase, values)
        end

        # Register the function with JuMP
        # JuMP.register(model, :evaluate_polynomial, N, (x...) -> evaluate_polynomial(LBase, x...); autodiff = true)
        JuMP.register(model, :LBase, N, (x...) -> evaluate_polynomial(LBase, x...); autodiff = true)

        # Define the objective function using @NLobjective
        @NLobjective(model, Max, evaluate_polynomial(LBase, x...))

        optimize!(model)

        # Extract the result if needed
        result = value.(x)
        println(result)
    end

    return
end


end