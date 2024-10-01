module GenerateLagrangePoly

# import Pkg
# Pkg.add("Ipopt")
# Pkg.add("JuMP")

include("GenerateLagrangeBases.jl")
using DynamicPolynomials
using JuMP
using CSDP
using SumOfSquares

function generate_lagrange_poly(Y::Vector{Float64}, LBases)
    poly = 0
    for (y, b) in zip(Y, LBases)
        poly += y*b
    end
    return poly
end

function compute_poisedness(center::Vector{Float64}, radius::Float64, LBases)
    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

    p = Vector{Float64}();
    for LBase in LBases
            model = SOSModel(solver)

            @variable(model, α)
            @objective(model, Min, α)

            x = DynamicPolynomials.variables(LBase)
            the_sum = sum((x[j] - center[j])^2 for j in eachindex(center))
            S = SumOfSquares.@set the_sum <= radius^2
            @constraint(model, c2, LBase <= α, domain = S, maxdegree = 2)

            optimize!(model)
            
            # FOR DEBUGGING
            #println(solution_summary(model))
            #θ = moment_matrix(model[:c2]);
            #η = atomic_measure(θ, 1e-3);    
            #x_opt = η.atoms[1].center

            SumOfSquares.unregister(model, :c2)
            append!(p, solution_summary(model).objective_value)
    end

    return p
end


end