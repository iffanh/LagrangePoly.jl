function generate_lagrange_poly(Y::Vector{Float64}, lpolys)
    poly = 0
    for (y, b) in zip(Y, lpolys)
        poly += y*b
    end
    return poly
end

function compute_poisedness(center::Vector{Float64}, radius::Float64, lpolys)
    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

    p = Vector{Float64}();
    for lpoly in lpolys
            model = SOSModel(solver)

            @variable(model, α)
            @objective(model, Min, α)

            x = DynamicPolynomials.variables(lpoly)
            the_sum = sum((x[j] - center[j])^2 for j in eachindex(center))
            S = SumOfSquares.@set the_sum <= radius^2
            @constraint(model, c2, lpoly <= α, domain = S, maxdegree = 2)

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


function find_point_that_maximizes_lagrange_base(lpoly, center::Vector{Float64}, radius::Float64)
    
    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
    model = SOSModel(solver)

    JuMP.@variable(model, α)
    JuMP.@objective(model, Min, α)

    x = DynamicPolynomials.variables(lpoly)
    the_sum = sum((x[j] - center[j])^2 for j in eachindex(center))
    S = SumOfSquares.@set the_sum <= radius^2
    @constraint(model, c2, lpoly <= α, domain = S, maxdegree = 2)
    optimize!(model)

    θ = moment_matrix(model[:c2]);
    η = atomic_measure(θ, 1e-3);    
    x_opt = η.atoms[1].center

    return x_opt
end

function model_improvement(lpolys, d::Matrix{Float64}, radius::Float64, Λth::Float64)

    center = d[1, :]
    for k in size(d, 1)
        
        Λs = compute_poisedness(center, radius, lpolys)
        Λ = max(Λs...)
            
        if Λ <= Λth
            break
        end

        imax = argmax(Λs)

        dn = find_point_that_maximizes_lagrange_base(lpolys[imax], center, radius)
        lpolys, d = GenerateLagrangeBases.update_lagrange_bases(lpolys, d, dn, imax)
    end

    return lpolys, d
end