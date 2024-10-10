function _generate_lagrange_poly(Y::Vector{Float64}, lpolys)
    poly = 0
    for (y, b) in zip(Y, lpolys)
        poly += y*b
    end
    return poly
end

function generate_lagrange_poly(n::Int, m::Int, d::Matrix{Float64}, Y::Vector{Float64}, x::Vector{Variable{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}})
    lpolys = generate_lagrange_bases(n, m, d, x)
    return _generate_lagrange_poly(Y, lpolys)
end

function generate_lagrange_poly(n::Int, m::Int, d::Matrix{Float64}, Y::Vector{Float64})
    lpolys = generate_lagrange_bases(n, m, d)
    return _generate_lagrange_poly(Y, lpolys)
end

function compute_poisedness(center::Vector{Float64}, radius::Float64, lpolys)
    solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

    p = Vector{Float64}();
    for lpoly in lpolys
            model = SOSModel(solver)

            JuMP.@variable(model, α)
            JuMP.@objective(model, Min, α)

            x = DynamicPolynomials.variables(lpoly)
            the_sum = sum((x[j] - center[j])^2 for j in eachindex(center))
            S = SumOfSquares.@set the_sum <= radius^2
            JuMP.@constraint(model, c2, lpoly <= α, domain = S, maxdegree = 2)

            JuMP.optimize!(model)
            
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

function model_improvement!(ℓs, d::Matrix{Float64}, Δ::Float64, Λth::Float64)

    center = d[1, :]

    for k in range(1,size(d)[1])
        
        Λs = compute_poisedness(center, Δ, ℓs)
        Λ = max(Λs...)
        
        if Λ <= Λth
            break
        end
        
        imax = argmax(Λs)

        dn = find_point_that_maximizes_lagrange_base(ℓs[imax], center, Δ)
        ℓs, d = update_lagrange_bases(ℓs, d, dn, imax)
    end

    # poisedness is invariant under linear translation
    shift = d[1,:] - center
    d = d .- shift'
    
    return ℓs, d
end