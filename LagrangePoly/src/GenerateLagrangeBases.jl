module GenerateLagrangeBases

# import Pkg;
# Pkg.add("DynamicPolynomials")

using DynamicPolynomials, LinearAlgebra, Combinatorics

function multinomial_coefficient(powers::Vector{Int})
    return prod(factorial(p) for p in powers)
end

function generate_monomials_dynamic(n::Integer, m::Integer, x)
    monomial_powers = exponent_combinations(n, m)
    
    monomials = []
    for powers in reverse(monomial_powers)
        coefficient = multinomial_coefficient(powers)
        monomial = prod(x[i]^powers[i] for i in 1:n) / coefficient
        push!(monomials, monomial)
    end
    
    return monomials
end

function exponent_combinations(n::Integer, m::Integer)
    comb = [tuple for tuple in multinomial(n, m)]
    return comb
end

function multinomial(n::Integer, m::Integer)
    nomials = []
    for stars in combinations(1:m+n-1, n-1)
        result = [stars[1]-1; diff(stars) .- 1; m + n - 1 - stars[end]]
        push!(result, m - sum(result))
        push!(nomials, result)
    end

    return nomials
end

function compute_basis_matrix(data_points, basis, x)
    basis_matrix = []
    for i in axes(data_points, 1)
        inps = data_points[i, :]
        values = Dict(x[j] => inps[j] for j in eachindex(inps))
        entry = [subs(b, values...) for b in basis]
        push!(basis_matrix, entry)
    end

    basis_matrix = hcat(basis_matrix...)
    basis_matrix = Array{Float64}(basis_matrix)

    # convert(Matrix(Float64), basis_matrix)
    # Float64.(basis_matrix)
    return basis_matrix
end

function build_lagrange_bases_frobenius(data_points, basis, x)
    M = compute_basis_matrix(data_points, basis, x)

    n = size(data_points, 2)
    M_l = M[1:n+1, :]' 
    M_q = M[n+2:end, :]' 
    
    Vb_l = basis[1:n+1]
    Vb_q = basis[n+2:end]

    A = M_q * M_q'
    F = vcat(hcat(A, M_l), hcat(M_l', zeros(size(M_l, 2), size(M_l, 2))))

    Finv = pinv(real(F))
    B = vcat(M_q * Vb_q, Vb_l)
    S = Finv * B
    
    polys = [S[i] for i in 1:size(A, 1)]

    return polys
end

function generate_lagrange_bases(n::Integer, m::Integer, d::Matrix{Float64})
    @polyvar x[1:n]
    linear_monomials = generate_monomials_dynamic(n, m-1, x)
    quadratic_monomials = generate_monomials_dynamic(n, m, x)

    monomials = []
    append!(monomials, [1])
    append!(monomials, linear_monomials)
    append!(monomials, quadratic_monomials)

    lpoly = build_lagrange_bases_frobenius(d, monomials, x)
    return lpoly, x
end

end