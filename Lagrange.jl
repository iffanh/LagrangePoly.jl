import Pkg; 
NDim = 3
NOrder = 2


using Combinatorics
using Symbolics
using LinearAlgebra

function multinomial_coefficient(powers::Vector{Int})
    # return factorial(sum(powers)) / prod(factorial(p) for p in powers)
    return prod(factorial(p) for p in powers)
end

function generate_monomials_symbolic(n::Integer, m::Integer, x::Vector{Num})
    # Generate all combinations of exponents with sum equal to m and n variables
    monomial_powers = exponent_combinations(n, m)
    
    # Generate monomials with normalization coefficients
    monomials = []
    for powers in reverse(monomial_powers)
        coefficient = multinomial_coefficient(powers)
        monomial = prod(x[i]^powers[i] for i in 1:n)/coefficient
        push!(monomials, monomial)
    end
    
    return monomials
end

# Helper function to generate exponent combinations that sum to m
function exponent_combinations(n::Integer, m::Integer)
    # Generate all combinations of n non-negative integers that sum up to exactly m
    comb = [tuple for tuple in multinomial(n, m)]
    return comb
end

# Helper function to compute multinomial (combinations of non-negative integers summing to m)
function multinomial(n::Integer, m::Integer)
    # We use the stars and bars approach to generate the partitions
    # We choose n-1 dividers (bars) and m stars, and iterate over all possibilities
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
        entry = substitute.(basis, Ref(values))
        push!(basis_matrix, entry)
    end

    basis_matrix = hcat(basis_matrix...)

    return basis_matrix
end

function build_lagrange_polynomials_frobenius(data_points, basis, x, is_gen=false)
    #=
    Responsible for generating the Lagrange polynomials in Frobenius norm sense.

    Args:
        basis: List of polynomial base
        data_points: matrix of interpolation set input Y {y1, y2, ..., yp}

    Returns:
        List of Lagrange polynomials, {lambda_1, lambda_1, ..., lambda_p}
    =#

    M = compute_basis_matrix(data_points, basis, x)
    println(typeof(M))

    # Divide between linear and quadratic terms
    n = size(data_points, 2)
    M_l = M[1:n+1, :]' # size = (p, n+1)
    M_q = M[n+2:end, :]' # size = (p, b - (n+1))
    
    Vb_l = basis[1:n+1]
    Vb_q = basis[n+2:end]

    # Build matrix F (eq 5.7)
    A = M_q * M_q'

    println(typeof(A), typeof(M_q), typeof(M_l))
    
    F = vcat(hcat(A, M_l), hcat(M_l', zeros(size(M_l, 2), size(M_l, 2))))

    Finv = pinv(real(F))
    
    # Construct (Eq 5.9) F x L = B
    B = vcat(M_q * Vb_q, Vb_l)
    S = Finv * B
    
    polys = [simplify(S[i], expand=true) for i in 1:size(A, 1)]

    return polys
end

# Example usage

n = 2 # Dimension (number of variables)
m = 2  # Order (total degree of polynomial)
x = Symbolics.variables(:x, 1:n)
linear_monomials = generate_monomials_symbolic(n, m-1, x)
quadratic_monomials = generate_monomials_symbolic(n, m, x)

monomials = []
append!(monomials, [1])
append!(monomials, linear_monomials)
append!(monomials, quadratic_monomials)

# Print all possible monomials
println(monomials)

data_points = [1.0 2.0; 4.0 5.0; 5.0 6.0; 7.0 8.0; 9.0 10.0]  # Example data points

lpolynomials = build_lagrange_polynomials_frobenius(data_points, monomials, x)
