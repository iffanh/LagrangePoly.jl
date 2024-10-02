using TestItemRunner

@run_package_tests verbose = true

@testitem "Simple case" begin
    using DynamicPolynomials
    using LagrangePoly


    n = 2; # number of variables
    m = 2; # polynomial degree

    # data_points = [1.0 2.0; 4.0 5.0; 5.0 6.0; 7.0 8.0; 9.0 10.0];  # Example data points
    data_points = [0.0 0.0; 
                1.0 0.0;
                0.0 1.0; 
                2.0 0.0;
                1.0 1.0; 
                0.0 2.0];


    myPoly = generate_lagrange_bases(n, m, data_points);

    @test typeof(myPoly) == Vector{DynamicPolynomials.Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, MultivariatePolynomials.Graded{MultivariatePolynomials.LexOrder}, Float64}} ## TODO: replace this with something meaningful
end


