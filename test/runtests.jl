using TestItemRunner

@run_package_tests verbose = true

@testitem "Simple case" begin
    using DynamicPolynomials
    using LagrangePoly


    n = 2; # number of variables
    m = 2; # polynomial degree

    data_points = [0.0 0.0; 
                1.0 0.0;
                0.0 1.0; 
                2.0 0.0;
                1.0 1.0; 
                0.0 2.0];


    myPoly = generate_lagrange_bases(n, m, data_points);

    @test typeof(myPoly) == Vector{DynamicPolynomials.Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, MultivariatePolynomials.Graded{MultivariatePolynomials.LexOrder}, Float64}} ## TODO: replace this with something meaningful

    # From Chapter 3 in Conn's book
    tol = 1e-8
    @test abs(myPoly[1].a[1] - 1) < tol 
    @test abs(myPoly[1].a[3] - -1.5) < tol 
    @test abs(myPoly[1].a[2] - -1.5) < tol 
    @test abs(myPoly[1].a[4] - 0.5) < tol 
    @test abs(myPoly[1].a[6] - 0.5) < tol 
    @test abs(myPoly[1].a[5] - 1) < tol 

    @test abs(myPoly[2].a[1] - 0) < tol
    @test abs(myPoly[2].a[3] - 2) < tol
    @test abs(myPoly[2].a[2] - 0) < tol
    @test abs(myPoly[2].a[4] - 0) < tol
    @test abs(myPoly[2].a[6] - -1) < tol
    @test abs(myPoly[2].a[5] - -1) < tol

    @test abs(myPoly[3].a[1] - 0) < tol
    @test abs(myPoly[3].a[3] - 0) < tol
    @test abs(myPoly[3].a[2] - 2) < tol
    @test abs(myPoly[3].a[4] - -1) < tol
    @test abs(myPoly[3].a[6] - 0) < tol
    @test abs(myPoly[3].a[5] - -1) < tol

    @test abs(myPoly[4].a[1] - 0) < tol
    @test abs(myPoly[4].a[3] - -0.5) < tol
    @test abs(myPoly[4].a[2] - 0) < tol
    @test abs(myPoly[4].a[4] - 0) < tol
    @test abs(myPoly[4].a[6] - 0.5) < tol
    @test abs(myPoly[4].a[5] - 0) < tol

    @test abs(myPoly[5].a[1] - 0) < tol
    @test abs(myPoly[5].a[3] - 0) < tol
    @test abs(myPoly[5].a[2] - 0) < tol
    @test abs(myPoly[5].a[4] - 0) < tol
    @test abs(myPoly[5].a[6] - 0) < tol
    @test abs(myPoly[5].a[5] - 1) < tol

    @test abs(myPoly[6].a[1] - 0) < tol
    @test abs(myPoly[6].a[3] - 0) < tol
    @test abs(myPoly[6].a[2] - -0.5) < tol
    @test abs(myPoly[6].a[4] - 0.5) < tol
    @test abs(myPoly[6].a[6] - 0) < tol
    @test abs(myPoly[6].a[5] - 0) < tol

end


