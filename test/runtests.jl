using TestItemRunner

@run_package_tests verbose = true

@testitem "Generation of full quadratic polynomial" begin
    using DynamicPolynomials
    using LagrangePoly


    n = 2; # number of variables
    m = 2; # polynomial degree

    # From Chapter 3 in Conn's book
    data_points = [0.0 0.0; 
                1.0 0.0;
                0.0 1.0; 
                2.0 0.0;
                1.0 1.0; 
                0.0 2.0];


    myPoly = generate_lagrange_bases(n, m, data_points);

    @test typeof(myPoly) == Vector{DynamicPolynomials.Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, MultivariatePolynomials.Graded{MultivariatePolynomials.LexOrder}, Float64}} ## TODO: replace this with something meaningful

    # test the bases
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

    # test the lagrange identity
    for i in eachindex(myPoly)
        for j in eachindex(myPoly)
            if i == j
                val = 1
            else
                val = 0
            end
            @test abs(myPoly[i](data_points[j,:]) - val) < tol
        end
    end 

    # test the poisedness
    data_points = [0.5 0.5; 
                    0.524 0.0006; 
                    0.032 0.323;
                    0.187 0.890; 
                    0.982 0.368;
                    0.774 0.918]

    myPoly = generate_lagrange_bases(n, m, data_points);
    p = LagrangePoly.compute_poisedness(data_points[1,:], 0.5, myPoly)
    @test abs(max(p...) - 1) < 1e-3
end

@testitem "Generation of linear polynomial" begin
    
    using DynamicPolynomials
    using LagrangePoly


    n = 2; # number of variables
    m = 2; # polynomial degree

    data_points = [0.0 1.0; 
                1.0 0.0];


    myPoly = generate_lagrange_bases(n, m, data_points); ## should generate only linear terms

    tol = 1e-8
    @test abs(myPoly[1].a[1] - 1/3) < tol
    @test abs(myPoly[1].a[2] - 2/3) < tol
    @test abs(myPoly[1].a[3] - -1/3) < tol

    @test abs(myPoly[2].a[1] - 1/3) < tol
    @test abs(myPoly[2].a[2] - -1/3) < tol
    @test abs(myPoly[2].a[3] - 2/3) < tol

    for i in eachindex(myPoly)
        for j in eachindex(myPoly)
            if i == j
                val = 1
            else
                val = 0
            end
            @test abs(myPoly[i](data_points[j,:]) - val) < tol
        end
    end 

end

