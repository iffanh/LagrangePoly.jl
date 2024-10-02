using LagrangePoly

include("GenerateLagrangeBases.jl");

n = 2; # number of variables
m = 2; # polynomial degree

# data_points = [1.0 2.0; 4.0 5.0; 5.0 6.0; 7.0 8.0; 9.0 10.0];  # Example data points
data_points = [0.5 0.5; 
               0.524 0.0006; 
               0.032 0.323;
               0.187 0.890; 
               0.982 0.368;
               0.774 0.918]

myPoly = GenerateLagrangeBases.generate_lagrange_bases(n, m, data_points);


new_point = [0.2, 0.2];
newPoly, new_data_points = GenerateLagrangeBases.update_lagrange_bases(myPoly, data_points, new_point, 3);



new_point = [0.6, 0.6];
newPoly, new_data_points = GenerateLagrangeBases.update_lagrange_bases(newPoly, new_data_points, new_point, 6);



# @polyvar x[1:n]
# println(x[1], x[2])
# linear_monomials = GenerateLagrangeBases.generate_monomials_dynamic(n, m-1, x)
# quadratic_monomials = GenerateLagrangeBases.generate_monomials_dynamic(n, m, x)


# data_points2 = [4.0 5.0; 5.0 6.0; 7.0 8.0; 9.0 10.0; 1.0 2.0]  # Example data points

# myPoly2 = GenerateLagrangeBases.generate_lagrange_bases(n, m, data_points2)


# data_points3 = [4.2 5.1; 5.3 6.2; 7.1 8.2; 9.3 10.2; 1.1 2.1]  # Example data points

# myPoly3 = GenerateLagrangeBases.generate_lagrange_bases(n, m, data_points3)


include("GenerateLagrangePoly.jl");

using SumOfSquares

Y = [1.0, 2.0, 3.0, 4.0, 5.0];
GenerateLagrangePoly.generate_lagrange_poly(Y, myPoly);


GenerateLagrangePoly.compute_poisedness(data_points[1,:], 0.5, myPoly);

p, d = GenerateLagrangePoly.model_improvement(newPoly, new_data_points, 0.5, 1.2);

