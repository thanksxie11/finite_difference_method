n = 400; % For 100 samples
% For example, if the bounds for the variables are:
% var1: [10, 20], var2: [0, 5], var3: [-5, 5], var4: [2, 8]
bounds = [60, 60;
          2000000, 2000000;
          4, 6;
          8.65, 10.5];
result = lhs_four_variables(n, bounds);
disp(result);
function samples = lhs_four_variables(n, bounds)
    % n: number of samples
    % bounds: 4x2 matrix, each row is the [lower_bound, upper_bound] for each variable
    
    % Generate Latin Hypercube Sample of four variables
    lhsample_matrix = lhsdesign(n,4);
    
    % Scale the LHS matrix according to bounds
    for i = 1:4
        lhsample_matrix(:,i) = bounds(i,1) + (bounds(i,2)-bounds(i,1))*lhsample_matrix(:,i);
    end
    
    samples = lhsample_matrix;
end
