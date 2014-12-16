function q7

n = 13;
start = 2;
res= zeros(1, n - start - 1);
vk = zeros(1, n);

% Compute all V(2^k) values
for k=2:n 
    vk(k) = v(2^k); 
end

% Compute ratios
for k = (start+2):n
    res(k-start-1) = abs(vk(k-1) - vk(k-2)) / abs(vk(k) - vk(k-1));
end

% Plot points
plot((start+2:n), res);

% Graph details
title('Graph of convergence of the finite-difference discretisation method')
xlabel('k'); % x-axis label
ylabel('|v(2^(k-1)) - v(2^(k-2))|/ |v(2^k) - v(2^(k-1))|'); % y-axis label
text(2.2 + start, res(1), strcat('k = 4, ratio =', num2str(res(1))));
text(4 + start, res(4) + 0.02, strcat('k = 7, ratio =', num2str(res(4))));
text(8.7 + start, res(10) - 0.01, strcat('k = 13, ratio =', num2str(res(10))));

end

function res = v(N)
% Code to calculate the solution to a 
% nonlinear BVP with Dirichlet and Robin BCs

% Define the interval over which solution is calculated.
a= 0; b= 1;

mu = 10;

% Define the grid spacing.
h=(b-a)/N;

% Define the grid. Note that there is no need to solve at x=a.
x = reshape(linspace(a+h,b,N),N,1);

% Define an initial guess for the solution 
% of the BVP (make sure it is a column vector)
U= ones(N, 1);

% Ensure that F is such that at least one iteration is done.
F=ones(N,1);
J=sparse(N,N);

% Loop 10 times
for j = 1:10

    %% Define F(u).
    % Boundary conditions.
    F(N) = (U(N) * ((2 * h * (U(N)^2)) + 3)) - (4 * (U(N-1))) + U(N-2);

    % Finite difference approximation to ODE at interior nodes.
    F(1) = (U(2) * (2 + h * exp(U(1)))) - (4 * U(1)) + (2 - (h * exp(U(1))));
    F(1) = F(1) - (2 * (h^2) * mu * sin(2 * pi * x(1)));

    for i=2:N-1
        F(i)= (U(i+1) * (2 + h * exp(U(i)))) - (4 * U(i)) + U(i-1)*(2 - h * exp(U(i)));
        F(i) = F(i) - (2 * (h^2) * mu * sin(2 * pi * x(i)));
    end;

    %% Define the Jacobian J.
    % First row corresponds to BC at x=0.
    J(1,:) = horzcat(((U(2) - 1) * h * exp(U(1))) - 4, 2 + (h * exp(U(1))), zeros(1, N-2));

    % Last row corresponds to BC at x=1.
    J(N,:) = horzcat(zeros(1, N-3), 1, -4, (4 * h * U(N)^2) + 3); 

    % Intermediate rows correspond to F(i)=...
    for i=2:(N-1)
        % Diagonal entries.
        J(i,i)= ((U(i+1) - U(i-1)) * h * exp(U(i))) - 4;
        % Left diagonal entries
        J(i, i-1) = 2 - (h * exp(U(i)));
        % Right diagonal entries
        J(i, i+1) = 2 + (h * exp(U(i)));
    end
    
    %% Having defined F and J, update the approximate solution to
    % the difference equations.  
    U=U-sparse(J)\F;

end

% Return U at X =1/2
res = U(N/2);

end
