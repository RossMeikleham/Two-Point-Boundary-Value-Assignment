% ID: 1107023
function res = bvp(mu)
% Code to calculate the solution to a 
% nonlinear BVP with Dirichlet and Robin BCs

% Define the interval over which solution is calculated.
a= 0; b= 1;

% Define N.
N = 100;

% Define the grid spacing.
h=(b-a)/N;

% Define the grid. Note that there is no need to solve at x=a.
x = reshape(linspace(a+h,b,N),N,1);

% Define an initial guess for the solution 
% of the BVP (make sure it is a column vector)
U= ones(N, 1);

% Define a tolerance for the termination of Newton-Raphson.
tol= 10^(-8);

% Ensure that F is such that at least one iteration is done.
F=ones(N,1);
J=zeros(N,N);

% Loop while the norm(F) is greater than tol.
while (norm(F)>tol )

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
    U=U-J\F;

end

% Return values of U at X = 1/4, X = 1/2, and X = 3/4
res = [U(N/4), U(N/2), U((3*N)/4)];

end
