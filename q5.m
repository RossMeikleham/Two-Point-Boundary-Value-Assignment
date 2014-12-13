% Code to calculate the solution to a 
% nonlinear BVP with Dirichlet and Robin BCs

% Some physical parameter.
mu = 10;

% Define the interval over which solution is calculated.
a= 0; b= 1;

% Define N.
N= 200;

% Define the grid spacing.
h=(b-a)/N;

% Define the grid. Note that there is no need to solve at x=a.
x = reshape(linspace(a+h,b,N),N,1);

% Define an initial guess for the solution 
% of the BVP (make sure it is a column vector)
U= zeros(N, 1);

% Define a tolerance for the termination of Newton-Raphson.
tol= 10^(-8); 
% Ensure that F is such that at least one iteration is done.
F=ones(N,1);
J=zeros(N,N);

% Store the initial guess in SOL.
SOL= U;

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
    %J(N,:) = horzcat(zeros(1, N-2),  2/(h^2), -2/(h^2) - (U(N)^2) * exp(U(N))*(U(N) + 3));
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

    %% Store the new approximation.
    SOL= [SOL,U];
    size(SOL)
end

%% Insert your plot commands here.
plots = []
labels = []
sz = size(SOL);
for i  = 1:(sz(2))    
    if (i > 4)
        p = plot(x, exp(SOL(:,i)));
        %text(0.5, SOL(i + (N/2) - 1), strcat('iteration  ', int2str(iterNo)));
        plots = [plots, p];
        %labels = [labels, strcat('iteration  ', int2str((i - 1)/N + 1))];
        hold all;
    end
end
%legend(plots, labels);
hold off

