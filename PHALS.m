function [U, obj_vec,opt_vec, time_vec]=PHALS(A,U0,maxIter, maxtime)
% PHALS: The PHALS method for solving the symmetric non-negative matrix factorization (SymNMF) problem.
%
% Inputs:
%   A        - Input data matrix (usually symmetric and non-negative).
%   U0       - Initial guess for the factor matrix U.
%   maxIter  - Maximum number of iterations.
%   maxtime  - Maximum allowable computation time (in seconds).
%
% Outputs:
%   U        - The final factor matrix obtained.
%   obj_vec  - Vector of objective values at each iteration.
%   opt_vec  - Vector of optimality gap values at each iteration.
%   time_vec - Vector of cumulative time taken at each iteration.
%

% Initialization
U=U0; % Initialize U with the input matrix U0
[~,r]=size(U); % Get the number of columns (rank) of U
tempA=sum(A(:).^2);% Compute the squared Frobenius norm of A
%tempA=norm(A,'fro')^2;
obj_vec=tempA-2*trace(U'*A*U)+norm(U'*U,'fro')^2; % Initial objective value
opt_vec=norm(U-max(U-4*(U*(U'*U) - A*U), 0),'inf'); % Initial optimality gap
time_vec=0; % Start tracking time

% Main iteration loop
for i=1:maxIter
    t1=tic; % Start timing the current iteration
    W=A*U; % Compute A*U
    Q=U'*U; % Compute U^T * U
    df=0; % Initialize the objective change for the iteration

    % Loop over each column of U
    for j=1:r
        % Compute the update direction for the j-th column
        d=max((W(:,j)-U*Q(:,j))/Q(j,j),-U(:,j));
        dd=d'*d; % Compute norm of d
        dU=d'*U; % Compute d^T * U

        % Define coefficients for the quartic polynomial for step size calculation
        p1=4*dd^2;
        p2=12*dU(j)*dd;
        p3=4*(dU(j)^2+dd*Q(j,j)-d'*(A*d)+dU*dU');
        p4=4*(-d'*W(:,j)+dU*Q(:,j));

        % Compute the step size
        lamb=Stepsize(p1,p2,p3,p4);

        % Update the j-th column of U using the computed step size
        U(:,j)=U(:,j)+lamb*d;

        % Update the upper triangle of Q (since Q = U^T * U is symmetric)
        Q(j,j+1:end)=Q(j,j+1:end)+lamb*dU(j+1:end);

        % Accumulate the objective change
        df=df+p1/4*lamb^4+p2/3*lamb^3+p3/2*lamb^2+p4*lamb;
    end

    % Update timing, objective value, and optimality gap for the current iteration
    time_vec=[time_vec,toc(t1)+time_vec(end)]; % Cumulative time
    obj_vec=[obj_vec,df+obj_vec(end)]; % Update objective value
    opt_vec=[opt_vec norm(U-max(U-4*(U*(U'*U) - A*U), 0),'inf')]; % Update optimality gap



    % Stop the algorithm if the cumulative time exceeds the maximum allowed time
    % Or the difference of optimality gap less than 1e-12
    % Or the relative decrease of the objctive functoin less than 1e-12
    if time_vec(end)>=maxtime || abs(opt_vec(end)-opt_vec(end-1))<1e-12 ...
            || abs(df/obj_vec(end-1))<1e-12
        break;
    end
end

% Normalize the objective values with respect to the Frobenius norm of A
obj_vec=sqrt(obj_vec./tempA);
end
