function [xdot] = CW_ode(t,x,u,n)
format long 
% Reshaping the input into an 6 by n matrix
[l,~] = size(x);
x = reshape(x,6,l/6);

% Computing the CW control and dynamic matrices
B = [zeros(3); eye(3)];
A = [zeros(3) eye(3);
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];

% Computing the derivatives for the ode
xdot = A*x + B*u;
xdot = xdot(:);

end