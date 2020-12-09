function [xdot] = TH_ode(t,x,u)
format long 
% global e_chief
% % Reshaping the input into an 6 by n matrix
% [l,~] = size(x);
% x = reshape(x,6,l/6);
% 
% % Computing the CW control and dynamic matrices
% B = [zeros(3); eye(3)];
% k = 1+e_chief*cos(t);
% A = [zeros(3) eye(3);
%     3/k  0  0   0   2  0;
%     0    0  0   -2  0  0;
%     0    0  -1  0   0  0];
% % Computing the derivatives for the ode
% xdot = A*x + B*u;
% xdot = xdot(:);


global e_chief
k = 1+e_chief*cos(t);
A = [zeros(3) eye(3);
    3/k  0  0   0   2  0;
    0    0  0   -2  0  0;
    0    0  -1  0   0  0];
B = [zeros(3); eye(3)];
[l,~] = size(x);
x = reshape(x,6,l/6);
xdot = A*x + B*u;
xdot = xdot(:);
end

