clear 
clc
addpath('/Users/hermann/Desktop/casadi-osx-matlabR2015a-v3.5.5')
import casadi.*

% opts = struct('main', true,...
%               'mex', true);
% %%
% clear
% 
% mex gen.c -largeArrayDims  % Matlab
% disp(gen('f', 3.14))


% X = SX.sym('X',6); mu = SX.sym('mu');
% r = X(1:3); v = X(4:6);
% Xdot = [v;- mu*r/norm(r)^3];
% 
% f = Function('f',...
%     {X,mu},...
%     {Xdot},...
%     {'X','mu'},{'Xdot'});
% f.print_dimensions
% 
% J0 = f.jacobian();
% J0.print_dimensions
% J1 = f.jac();
% J1.print_dimensions

t     = SX.sym('t'); Penalty = SX.sym('Penalty');
rho1  = SX.sym('rho1',6); rho2 = SX.sym('rho2',6);
drho1 = SX.sym('rho1',6); drho2 = SX.sym('rho2',6);
Xaug  = [rho1;rho2];
dXaug = [drho1;drho2];
AgentNum = 2; M = 3; N = 6;
e_chief = 0.25;
k = 1+e_chief*cos(t);
A = [zeros(3) eye(3);
    3/k  0  0   0   2  0;
    0    0  0   -2  0  0;
    0    0  -1  0   0  0];
for i = 1:AgentNum
    
    if i == 1
        AMat_Aug = A;
        % [Fc F], note: the F_bar matrix in paper
        F = eye(N); F_Aug = F;
        % penalty matrix (D matrix in the paper)
        D = diag([Penalty*ones(1,N-M) ones(1,M)]); D_Aug = D;
    else
        AMat_Aug = blkdiag(AMat_Aug,A);
        F_Aug = blkdiag(F_Aug,F);
        D_Aug = blkdiag(D_Aug,D);
    end
end
% drift term f(x) 
f = AMat_Aug*Xaug;
% Metric, before adding the state constraint barrier functions
H = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
% -------------------- state constraint barrier function ------------------
B = [];
b = 0;
B = [B b]; % attach the barrier function to B

% -------------------- Metric and curve Length ----------------------------
% the metric with state constraint barier functions
G = (sum(B)+1)*H;
% actuated curve length
L = Function('L',{Xaug,dXaug,Penalty,t},...
    {(dXaug - f).' * G * (dXaug - f)},...
    {'X','dX','k','t'},...
    {'Curve Length'});

CasadigetL = L.jac();

CasadigetL.print_dimensions






