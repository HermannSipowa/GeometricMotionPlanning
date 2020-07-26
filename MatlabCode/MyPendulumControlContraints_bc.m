function [pl,ql,pr,qr] = MyPendulumControlContraints_bc(xl,ul,xr,ur,t) % Boundary condition for GHF
global X0 Xf
n  = size(X0);
pl = ul - X0;
ql = zeros(n);
pr = ur - Xf;
qr = zeros(n);
end
% --------------------------------------------------------------------------