function [pl,ql,pr,qr] = TahoeBC(xl,ul,xr,ur,t) % Boundary condition for GHF
global Xrel0 Xrelf
n  = size(Xrel0);
pl = ul - Xrel0;
ql = zeros(n);
pr = ur - Xrelf;
qr = zeros(n);
end
% --------------------------------------------------------------------------