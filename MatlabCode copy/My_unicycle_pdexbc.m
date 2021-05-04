function [pl,ql,pr,qr] = My_unicycle_pdexbc(xl,ul,xr,ur,t) % Boundary condition for GHF

pl = ul-[-1;0;0];
ql = [0;0;0];
pr = ur-[1;0;0];
qr = [0;0;0];
end