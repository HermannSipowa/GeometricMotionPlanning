function [c,f,s] = mypdexpde(x,t,u,DuDx,k1,k2,k3,N) % Define PDE; right-hand-side of AGHF

XChief = full([Chief_x(x);  Chief_y(x);  Chief_z(x);...
                 Chief_dx(x); Chief_dy(x); Chief_dz(x)]);
LL = full(L(u,DuDx,k1,k2,k3,x,XChief));
CasADiresult = full(dLdx(u,DuDx,k1,k2,k3,x,XChief,LL))';
CasADir = [CasADiresult(1:13), CasADiresult(14:26)]; % [dL/dx, dL/dxdot]

f = CasADir(:,2);  % dL/dxdot
s = -CasADir(:,1); % -dL/dx
c = ones(N,1);

end