function [c,f,s] = MyPendulumControlContraints_pde(x,t,u,DuDx)
% Define PDE; Evaluate right-hand-side of GHF

global Mass g Penalty tMin tMax rMin rMax ...
       tClearance rClearance z1 z2 rr ...
       Alpha U1_max U2_max
r        = u(1);
theta    = u(2);
rdot     = u(3);
thetadot = u(4);
ctrl1    = u(5);
ctrl2    = u(6);
n        = size(u,1);

%% Defining the Obstacle(s)
Tmin = tMin + tClearance; % Minimum angle
Tmax = tMax - tClearance; % Minimum angle
Rmin = rMin + rClearance; % Minimum radius
Rmax = rMax - rClearance; % Minimum radius
R    = rr + rClearance; 
tip  = [r.*cos(theta); r.*sin(theta)];
obs1 = [z1(1).*cos(z1(2)); z1(1).*sin(z1(2))];
obs2 = [z2(1).*cos(z2(2)); z2(1).*sin(z2(2))];
l5   = norm(tip-obs1);
l6   = norm(tip-obs2); 


if theta>=Tmin
    lambda1  = 0;
    plambda1 = zeros(n,1);
else
    l1       = theta-Tmin;
    lambda1  = Alpha*((l1^2-tClearance^2)/l1^2)^2;
    plambda1 = Alpha*[0; 2*tClearance^2/(theta-Tmin)^3; 0; 0; 0; 0];
end

if theta<=Tmax
    lambda2  = 0;
    plambda2 = zeros(n,1);
else
    l2       = theta-Tmax;
    lambda2  = Alpha*((l2^2-tClearance^2)/l2^2)^2;
    plambda2 = Alpha*[0; 2*tClearance^2/(theta-Tmax)^3; 0; 0; 0; 0];
end

if r>=Rmin
    lambda3  = 0;
    plambda3 = zeros(n,1);
else
    l3       = r-Rmin;
    lambda3  = Alpha*((l3^2-rClearance^2)/l3^2)^2;
    plambda3 = Alpha*[2*rClearance^2/(r-Rmin)^3; 0; 0; 0; 0; 0];
end

if r<=Rmax
    lambda4  = 0;
    plambda4 = zeros(n,1);
else
    l4       = r-Rmax;
    lambda4  = Alpha*((l4^2-rClearance^2)/l4^2)^2;
    plambda4 = Alpha*[2*rClearance^2/(r-Rmax)^3; 0; 0; 0; 0; 0];
end

if l5>=R
    lambda5  = 0;
    plambda5 = zeros(n,1);
else
    lambda5  = Alpha*((l5^2-R^2)/l5^2)^2;
    plambda5 = Alpha*4*R^2*(l5^2-R^2)/l5^6*[r-z1(1)*cos(theta-z1(2));
                                      r*z1(1)*sin(theta-z1(2)); 0; 0; 0; 0];
end

if l6>=R
    lambda6  = 0;
    plambda6 = zeros(n,1);
else
    lambda6  = Alpha*((l6^2-R^2)/l6^2)^2;
    plambda6 = Alpha*4*R^2*(l6^2-R^2)/l6^6*[r-z2(1)*cos(theta-z2(2));
                                      r*z2(1)*sin(theta-z2(2)); 0; 0; 0; 0];
end

if U1_max>=ctrl1
    lambda7  = 0;
    plambda7 = zeros(n,1);
else
    lambda7  = Alpha/(U1_max^2 - ctrl1^2);
    plambda7 = Alpha*[zeros(4,1); 2*ctrl1/(U1_max^2 - ctrl1^2)^2; 0];
end

if U2_max>=ctrl2
    lambda8  = 0;
    plambda8 = zeros(n,1);
else
    lambda8  = Alpha/(U2_max^2 - ctrl2^2);
    plambda8 = Alpha*[zeros(5,1); 2*ctrl2/(U2_max^2 - ctrl2^2)^2];
end

lambda  = 1+lambda7+lambda8+lambda1+lambda2+lambda3+lambda4+lambda5+lambda6;
plambda = plambda7+plambda8+plambda1+plambda2+plambda3+plambda4+plambda5+plambda6; 


%% Definint the Metric
H = diag([Penalty,Penalty,Penalty,Penalty,1,1]); % Riemannian metric
G = lambda*H;


%% Partial derivatives of G
pG(:,:,1) = plambda(1)*H; % pG_pr
pG(:,:,2) = plambda(2)*H; % pG_ptheta
pG(:,:,3) = plambda(3)*H; % pG_prdot
pG(:,:,4) = plambda(4)*H; % pG_pthetadot
pG(:,:,5) = plambda(5)*H; % pG_pctrl1
pG(:,:,6) = plambda(6)*H; % pG_pctrl2


%% Evaluate christoffels symboles
invG  = inv(G);
Chris = zeros(size(pG));
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                Chris(i,j,k) = Chris(i,j,k)...
                               +1/2*invG(i,l)*(pG(l,j,k)...
                               +pG(l,k,j)-pG(j,k,l));
            end
        end
    end
end


%% fD and its partal derivative
Fd = [rdot; thetadot; r*thetadot^2 + g*cos(theta)/Mass + ctrl1/Mass;
      -2*rdot*thetadot/r + (g*sin(theta) + ctrl2)/(Mass*r); 0; 0];
pFd = [0, 0, 1, 0, 0, 0;
       0, 0, 0, 1, 0, 0;
       thetadot^2, -g*sin(theta)/Mass, 0, 2*r*thetadot, 1/Mass, 0;
       2*rdot*thetadot/r^2 - (g*sin(theta) + ctrl2)/(Mass*r^2),...
       g*cos(theta)/(Mass*r), -2*thetadot/r, -2*rdot/r, 0, 1/Mass;
       zeros(2,n)];


%% Computing the coeficients to the PDE
s       = zeros(n,1);
MatProd = zeros(n,1);
for i = 1:n
    MatProd(i) = (DuDx-Fd)'*pG(:,:,i)*Fd;
end
rvec = invG*(pFd'*G*(DuDx-Fd) + 1/2*MatProd);
c = ones(n,1);
f = DuDx-Fd;
for i = 1:n
    s(i) = (DuDx-Fd)'*squeeze(Chris(i,:,:))*(DuDx-Fd)+rvec(i);
end

end
% --------------------------------------------------------------------------