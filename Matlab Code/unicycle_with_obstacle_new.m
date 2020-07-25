close all
clear all
clc
start_up
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-20);

tmax=4;                                     % Integration time for PDE
xmax=1;
tpoints=500;
xpoints=500;

m = 0;
x = linspace(0,xmax,xpoints);                  % discretization of the curve
t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of time interval
global z1 z2 rr % position of the obstacles in R2, rr= radius of obstacles

z1=[0.5;0];z2=[-0.5;0];
rr=.03;

sol = pdepe(m,@mypdexpde,@mypdexic,@mypdexbc,x,t);  % Solve GHF
dataPath = 'Outputs_Data/';
save([dataPath,'sol','.mat'], 'sol','-v7.3')

%% PDE, initial condition and boundary condition functions used
function [c,f,s] = mypdexpde(x,t,u,DuDx)
% Define PDE; Evaluate right-hand-side of GHF


k=180;     % Penalty in the constrained directions


global z1 z2 rr
R = rr+.15;              % Safety zone around obstacle
l1=norm(u(1:2)-z1);
l2=norm(u(1:2)-z2);      % Barrier function lambda


if l1>=R
    lambda1=0;
    plambda1=[0;0;0];
else
    lambda1=((l1^2-R^2)/l1^2)^2;
    plambda1=4*R^2*(l1^2-R^2)/l1^6*[u(1)-z1(1);u(2)-z1(2);0];
end

if l2>=R
    lambda2=0;
    plambda2=[0;0;0];
else
    lambda2=((l2^2-R^2)/l2^2)^2;
    plambda2=4*R^2*(l2^2-R^2)/l2^6*[u(1)-z2(1);u(2)-z2(2);0];
end

lambda=lambda1+lambda2+1;
plambda=plambda1+plambda2;

% Metric
G=[cos(u(3))^2+k*sin(u(3))^2,(1-k)*sin(u(3))*cos(u(3)),0;
    (1-k)*sin(u(3))*cos(u(3)),sin(u(3))^2+k*cos(u(3))^2,0;
    0, 0, 1];


% Partial derivatives of G
pG=zeros(3);
pG(:,:,2)=zeros(3);
pG(:,:,3)=[(k-1)*sin(2*u(3)),(1-k)*cos(2*u(3)),0;
    (1-k)*cos(2*u(3)),(1-k)*sin(2*u(3)),0;
    0,0,0];


pG(:,:,1)=plambda(1)*G;
pG(:,:,2)=plambda(2)*G;
pG(:,:,3)=lambda*pG(:,:,3);

invG=inv(lambda*G);
% Evaluate christoffels symboles
Chris=zeros(size(pG));
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                Chris(i,j,k)=Chris(i,j,k)+...
                    0.5*invG(i,l)*(pG(l,j,k)+pG(l,k,j)-pG(j,k,l));
            end
        end
    end
end
c = [1;1;1];
f = DuDx;
s = DuDx'*squeeze(Chris(1,:,:))*DuDx*[1;0;0]+...
    DuDx'*squeeze(Chris(2,:,:))*DuDx*[0;1;0]+...
    DuDx'*squeeze(Chris(3,:,:))*DuDx*[0;0;1];

end
% --------------------------------------------------------------------------

function u0 = mypdexic(x)  % Initial condition for GHF

u0=[2*x-1;1*sin(2*pi*x);0];

end
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t) % Boundary condition for GHF

pl = ul-[-1;0;0];
ql = [0;0;0];
pr = ur-[1;0;0];
qr = [0;0;0];
end