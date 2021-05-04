close all
clear all
clc
start_up
options = odeset('RelTol',2.22045e-5,'AbsTol',2.22045e-6);

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

tic
sol = pdepe(m,@My_unicycle_pdexpde,@My_unicycle_pdexic,@My_unicycle_pdexbc,x,t);  % Solve GHF
toc

% dataPath = 'Outputs_Data/';
% save([dataPath,'sol','.mat'], 'sol','-v7.3')
% Data = load('Outputs_Data/sol.mat');
% sol = Data.sol;

%% The solution is of form sol(t,x,i)
u1 = sol(:,:,1);    % x position
u2 = sol(:,:,2);    % y position
u3 = sol(:,:,3);    % angle theta

%%
tpointsCorase = 100;
xpointsCorase = 100;
Tmesh = [0 logspace(-4,log10(tmax),tpointsCorase-1)]; % discretization of time 
Xmesh = linspace(0,xmax,xpointsCorase); % discretization of the curve
[X,Y] = meshgrid(Tmesh,Xmesh);
u1Coarse = interp2(t,x,sol(:,:,1),X,Y); % radius r
u2Coarse = interp2(t,x,sol(:,:,2),X,Y); % angle theta
u3Coarse = interp2(t,x,sol(:,:,3),X,Y); % raduis rate rdot

figure
surf(t,x,u1)
hold on
surf(X,Y,u1Coarse)
title('u1(t,s)')
xlabel('Time t')
ylabel('Homotopy s')

figure
surf(t,x,u2)
hold on
surf(X,Y,u2Coarse)
title('u2(t,s)')
xlabel('Time t')
ylabel('Homotopy s')

figure
surf(t,x,u3)
hold on
surf(X,Y,u3Coarse)
title('u3(t,s)')
xlabel('Time t')
ylabel('Homotopy s')
%%
% 3d homotopy plot without plotting obstacles
figure
X=u1(1,:);
Y=u2(1,:);
Z=u3(1,:);
h1=plot3(X,Y,Z,'k','LineWidth',2);
axis([-1.5,1.5,-1.5,1.5,-1.5,1.5]);
grid on;
title('3D configuration space curve');
xlabel('x');
ylabel('y');
zlabel('$\theta$');
hold on;
pause;

for i=1:tpoints
    X=u1(i,:);h1.XDataSource='X';
    Y=u2(i,:);h1.YDataSource='Y';
    Z=u3(i,:);h1.ZDataSource='Z';
    
    refreshdata(h1,'caller');
    drawnow;
end

% Show homotopy in 2D with obstacles
figure
X=u1(1,:);
Y=u2(1,:);
h1=plot(X,Y,'k','LineWidth',2);
hold on;
j=0:0.1:2*pi;
plot(z1(1)+rr*cos(j),z1(2)+rr*sin(j),'r:');
plot(z2(1)+rr*cos(j),z2(2)+rr*sin(j),'r:');
scatter(z1(1),z1(2),'r','filled');
scatter(z2(1),z2(2),'r','filled');
axis([-1.5,1.5,-1.0,1.0]);
grid ON;
title('xy-plane trajectory and unicycle simulation');
xlabel('x');
ylabel('y');
pause;

for i=1:tpoints
    X=u1(i,:);h1.XDataSource='X';
    Y=u2(i,:);h1.YDataSource='Y';
    
    refreshdata(h1,'caller');
    drawnow;
    
end


% Extract Controls
bu        = zeros(2,xpoints);
[p1, dp1] = pdeval(m,x,sol(end,:,1),x);
[p2, dp2] = pdeval(m,x,sol(end,:,2),x);
[p3, dp3] = pdeval(m,x,sol(end,:,3),x);

for i = 1 : xpoints
    projB   = [cos(p3(i)) sin(p3(i)) 0; 0 0 1];
    bu(:,i) = projB*[dp1(i);dp2(i);dp3(i)];
end



% Find x^* by integrating ODE
xstar       =   zeros(3,xpoints);
xstar(1,1)  =   p1(1); %Initial state
xstar(2,1)  =   p2(1);
xstar(3,1)  =   p3(1);

for i = 1: xpoints-1
    cB = [cos(xstar(3,i)) 0;
        sin(xstar(3,i)) 0;
        0 1];
    xstar(:,i+1) = xstar(:,i) + cB*bu(:,i)*(x(i+1)-x(i));
end


% Plot the trajectory x^*
% Negligible error so omitted
% Drawing of unicycle
w=0.3; % width
l=0.1; % height
R=[-w/2 w/2 w/2 -w/2 -w/2; -l/2 -l/2 l/2 l/2 -l/2];


xc=xstar(1,1);yc=xstar(2,1);theta=xstar(3,1);
XY=[cos(theta) -sin(theta);sin(theta) cos(theta)]*R;  %rotate angle alpha
h=patch(XY(1,:)+xc,XY(2,:)+yc,[0.4,0.4,1]);

plot(z1(1)+rr*cos(j),z1(2)+rr*sin(j),'r:');
plot(z2(1)+rr*cos(j),z2(2)+rr*sin(j),'r:');
scatter(z1(1),z1(2),'r','filled');
scatter(z2(1),z2(2),'r','filled');

hold on;
% axis([-0.5,1.5,-0.5,1.5]);
h2=quiver(xc,yc,cos(theta),sin(theta),'r');
set(h2,'AutoScaleFactor',0.15,'MaxHeadSize',2)
grid on
pause();

% Show unicycle simulation
for i=1:xpoints
    xc=sol(tpoints,i,1);yc=sol(tpoints,i,2);
    xc=xstar(1,i);yc=xstar(2,i);theta=xstar(3,i);
    R=[-w/2 w/2 w/2 -w/2 -w/2; -l/2 -l/2 l/2 l/2 -l/2];
    XY=[cos(theta) -sin(theta);sin(theta) cos(theta)]*R;
    h.XData=XY(1,:)+xc;h.YData=XY(2,:)+yc;
    h2.XData=xc;h2.YData=yc;h2.UData=cos(theta);h2.VData=sin(theta);
    refreshdata(h);
    refreshdata(h2);
    drawnow
    
end

% Plotting the required control
figure
subplot(2,1,1)
plot(t,bu(1,:))
subplot(2,1,2)
plot(t,bu(2,:))
% --------------------------------------------------------------------------


% %% PDE, initial condition and boundary condition functions used
% function [c,f,s] = My_unicycle_pdexpde(x,t,u,DuDx)
% % Define PDE; Evaluate right-hand-side of GHF
% 
% 
% k=180;     % Penalty in the constrained directions
% 
% 
% global z1 z2 rr
% R = rr+.15;              % Safety zone around obstacle
% l1=norm(u(1:2)-z1);
% l2=norm(u(1:2)-z2);      % Barrier function lambda
% 
% 
% if l1>=R
%     lambda1=0;
%     plambda1=[0;0;0];
% else
%     lambda1=((l1^2-R^2)/l1^2)^2;
%     plambda1=4*R^2*(l1^2-R^2)/l1^6*[u(1)-z1(1);u(2)-z1(2);0];
% end
% 
% if l2>=R
%     lambda2=0;
%     plambda2=[0;0;0];
% else
%     lambda2=((l2^2-R^2)/l2^2)^2;
%     plambda2=4*R^2*(l2^2-R^2)/l2^6*[u(1)-z2(1);u(2)-z2(2);0];
% end
% 
% lambda=lambda1+lambda2+1;
% plambda=plambda1+plambda2;
% 
% % Metric
% G=[cos(u(3))^2+k*sin(u(3))^2,(1-k)*sin(u(3))*cos(u(3)),0;
%     (1-k)*sin(u(3))*cos(u(3)),sin(u(3))^2+k*cos(u(3))^2,0;
%     0, 0, 1];
% 
% 
% % Partial derivatives of G
% pG=zeros(3);
% pG(:,:,2)=zeros(3);
% pG(:,:,3)=[(k-1)*sin(2*u(3)),(1-k)*cos(2*u(3)),0;
%     (1-k)*cos(2*u(3)),(1-k)*sin(2*u(3)),0;
%     0,0,0];
% 
% 
% pG(:,:,1)=plambda(1)*G;
% pG(:,:,2)=plambda(2)*G;
% pG(:,:,3)=lambda*pG(:,:,3);
% 
% invG=inv(lambda*G);
% % Evaluate christoffels symboles
% Chris=zeros(size(pG));
% for i=1:3
%     for j=1:3
%         for k=1:3
%             for l=1:3
%                 Chris(i,j,k)=Chris(i,j,k)+...
%                     0.5*invG(i,l)*(pG(l,j,k)+pG(l,k,j)-pG(j,k,l));
%             end
%         end
%     end
% end
% c = [1;1;1];
% f = DuDx;
% s = DuDx'*squeeze(Chris(1,:,:))*DuDx*[1;0;0]+...
%     DuDx'*squeeze(Chris(2,:,:))*DuDx*[0;1;0]+...
%     DuDx'*squeeze(Chris(3,:,:))*DuDx*[0;0;1];
% 
% end
% % --------------------------------------------------------------------------
% 
% function u0 = My_unicycle_pdexic(x)  % Initial condition for GHF
% 
% u0=[2*x-1;1*sin(2*pi*x);0];
% 
% end
% % --------------------------------------------------------------------------
% 
% function [pl,ql,pr,qr] = My_unicycle_pdexbc(xl,ul,xr,ur,t) % Boundary condition for GHF
% 
% pl = ul-[-1;0;0];
% ql = [0;0;0];
% pr = ur-[1;0;0];
% qr = [0;0;0];
% end