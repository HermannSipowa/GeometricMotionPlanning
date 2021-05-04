%% Pendulum Controlled Swing (with gravity and obstacle)
% Procedure:
% 1. Passing initial condition and boundary condition to the PDE
% 2. Solve PDE
% 3. The solution at tmax approximates the steady state solution
% 4. Extract the controls from sol(tpoints,:,:)
% 5. Simulate the motion of unicycle using extracted controls
%--------------------------------------------------------------------------

close all
clear all
clc
start_up
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-20);
global X0 Xf Mass g Penalty tMin tMax rMin rMax ...
       tClearance rClearance z1 z2 rr Alpha Progress tspan
IntTime = 1;
S_Depth = 10;     % Integration time for PDE
S_Int   = 5E3;   % Number of points in the time discretization
T_Int   = 5E3;   % Number of points in the spcial discretization
Mass    = 2;     % Mass of the ball
g       = 9.81;  % Gravitational acceleration
Alpha   = 1;
Penalty = 1E4;   % Penalty for moving in unfeasible/constrained directions
X0      = [7 -4*pi/9 -0.5 -0.01]'; % Initial state
Xf      = [9 pi/3 0 0]'; % Final state
m       = 0; % Integrator parameter (due to cartesian coordinates)
x       = linspace(0,IntTime,T_Int); % discretization of the curve
t       = [0 logspace(-4,log10(S_Depth),S_Int-1)]; % discretization of time 
tMin    = -pi/2;
tMax    = pi/2;
rMin    = 3;
rMax    = 10;
z1      = [7; -pi/4]; % Location of obstacle 1
z2      = [8; pi/5];  % Location of obstacle 2
rr      = 0.75;        % Radius of the obstacles
c1      = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4      = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('LightSlateGray');
tClearance = pi/36;
rClearance = 0.3;
Progress   = 0;
tspan = x;


%% Solving the Geodesic (Parabolic PDE, the solution is of form sol(t,x,u)
tStart = tic;
system('caffeinate -dims &');
sol = pdepe(m,@mypendulumpde,@mypdexic,@mypdexbc,x,t); 
tEnd = toc(tStart);
system('killall caffeinate');
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
u1 = sol(:,:,1); % radius r
u2 = sol(:,:,2); % angle theta
u3 = sol(:,:,3); % raduis rate rdot
u4 = sol(:,:,4); % angular rate thetadot

%%  Extract Controls
ctrl      = zeros(2,length(x));
[p1, dp1] = pdeval(m,x,sol(end,:,1),x);
[p2, dp2] = pdeval(m,x,sol(end,:,2),x);
[p3, dp3] = pdeval(m,x,sol(end,:,3),x);
[p4, dp4] = pdeval(m,x,sol(end,:,4),x);

for i = 1:length(x)
    r         = p1(i);
    theta     = p2(i);
    rdot      = p3(i);
    thetadot  = p4(i);
    Fbar      = [1 0 0 0; 0 1 0 0; 0 0 1/Mass 0; 0 0 0 1/(r*Mass)];
    Xdot      = [dp1(i) dp2(i) dp3(i) dp4(i)]';
    fd        = [rdot; thetadot; r*thetadot^2 + g*cos(theta)/Mass;
                 -2*rdot*thetadot/r + g*sin(theta)/(Mass*r)];
    ctrl(:,i) = [zeros(2) eye(2)]/Fbar*(Xdot-fd);
end


%% Find x^* by integrating ODE
tic
[~, xstar] = ode113(@(t,X)Dynamics(t,X,ctrl),x,X0);
toc

%% Interpolate the results to a coarser grid
tpointsCorase = 500;
xpointsCorase = 1000;
Tmesh = [0 logspace(-4,log10(S_Depth),tpointsCorase-1)]; % discretization of time
Xmesh = linspace(0,IntTime,xpointsCorase); % discretization of the curve
[X,Y] = meshgrid(Xmesh,Tmesh);
[xgrid,tgrid] = meshgrid(x,t);

xstar = interp1(x,xstar,Xmesh);
u1    = interp2(xgrid,tgrid,u1,X,Y,'spline');
u2    = interp2(xgrid,tgrid,u2,X,Y,'spline');
u3    = interp2(xgrid,tgrid,u3,X,Y,'spline');
u4    = interp2(xgrid,tgrid,u4,X,Y,'spline');


%% Plotting the results
% plotting the configuration manifold
%-------------------------------------------
[r,th] = meshgrid(0:.1:11,0:pi/30:(2*pi));
z_coord     = r;
x_coord  = r.*cos(th);
y_coord  = r.*sin(th);

figure
camroll(-90)
hold on;
contourf(x_coord,y_coord,z_coord,'edgecolor','none');
axis off
axis image
colormap(flipud(cmap(c6,100,30,30)))

% Delimiting the boundering of the reacheable subspace
%-------------------------------------------
rf  = linspace(rMin,rMax,50);
x_f = rf.*cos(tMin);
y_f = rf.*sin(tMin);
p   = plot(x_f,y_f,'--');
p.Color = c4;
p.MarkerSize = 8;

rf  = linspace(rMin,rMax,50);
x_f = rf.*cos(tMax);
y_f = rf.*sin(tMax);
p   = plot(x_f,y_f,'--');
p.Color = c4;
p.MarkerSize = 8;

rf  = rMax; angle = linspace(tMin,tMax,50);
x_f = rf.*cos(angle);
y_f = rf.*sin(angle);
p   = plot(x_f,y_f,'--');
p.Color = c4;
p.MarkerSize = 8;

rf  = rMin; angle = linspace(tMin,tMax,50);
x_f = rf.*cos(angle);
y_f = rf.*sin(angle);
p   = plot(x_f,y_f,'--');
p.Color = c4;

% Plotting the obstacles in the configuration space
%-------------------------------------------
j = 0:0.1:2*pi;
plot(z1(1).*cos(z1(2))+rr*cos(j),z1(1).*sin(z1(2))+rr*sin(j),'r:');
plot(z2(1).*cos(z2(2))+rr*cos(j),z2(1).*sin(z2(2))+rr*sin(j),'r:');
scatter(z1(1).*cos(z1(2)),z1(1).*sin(z1(2)),'r','filled');
scatter(z2(1).*cos(z2(2)),z2(1).*sin(z2(2)),'r','filled');

% Ploting the initial condition
%-------------------------------------------
r0 = X0(1); theta0 = X0(2);
x_0 = r0.*cos(theta0);
y_0 = r0.*sin(theta0);
pi   = plot(x_0,y_0);
pi.Color = 'magenta';
pi.Marker = '*';
pi.MarkerSize = 8;

% Plotting the end condition
%-------------------------------------------
rf  = Xf(1); thetaf = Xf(2);
x_f = rf.*cos(thetaf);
y_f = rf.*sin(thetaf);
pf   = plot(x_f,y_f);
pf.Color = c4;
pf.Marker = '+';
pf.MarkerSize = 8;

h = legend([pi,pf],{'Initial conf','Final conf'});
rect = [0.48, 0.46, 0.08, 0.08];
set(h, 'Position', rect)
box on

%%
% Show homotopy in 2D with obstacles
%-------------------------------------------
r = u1(1,:); theta = u2(1,:);
X_path = r.*cos(theta);
Y_path = r.*sin(theta);
h1 = plot(X_path,Y_path,'k:.','LineWidth',2);
h = legend([pi,pf,h1],{'Initial conf','Final conf','Trajectory'});
pause;

for i=1:length(Tmesh)   
    
    r = u1(i,:); theta = u2(i,:);
    X_path = r.*cos(theta); h1.XDataSource = 'X_path';
    Y_path = r.*sin(theta); h1.YDataSource = 'Y_path';
    
    refreshdata(h1,'caller');
    drawnow;
    F(i) = getframe(gcf);
end

video = VideoWriter('Homotopy.avi','MPEG-4');
open(video)
writeVideo(video,F)
close(video)

%%
% Show homotopy in 2D with obstacles
%-------------------------------------------
% r = xstar(1,1); theta = xstar(1,2);
% X_traj = r.*cos(theta);
% Y_traj = r.*sin(theta);
% h2 = plot(X_traj,Y_traj,'*','Color',c1,'LineWidth',2);
% pause;
% 
% for i=1:length(Xmesh)
%     
%     r = xstar(i,1); theta = xstar(i,2);
%     X_traj = r.*cos(theta); h2.XDataSource = 'X_traj';
%     Y_traj = r.*sin(theta); h2.YDataSource = 'Y_traj';
%     
%     refreshdata(h2,'caller');
%     drawnow;
% end

%%
% Plotting the required control
%-------------------------------------------
figure
subplot(2,1,1)
plot(x,ctrl(1,:))
subplot(2,1,2)
plot(x,ctrl(2,:))



% ------------------------------------------------------------------------
% The followings are the PDE, initial condition and boundary condition for
% the Matlab commend pdepe:

function [c,f,s] = mypendulumpde(x,t,u,DuDx)
% Define PDE; Evaluate right-hand-side of GHF

global Mass g Penalty tMin tMax rMin rMax ...
    tClearance rClearance z1 z2 rr Alpha Progress

if t > Progress
    Progress = t
end
r        = u(1);
theta    = u(2);
rdot     = u(3);
thetadot = u(4);

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


if l5>=R
    lambda5  = 0;
    plambda5 = zeros(4,1);
else
    lambda5  = Alpha*((l5^2-R^2)/l5^2)^2;
    plambda5 = Alpha*4*R^2*(l5^2-R^2)/l5^6*[r-z1(1)*cos(theta-z1(2));
                                      r*z1(1)*sin(theta-z1(2)); 0; 0];
end

if l6>=R
    lambda6  = 0;
    plambda6 = zeros(4,1);
else
    lambda6  = Alpha*((l6^2-R^2)/l6^2)^2;
    plambda6 = Alpha*4*R^2*(l6^2-R^2)/l6^6*[r-z2(1)*cos(theta-z2(2));
                                      r*z2(1)*sin(theta-z2(2)); 0; 0];
end

if theta>=Tmin
    lambda1  = 0;
    plambda1 = zeros(4,1);
else
    l1       = theta-Tmin;
    lambda1  = Alpha*((l1^2-tClearance^2)/l1^2)^2;
    plambda1 = Alpha*[0; 2*tClearance^2/(theta-Tmin)^3; 0; 0];
end

if theta<=Tmax
    lambda2  = 0;
    plambda2 = zeros(4,1);
else
    l2       = theta-Tmax;
    lambda2  = Alpha*((l2^2-tClearance^2)/l2^2)^2;
    plambda2 = Alpha*[0; 2*tClearance^2/(theta-Tmax)^3; 0; 0];
end

if r>=Rmin
    lambda3  = 0;
    plambda3 = zeros(4,1);
else
    l3       = r-Rmin;
    lambda3  = Alpha*((l3^2-rClearance^2)/l3^2)^2;
    plambda3 = Alpha*[2*rClearance^2/(r-Rmin)^3; 0; 0; 0];
end

if r<=Rmax
    lambda4  = 0;
    plambda4 = zeros(4,1);
else
    l4       = r-Rmax;
    lambda4  = Alpha*((l4^2-rClearance^2)/l4^2)^2;
    plambda4 = Alpha*[2*rClearance^2/(r-Rmax)^3; 0; 0; 0];
end

lambda  = 1 +lambda1+lambda2+lambda3+lambda4+lambda5+lambda6;
plambda = plambda1+plambda2+plambda3+plambda4+plambda5+plambda6; 


%% Definint the Metric
H = [Penalty 0 0 0;
     0 Penalty 0 0;
     0 0 Mass^2  0;
     0 0 0 (Mass*r)^2]; % Riemannian metric
G = lambda*H;


%% Partial derivatives of G
pH_pr     = [zeros(3,4); 0, 0, 0, 2*r*Mass^2];
pG(:,:,1) = plambda(1)*H + lambda*pH_pr; % pG_pr
pG(:,:,2) = plambda(2)*H;                % pG_ptheta
pG(:,:,3) = plambda(3)*H;                % pG_prdot
pG(:,:,4) = plambda(4)*H;                % pG_pthetadot


%% Evaluate christoffels symboles
invG  = inv(G);
n     = size(invG,1);
Chris = zeros(size(pG));
for i = 1:n
    for j = 1:n
        for k = 1:n
            for l = 1:n
                Chris(i,j,k) = Chris(i,j,k)...
                               +1/2*invG(i,l)*(pG(l,j,k)...
                               +pG(l,k,j)-pG(j,k,l));
            end
        end
    end
end


%% fD and its partal derivative
Fd = [rdot; thetadot; r*thetadot^2 + g*cos(theta);
      -2*rdot*thetadot/r + g*sin(theta)/r];
pFd = [0, 0, 1, 0;
       0, 0, 0, 1;
       thetadot^2, -g*sin(theta), 0, 2*r*thetadot;
       2*rdot*thetadot/r^2 - g*sin(theta)/r^2,...
       g*cos(theta)/r, -2*thetadot/r, -2*rdot/r];


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

function u0 = mypdexic(x) % Initial condition for GHF
global X0 Xf
u0 = [X0(1) + (Xf(1)-X0(1))*sin(5*pi*x/2)
      X0(2) + (Xf(2)-X0(2))*sin(pi*x/2);
      X0(3)*cos(pi*x/2) + Xf(3)*sin(pi*x/2);
      X0(4)*cos(pi*x/2) + Xf(4)*sin(pi*x/2)];

end
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t) % Boundary condition for GHF
global X0 Xf
n  = size(X0);
pl = ul - X0;
ql = zeros(n);
pr = ur - Xf;
qr = zeros(n);
end
% --------------------------------------------------------------------------

function [xdot] = Dynamics(t,X,U)
global Mass g tspan
idx = find(tspan <= t, 1, 'last' );
r        = X(1);
theta    = X(2);
rdot     = X(3);
thetadot = X(4);
fd = [rdot; thetadot; r*thetadot^2 + g*cos(theta);
    -2*rdot*thetadot/r + g*sin(theta)/r];
f1 = [0; 0; 1/Mass; 0];
f2 = [0; 0; 0; 1/(r*Mass)];

xdot = fd + f1*U(1,idx) + f2*U(2,idx);

end