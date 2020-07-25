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
format short e
start_up
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-20);
global X0 Xf Mass g Penalty tMin tMax rMin rMax ...
       tClearance rClearance z1 z2 rr Alpha ...
       U1_max U2_max tmax xmax 
tmax    = 3;     % Integration time for PDE
xmax    = 5;
tpoints = 10000;    % Number of points in the time discretization
xpoints = 10000;    % Number of points in the spcial discretization
% dx      = 5E-3;
% x       = 0:dx:1;
U1_max  = 100;
U2_max  = 100;
Mass    = 2;      % Mass of the ball
g       = 9.80;   % Gravitational acceleration
Alpha   = 1E4;
Penalty = 1E5;    % Penalty for moving in unfeasible/constrained directions
X0      = [7, -4*pi/9, -0.5, -0.01, 0, 0]'; % Initial state
Xf      = [9, pi/3, 0, 0, 0, 0]'; % Final state
m       = 0; % Integrator parameter (due to cartesian coordinates)
x       = linspace(0,xmax,xpoints); % discretization of the curve
t       = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of time 
tMin    = -pi/2;
tMax    = pi/2;
rMin    = 3;
rMax    = 10;
z1      = [7; -pi/4]; % Location of obstacle 1
z2      = [8; pi/5];  % Location of obstacle 2
rr      = 0.75;       % Radius of the obstacles
c1      = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4      = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('LightSlateGray');
tClearance = pi/36;
rClearance = 0.3;


%% Solving the Geodesic (Parabolic PDE, the solution is of form sol(t,x,u)
tStart = tic;
sol = pdepe(m,@mypendulumpde,@mypdexic,@mypdexbc,x,t); 
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

% u1 = sol(:,:,1); % radius r
% u2 = sol(:,:,2); % angle theta
% u3 = sol(:,:,3); % raduis rate rdot
% u4 = sol(:,:,4); % angular rate thetadot

%%  Extract Controls
% ctrl      = zeros(2,length(x));
% [p1, dp1] = pdeval(m,x,sol(end,:,1),Xmesh);
% [p2, dp2] = pdeval(m,x,sol(end,:,2),Xmesh);
% [p3, dp3] = pdeval(m,x,sol(end,:,3),Xmesh);
% [p4, dp4] = pdeval(m,x,sol(end,:,4),Xmesh);
% [ctrl(1,:), dp5] = pdeval(m,x,sol(end,:,4),Xmesh);
% [ctrl(2,:), dp6] = pdeval(m,x,sol(end,:,4),Xmesh);

% for i = 1:length(x)
%     r         = p1(i);
%     theta     = p2(i);
%     rdot      = p3(i);
%     thetadot  = p4(i);
%     Fbar      = [1 0 0 0; 0 1 0 0; 0 0 1/Mass 0; 0 0 0 1/(r*Mass)];
%     Xdot      = [dp1(i) dp2(i) dp3(i) dp4(i)]';
%     fd        = [rdot; thetadot; r*thetadot^2 + g*cos(theta)/Mass;
%                  -2*rdot*thetadot/r + g*sin(theta)/(Mass*r)];
%     ctrl(:,i) = [zeros(2) eye(2)]/Fbar*(Xdot-fd);
% end


%% Find x^* by integrating ODE
ctrl(1,:) = sol(end,:,5); % linear control
ctrl(2,:) = sol(end,:,6); % angular control
n          = 4;
xstar      = zeros(n,length(x));
xstar(:,1) = X0(1:4);
F_xy       = @(t,X,U) Dynamics(t,X,U);

for i = 1:(length(x)-1)
    dx = x(i+1) - x(i);
    k_1 = F_xy(x(i),xstar(:,i),ctrl(:,i));
    k_2 = F_xy(x(i)+0.5*dx,xstar(:,i)+0.5*dx*k_1,ctrl(:,i));
    k_3 = F_xy((x(i)+0.5*dx),(xstar(:,i)+0.5*dx*k_2),ctrl(:,i));
    k_4 = F_xy((x(i)+dx),(xstar(:,i)+k_3*dx),ctrl(:,i));
    
    xstar(:,i+1) = xstar(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dx;  
end


%% Plotting the results
tpointsCorase = 100;
xpointsCorase = 100;
Tmesh = [0 logspace(-4,log10(tmax),tpointsCorase-1)]; % discretization of time 
Xmesh = linspace(0,xmax,xpointsCorase); % discretization of the curve
[X,Y] = meshgrid(Tmesh,Xmesh);

u1 = interp2(t,x,sol(:,:,1),X,Y); % radius r
u2 = interp2(t,x,sol(:,:,2),X,Y); % angle theta
u3 = interp2(t,x,sol(:,:,3),X,Y); % raduis rate rdot
u4 = interp2(t,x,sol(:,:,4),X,Y); % angular rate thetadot

figure
surf(X,Y,u1)
title('u1(t,s)')
xlabel('Time t')
ylabel('Homotopy s')

figure
surf(X,Y,u2)
title('u2(t,s)')
xlabel('Time t')
ylabel('Homotopy s')

figure
surf(X,Y,u3)
title('u3(t,s)')
xlabel('Time t')
ylabel('Homotopy s')

figure
surf(X,Y,u4)
title('u4(t,s)')
xlabel('Time t')
ylabel('Homotopy s')

%%
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
p   = plot(x_0,y_0);
p.Color = 'magenta';
p.Marker = '*';
p.MarkerSize = 8;

% Plotting the end condition
%-------------------------------------------
rf  = Xf(1); thetaf = Xf(2);
x_f = rf.*cos(thetaf);
y_f = rf.*sin(thetaf);
p   = plot(x_f,y_f);
p.Color = c4;
p.Marker = '+';
p.MarkerSize = 8;

%%
% Show homotopy in 2D with obstacles
%-------------------------------------------
r = u1(1,:); theta = u2(1,:);
X_path = r.*cos(theta);
Y_path = r.*sin(theta);
h1 = plot(X_path,Y_path,'k:.','LineWidth',2);
pause;

for i=1:tpointsCorase  
    
    r = u1(i,:); theta = u2(i,:);
    X_path = r.*cos(theta); h1.XDataSource = 'X_path';
    Y_path = r.*sin(theta); h1.YDataSource = 'Y_path';
    
    refreshdata(h1,'caller');
    drawnow;
end

%%
% Show homotopy in 2D with obstacles
%-------------------------------------------
r = xstar(1,1); theta = xstar(2,1);
X_traj = r.*cos(theta);
Y_traj = r.*sin(theta);
h2 = plot(X_traj,Y_traj,'*','Color',c1,'LineWidth',2);
pause;

for i=1:length(x)
    
    r = xstar(1,i); theta = xstar(2,i);
    X_traj = r.*cos(theta); h2.XDataSource = 'X_traj';
    Y_traj = r.*sin(theta); h2.YDataSource = 'Y_traj';
    
    refreshdata(h2,'caller');
    drawnow;
end

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

function u0 = mypdexic(x) % Initial condition for GHF
global X0 Xf U1_max U2_max xmax
u0 = [X0(1) + (Xf(1)-X0(1))*sin(5*pi*x/(2*xmax))
      X0(2) + (Xf(2)-X0(2))*sin(pi*x/(2*xmax));
      X0(3)*cos(pi*x/(2*xmax)) + Xf(3)*sin(pi*x/(2*xmax));
      X0(4)*cos(pi*x/(2*xmax)) + Xf(4)*sin(pi*x/(2*xmax));
      sin(x*pi/xmax);
      sin(x*pi/xmax)];

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
global Mass g
r        = X(1);
theta    = X(2);
rdot     = X(3);
thetadot = X(4);
fd = [rdot; thetadot; r*thetadot^2 + g*cos(theta)/Mass;
      -2*rdot*thetadot/r + g*sin(theta)/(Mass*r)];
f1 = [0; 0; 1/Mass; 0];
f2 = [0; 0; 0; 1/(r*Mass)];

xdot = fd + f1*U(1) + f2*U(2);

end
