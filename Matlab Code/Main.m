%% Example 1: 2D Surface in 3D Space
%-------------------------------------------
close all
clear all
clc
start_up

c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('LightSlateGray');

Smap  = 10;
Straj = 15;
basis = 20;
xgrid = linspace(-Smap, Smap, 100);
ygrid = linspace(-Smap, Smap, 100);
[xi, yi] = meshgrid(xgrid, ygrid);

a = .01; b = -.01; c = .03; d = -.02; e = .04; f = 0;
zi = a.*xi.^2 + b.*yi.^2 + c.*xi.*yi + d.*xi + e.*yi + f;
zfunc = @(x,y) a.*x.^2 + b.*y.^2 + c.*x.*y + d.*x + e.*y + f;
zgrid = zfunc(xgrid,ygrid);

%% Solve BVP
%-------------------------------------------
X0 = [-8,-7,0,0]';
n = 4;
Smap    = 3;
opts    = bvpset('RelTol',2.53447e-07,'AbsTol',2.22045e-20,'Stats','on'); % 'FJacobian',@jac,
xmesh   = linspace(-Smap, Smap, 10);
solinit = bvpinit(xmesh, X0);
sol5c   = bvp5c(@Pendulum_ode, @Pendulum_bcfcn, solinit, opts);



%% Ploting the time history of the system's states
%-------------------------------------------
ex_func = @(x,y) [1, 0, 2.*a.*x + c.*y + d] / norm([1, 0, 2.*a.*x + c.*y + d]);
ey_func = @(x,y) [0, 1, 2.*b.*y + c.*x + e] / norm([0, 1, 2.*b.*y + c.*x + e]);


X = X0(1);
Y = X0(2);
Z = zfunc(X,Y);
ex = ex_func(X,Y);
ey = ey_func(X,Y);

surf(xi,yi,zi, 'EdgeAlpha', 0)
hold on

plt1 = scatter3(sol5c.y(1,1),sol5c.y(2,1),zfunc(sol5c.y(1,1),sol5c.y(2,1)),80);
plt1.MarkerEdgeColor = c2;
plt1.MarkerFaceColor = c2;
plt1.Marker = 'o';

plt2 = scatter3(sol5c.y(1,end),sol5c.y(2,end),zfunc(sol5c.y(1,end),sol5c.y(2,end)),80);
plt2.MarkerEdgeColor = c3;
plt2.MarkerFaceColor = c3;
plt2.Marker = 'o';

plt3 = plot3(sol5c.y(1,:),sol5c.y(2,:),zfunc(sol5c.y(1,:),sol5c.y(2,:)),'c');

plt4 = quiver3(X,Y,Z,ex(1),ex(2),ex(3),2);
plt4.LineWidth = 1;
plt4.Color = 'm';
hold on
plt5 = quiver3(X,Y,Z,ey(1),ey(2),ey(3),2);
plt5.LineWidth = 1;
plt5.Color = 'k';

r_unit = [X,Y,Z]/norm([X,Y,Z]);

colormap copper
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'view',[58 23])
axis equal
legend([plt1,plt2,plt3,plt4,plt5],{'Start','Goal','traj','$\hat{\textbf{\textit{u}}}$','$\hat{\textbf{\textit{v}}}$'})

%% Ploting the time history of the system's states
%-------------------------------------------
figure
YLabel = {'$r(t)$','$\theta(t)$','$\dot{r}(t)$','$\dot{\theta}(t)$'};
for i = 1:n
subplot(2,2,i)
plot(sol5c.x,sol5c.y(i,:),'r')
xlabel('t')
ylabel(YLabel(i))
grid on
end
sgt = sgtitle('Example Toy Problem');
sgt.FontSize = 20;

%% Example 2: Pendulum
%-------------------------------------------

%% Solve BVP
%-------------------------------------------
Smap    = 3;
n       = 4;
opts    = bvpset('FJacobian',@jac,'RelTol',2.53447e-07,'AbsTol',2.22045e-20,'Stats','on'); % 'FJacobian',@jac,
xmesh   = linspace(-Smap, Smap, 10);
solinit = bvpinit(xmesh, ones(4,1));
sol5c   = bvp5c(@bvpfcn, @bcfcn, solinit, opts);


%% Generating plot of the manifold
%-------------------------------------------
[r,t] = meshgrid(0:.1:8,0:pi/30:(2*pi));
z     = r;
% Convert to Cartesian
x = r.*cos(t);
y = r.*sin(t);

figure
hold on;
contourf(x,y,z,'edgecolor','none');
% Turn off axes and set square aspect ratio
axis off
axis image
colormap(flipud(cmap(c6,100,30,30)))

% Plotting the generated trajectory
%-------------------------------------------
x = sol5c.y(1,:).*cos(sol5c.y(2,:));
y = sol5c.y(1,:).*sin(sol5c.y(2,:));
g = plot(x,y,'Color',c2);

% Ploting the initial condition
%-------------------------------------------
r0 = 3; theta0 = -pi/4;
x_0 = r0.*cos(theta0);
y_0 = r0.*sin(theta0);
p   = plot(x_0,y_0);
p.Color = 'magenta';
p.Marker = '*';
p.MarkerSize = 8;

% Plotting the end condition
%-------------------------------------------
rf  = 5; thetaf = 0;
x_f = rf.*cos(thetaf);
y_f = rf.*sin(thetaf);
p   = plot(x_f,y_f);
p.Color = c4;
p.Marker = '+';
p.MarkerSize = 8;


%% Ploting the time history of the system's states
%-------------------------------------------
figure
YLabel = {'$r(t)$','$\theta(t)$','$\dot{r}(t)$','$\dot{\theta}(t)$'};
for i = 1:n
subplot(2,2,i)
plot(sol5c.x,sol5c.y(i,:),'r')
xlabel('t')
ylabel(YLabel(i))
grid on
end
sgt = sgtitle('Example Toy Problem');
sgt.FontSize = 20;


%% Used function
%-------------------------------------------
function xdot = bvpfcn(t,X) % equation to solve
r = X(1);
rdot = X(3);
thetadot = X(4);
xdot = [rdot, thetadot, r*thetadot^2, -2*rdot*thetadot/r]';
end 
%-------------------------------------------
function res = bcfcn(ya,yb) % boundary conditions
r0 = 3; theta0 = -pi/4; r_End = 5; theta_End = 0;
Start = [ r0; theta0]; End = [r_End; theta_End]; 

res = [ya(1:2) - Start
       yb(1:2) - End];
end
%-------------------------------------------
function dfdy = jac(t,X)
r = X(1); rdot = X(3); thetadot = X(4);
dfdy = [0   0   1   0
        0   0   0   1
        thetadot^2   0   0   2*thetadot*r
        2*thetadot*r/r^2   0   -2*rdot/r   -2*thetadot/r];
end


%-------------------------------------------
%-------------------------------------------
function xdot = Pendulum_ode(t,X) % equation to solve

u = X(1); v = X(2); udot = X(3); vdot = X(4);
a = .01; b = -.01; c = .03; d = -.02; e = .04;
denominator = 1 + d^2 + e^2 + 2*c*u*(e + 2*(a + b)*v) + 2*d*(2*a*u + c*v)...
                + c^2*(u^2 + v^2) + 4*(a^2*u^2 + b*v*(e + b*v));
            
xdot(1) = udot;
xdot(2) = vdot;
xdot(3) = 2*(d + 2*a*u + c*v)*(a*udot^2 + vdot*(c*udot + b*vdot)) / denominator;
xdot(4) = 2*(e + c*u + 2*b*v)*(a*udot^2 + vdot*(c*udot + b*vdot)) / denominator;

end
%-------------------------------------------
function res = Pendulum_bcfcn(ya,yb) % boundary conditions
x0 = -8; y0 = -7; x_End = 6; y_End = 7;
Start = [ x0; y0]; End = [x_End; y_End]; 

res = [ya(1:2) - Start
       yb(1:2) - End];
end
%-------------------------------------------


%% -------------------------------------------------------------------------
%% -------------------------------------------------------------------------
%% -------------------------------------------------------------------------
%% -------------------------------------------------------------------------
%% -------------------------------------------------------------------------
%% -------------------------------------------------------------------------
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
