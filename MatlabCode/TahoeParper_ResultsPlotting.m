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
Statestar_f = rf.*sin(tMin);
p   = plot(x_f,Statestar_f,'--');
p.Color = c4;
p.MarkerSize = 8;

rf  = linspace(rMin,rMax,50);
x_f = rf.*cos(tMax);
Statestar_f = rf.*sin(tMax);
p   = plot(x_f,Statestar_f,'--');
p.Color = c4;
p.MarkerSize = 8;

rf  = rMax; angle = linspace(tMin,tMax,50);
x_f = rf.*cos(angle);
Statestar_f = rf.*sin(angle);
p   = plot(x_f,Statestar_f,'--');
p.Color = c4;
p.MarkerSize = 8;

rf  = rMin; angle = linspace(tMin,tMax,50);
x_f = rf.*cos(angle);
Statestar_f = rf.*sin(angle);
p   = plot(x_f,Statestar_f,'--');
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
Statestar_f = rf.*sin(thetaf);
p   = plot(x_f,Statestar_f);
p.Color = c4;
p.Marker = '+';
p.MarkerSize = 8;

%%
% Show homotopy in 2D with obstacles
%-------------------------------------------
r = xrel(1,:); theta = yrel(1,:);
X_path = r.*cos(theta);
Y_path = r.*sin(theta);
h1 = plot(X_path,Y_path,'k:.','LineWidth',2);
pause;

for i=1:tpoints
    
    r = xrel(i,:); theta = yrel(i,:);
    X_path = r.*cos(theta); h1.XDataSource = 'X_path';
    Y_path = r.*sin(theta); h1.YDataSource = 'Y_path';
    
    refreshdata(h1,'caller');
    drawnow;
end

%%
% Show homotopy in 2D with obstacles
%-------------------------------------------
r = Statestar(1,1); theta = Statestar(2,1);
X_traj = r.*cos(theta);
Y_traj = r.*sin(theta);
h2 = plot(X_traj,Y_traj,'*','Color',c1,'LineWidth',2);
pause;

for i=1:length(time)
    
    r = Statestar(1,i); theta = Statestar(2,i);
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
plot(time,Tau_ctrl(1,:))
subplot(2,1,2)
plot(time,Tau_ctrl(2,:))
