%% Pendulum Controlled Swing (with gravity and obstacle)
% Procedure:
% 1. Passing initial condition and boundary condition to the PDE
% 2. Solve PDE
% 3. The solution at tmax approximates the steady state solution
% 4. Extract the controls from sol(tpoints,:,:)
% 5. Simulate the motion of unicycle using extracted controls
%--------------------------------------------------------------------------
clear all
close all
clc
global Xrel0 Xrelf Penalty ae mu_Earth JD AU Target Chasser X_chief IntTime tspan
tmax     = 100;    % Integration time for PDE
tpoints  = 2000;   % Number of points in the time discretization
m        = 0;      % Integrator parameter (due to cartesian coordinates)\
space    = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of the curve
Penalty  = 1E4;    % Penalty for moving in unfeasible/constrained directions
ae       = 6378.136;
mu_Earth = 3.986004415E5;
JD       = 2456296.25;
AU       = 149597870.7;
% Defining the spacecraft physical characteristic (Eefer Spacecraft class definition)
Chasser  = Spacecraft([30 30 4 116 40 0.9 0.9 3000 6000]); 
Target   = Spacecraft([10 6 2 40 10 0.2 0.5 6 12]);       

%% Inital condition in Orbital Elements

% Target initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a1      = 1.5E+4; % Semi-major axis in Km
e1      = 0.0;    % Eccentricity
inc1    = 45;     % Inclination in deg0
BigOmg1 = 20;     % RAAN in deg
LitOmg1 = 90;     % AOP in deg
M1      = 0;      % Mean anomaly in deg

% Chaser initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a2      = a1+10;               % Semi-major axis in Km
e2      = e1+1e-3;              % Eccentricity
inc2    = inc1+2e-3;           % Inclination in deg
BigOmg2 = BigOmg1;              % RAAN in deg
LitOmg2 = LitOmg1;              % AOP in deg
M2      = M1;                   % True anomaly in rad

% Desired trajectory initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a3      = a1;           % Semi-major axis in Km
e3      = e1+0.4/a1;   % Eccentricity
inc3    = inc1;         % Inclination in deg
BigOmg3 = BigOmg1;      % RAAN in deg
LitOmg3 = LitOmg1;      % AOP in deg
M3      = M1;           % True anomaly in rad


%% Inital condition in coratesian corditnates

% integrting the chief/reference trajectory
Period  = 2*pi*sqrt(a1^3/mu_Earth);
IntTime = Period;
tspan   = linspace(0,IntTime,10000);
options = odeset('RelTol',2.22045e-13,'AbsTol',2.22045e-30);

% Constructing the chief inital conditions from the chief's OEs
inc1 = deg2rad(inc1); BigOmg1 = deg2rad(BigOmg1); LitOmg1 = deg2rad(LitOmg1); M1 = deg2rad(M1);
COE_c = [a1,e1,inc1,BigOmg1,LitOmg1,M1];
[Position_target,Velocity_target]  = COEstoRV_MeanAnonaly(COE_c,mu_Earth,1);
X0_chief = [Position_target; Velocity_target];
[~, X_chief] = ode113(@(t,X)ChiefMotionODE(t,X,Target),tspan,X0_chief,options);

TN = DCM(Position_target,Velocity_target); rt_norm = norm(Position_target);
h_vec = cross(Position_target,Velocity_target); h_norm = norm(h_vec);
eh = h_vec/h_norm;
U_eci_Target  = F_CanonBall(tspan(1),Position_target,Target); % SRP force on the Target
N_nudot = h_vec/rt_norm^2 + dot(U_eci_Target,eh)*Position_target/h_norm;

% Constructing the deputy inital conditions from the deputy's OEs
inc2 = deg2rad(inc2); BigOmg2 = deg2rad(BigOmg2); LitOmg2 = deg2rad(LitOmg2); M2 = deg2rad(M2);
COE_d = [a2,e2,inc2,BigOmg2,LitOmg2,M2];
[Position_chaser0,Velocity_chaser0]  = COEstoRV_MeanAnonaly(COE_d,mu_Earth,1);
qr_chasser0 = [1 1 1 9]';
qr_chasser0 =  qr_chasser0/norm(qr_chasser0);
Omega_chaser0 = 0.2*[deg2rad(-.5) deg2rad(.4) deg2rad(.1)]'; % Angular velocity in inertial frame
CN = Quaternion_to_DCM(qr_chasser0); Omega_chaser_B0 = CN*Omega_chaser0;
NR_rel0 = Position_chaser0 - Position_target; NV_rel0 = Velocity_chaser0 - Velocity_target;
TR_rel0 = TN*NR_rel0; TV_rel0 = TN*(NV_rel0 - cross(N_nudot,NR_rel0));
Xrel0 = [qr_chasser0; Omega_chaser_B0; TR_rel0; TV_rel0];

% Constructing the deputy final conditions from the deputy's OEs
inc3 = deg2rad(inc3) + .6/a1; BigOmg3 = deg2rad(BigOmg3); LitOmg3 = deg2rad(LitOmg3); M3 = deg2rad(M3);
COE_r = [a3,e3,inc3,BigOmg3,LitOmg3,M3];
[Position_Desired,Velocity_Desired]  = COEstoRV_MeanAnonaly(COE_r,mu_Earth,1);
qr_chasserf = [1 0 0 0]';
qr_chasserf = qr_chasserf/norm(qr_chasserf);
Omega_chaserf = zeros(3,1);
CN = Quaternion_to_DCM(qr_chasserf); Omega_chaser_Bf = CN*Omega_chaserf;
NR_relf = Position_Desired - Position_target; NV_relf = Velocity_Desired - Velocity_target;
TR_relf = TN*NR_relf; TV_relf = TN*(NV_relf - cross(N_nudot,NR_relf));
Xrelf = [qr_chasserf; Omega_chaser_Bf; TR_relf; TV_relf];

%% Solving the Geodesic (Parabolic PDE, the solution is of form sol(t,x,u)
tic
sol = pdepe(m,@TahoePDE,@TahoeIC,@TahoeBC,tspan,space);
toc

%%  Extract Controls
n = length(Xrel0);
nn = length(tspan);
Tau_ctrl  = zeros(3,nn);
X = zeros(n,nn);
Xdot = zeros(n,nn);
for i = 1:n
    [X(i,:), Xdot(i,:)] = pdeval(m,tspan,sol(end,:,i),tspan);
end
for i = 1:nn
    Iinv = inv(diag(Chasser.Moment_Of_Inertia_Calculator()));
    fd   = RelativeMotionODE(i,X(:,i),Target,Chasser);
    Fbar = [zeros(4,3);Iinv;zeros(6,3)];
    Tau_ctrl(:,i) = (Fbar'*Fbar)\Fbar'*(Xdot(:,i)-fd);
end

%% Find x^* by integrating ODE
tic
[~, Xstar_f] = ode113(@(t,Aug_X)ControlledRelativeMotionODE(t,Aug_X,Tau_ctrl,Target,Chasser),tspan,Xrel0,options);
toc
%% Saving results
-------------------------------------------
dataPath = 'Outputs_Data/';
tpointsCorase = 500;
xpointsCorase = 500;
Tmesh = [0 logspace(-4,log10(tmax),tpointsCorase-1)]; % discretization of time
Xmesh = linspace(0,xmax,xpointsCorase); % discretization of the curve
[X,Y] = meshgrid(Tmesh,Xmesh);

for i = 1:length(X0)
    data = sprintf('State%d', i);
    udata = interp2(t,x,sol(:,:,i),X,Y);
    save([dataPath,data,'.mat'], 'udata','-v7.3')
end
Xstar = interp1(x,Xstar_f',X)';
save([dataPath,'Xstar','.mat'], 'Xstar','-v7.3')
save([dataPath,data,'.mat'], 'u1','-v7.3')

%% Visualizing the results
%-------------------------------------------
for k = 1
%     close all
%     clc
%     start_up
%     c1 = rgb('Black'); c2 = rgb('Cyan'); c3 = rgb('Goldenrod');
%     c4 = rgb('RosyBrown'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
%     
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~% sol(end,:,i)
%     % Plotting the 3D homotopies
%     figure('Renderer', 'painters', 'Position', [15 15 900 600])
%     h1 = plot3(sol(1,:,8),sol(1,:,9),sol(1,:,10),'r');
%     hold on
%     h2 = plot3(sol(1,1,8),sol(1,1,9),sol(1,1,10),'o','color', c2,...
%         'LineWidth',2,...
%         'MarkerEdgeColor',c2,...
%         'MarkerFaceColor',c2,...
%         'MarkerSize',5);
%     h3 = plot3(sol(1,end,8),sol(1,end,9),sol(1,end,10),'go',...
%         'LineWidth',2,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',c1,...
%         'MarkerSize',5);
%     h5 = plot3(0,0,0,'bo',...
%         'LineWidth',2,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',c5,...
%         'MarkerSize',15);
%     %arrow3D([.5e4,.5e4,.5e4] ,8e3*XSUN_LVLH,c3)
%     h = legend([h1, h5, h2, h3],{'Relative','chief','$\delta X_0$','$\delta X_f$'});
%     xlabel('X [km]')
%     ylabel('Y [km]')
%     zlabel('Z [km]')
%     for i=1:tpoints
%         X = sol(i,:,8) ;h1.XDataSource='X';
%         Y = sol(i,:,9) ;h1.YDataSource='Y';
%         Z = sol(i,:,10) ;h1.ZDataSource='Z';
%         
%         refreshdata(h1,'caller');
%         drawnow;
%         
%     end
%     
%     
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~%
%     % Plotting the 3D trajectory
%     for ii =1
%         c1 = rgb('Crimson'); c2 = rgb('DarkSlateGray');
%         figure('Renderer', 'painters', 'Position', [15 15 900 600])
%         for kk = 1:4
%             subplot(2,2,kk)
%             for jj=1
%                 h1 = plot3(Xstar_f(:,8),Xstar_f(:,9),Xstar_f(:,10),'r');
%                 hold on
%                 h2 = plot3(Xstar_f(1,8),Xstar_f(1,9),Xstar_f(1,10),'ko','color', c3,...
%                     'LineWidth',2,...
%                     'MarkerEdgeColor','b',...
%                     'MarkerFaceColor',c3,...
%                     'MarkerSize',5);
%                 h3 = plot3(Xstar_f(end,8),Xstar_f(end,9),Xstar_f(end,10),'go',...
%                     'LineWidth',2,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor',c1,...
%                     'MarkerSize',5);
%                 h5 = plot3(0,0,0,'bo',...
%                     'LineWidth',2,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor',c2,...
%                     'MarkerSize',15);
%                 % arrow3D([.5e4,.5e4,.5e4] ,8e3*XSUN_LVLH,c3)
%                 h = legend([h1, h5, h2, h3],{'Relative','chief','$\delta X_0$','$\delta X_f$'});
%                 rect = [0.43, 0.435, 0.125, 0.12];
%                 set(h, 'Position', rect)
%                 box on
%                 grid on
%                 xlabel('X [km]')
%                 ylabel('Y [km]')
%                 zlabel('Z [km]')
%             end
%             
%             if kk == 1
%                 view(70,10)
%             elseif kk==2
%                 view(0,90)  % XY
%             elseif kk==3
%                 view(0,0)   % XZ
%             else
%                 view(90,0)  % YZ
%             end
%             
%         end
%         sgt = sgtitle('Relative Trajectory','interpreter','tex');
%         sgt.FontSize = 30;
%     end
%     
%     
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~%
%     % Plotting the difference between ECI Cartesian coordinates of the deputy
%     for jj=1
%         
%         DQ1 = [Xstar_f(:,8) Xstar_f(:,11) Xstar_f(:,9) Xstar_f(:,12) Xstar_f(:,10) Xstar_f(:,13)];
%         figure('Renderer', 'painters', 'Position', [15 15 900 600])
%         YLabel={'$\delta x [km]$',...
%             '$\delta \dot{x} [km/s]$',...
%             '$\delta y [km]$',...
%             '$\delta \dot{y} [km/s]$',...
%             '$\delta z [km]$',...
%             '$\delta \dot{z} [km/s]$'};
%         for i=1:6
%             subplot(3,2,i)
%             plot1 = plot(tspan/86400,DQ1(:,i),'color', c4);
%             ylabel(YLabel(i))
%             xlabel('days')
%             grid on
%         end
%         sgt = sgtitle('Relative Position','interpreter','tex');
%         sgt.FontSize = 30;
%     end
%     
%     
%     for jj=1
%         figure('Renderer', 'painters', 'Position', [15 15 900 600])
%         YLabel={'$\tau_{_{\!1}}$','$\tau_{_{\!2}}$','$\tau_{_{\!3}}$'};
%         for i=1:3
%             subplot(3,2,i)
%             plot(tspan/86400,Tau_ctrl(i,:),'color', c6)
%             ylabel(YLabel(i))
%             xlabel('days')
%             grid on
%             % legend('Classical Dynamics', 'Dual Quaternions')
%         end
%         sgt = sgtitle('Acceleration and Torque on the Flat-Plate','interpreter','tex');
%         sgt.FontSize = 30;
%     end
%     
%     
%     for jj=1
%         Q1 = Xstar_f(:,1:4);
%         figure('Renderer', 'painters', 'Position', [15 15 900 600])
%         YLabel={'$q_{_{\!1}}$','$q_{_{\!2}}$','$q_{_{\!3}}$','$q_{_{\!4}}$'};
%         for i=1:4
%             subplot(2,2,i)
%             plot(tspan/86400,Q1(:,i),'color', c5) % -Q2(:,i)
%             ylabel(YLabel(i))
%             xlabel('days')
%             grid on
%             % legend('Classical Dynamics', 'Dual Quaternions')
%         end
%         sgt = sgtitle('Deputy Attitude (Quaternions)','interpreter','tex');
%         sgt.FontSize = 30;
%     end
%     
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~%
%     % Plotting the angular velocy
%     for jj=1
%         Q2 = Xstar_f(:,5:7);
%         figure('Renderer', 'painters', 'Position', [15 15 900 600])
%         YLabel={'$\omega_{_{\!1}}$','$\omega_{_{\!2}}$','$\omega_{_{\!3}}$'};
%         for i=1:3
%             subplot(3,1,i)
%             % plot(time/86400,omega_chasser1(:,i)-omega_target1(:,i))
%             plot(tspan/86400,Q2(:,i),'color', c5)
%             ylabel(YLabel(i))
%             xlabel('days')
%             grid on
%             % legend('Classical Dynamics', 'Dual Quaternions')
%         end
%         sgt = sgtitle('Deputy Angular Velocity','interpreter','tex');
%         sgt.FontSize = 30;
%     end
%     
end % Plots for relative motions (Chaser and Target) with classical dynamics

