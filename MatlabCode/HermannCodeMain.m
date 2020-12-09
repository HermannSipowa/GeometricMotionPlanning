%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hermann Kaptui Sipowa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% clear all
close all
start_up
format long e
% clc
global mu_dot tspan Period PysicalTime

%% Calculating the boundaries conditions
%=======================================%

mm         = 6; % Size of the state space
mu_Earth   = 3.986004415E5; % Earth gravitational parameter
Num_Agent  = 2; % Number of agents in the system
Num_Agent1 = 1; % Number of agents in the first final-orbit
Num_Agent2 = 1; % Number of agents in the second final-orbit

% Specifying the chief's orbit
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a_chief      = 1.42e4;   % Semi-major axis in Km
e_chief      = 0.0;      % Eccentricity
inc_chief    = 50;       % Inclination in deg0
BigOmg_chief = 10;       % RAAN in deg
LitOmg_chief = 10;       % AOP in deg
M_chief      = 0;        % True anomaly in deg


% Computing the nearly non-singular orbital elements
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
LitOmg_chief = deg2rad(LitOmg_chief); M_chief = deg2rad(M_chief);
COE = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief];
[Position_target,Velocity_target]  = COEstoRV(COE,mu_Earth); f_chief = M_chief;
q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
theta_chief = f_chief+LitOmg_chief;

% Specify the initial relative-orbital element of the deputies
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dela = 0;
dele = 1/(2*a_chief);
deli = 1/(4*a_chief);
delLitOmg = 0;
delBigOmg = 0;
delM = 0;

delq1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
delq2 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
deltheta = delLitOmg + delM; % Relative true latitude in rad
delCOE = [dela, deltheta, deli, delq1, delq2, delBigOmg]';


% Integrating the deputies' initial trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Period = 2*pi*sqrt(a_chief^3/mu_Earth);
IntTime = 1*Period;
tspan  = linspace(0,IntTime,10000);
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-30);
mu_dot = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
OE_chief = [a_chief, theta_chief, inc_chief, q1_chief, q2_chief, BigOmg_chief];
AMap = ForwardMapping(OE_chief, mu_Earth); % Linear mapping matrix
Xo = AMap*delCOE;
u = zeros(3,1);
[time, X_nom] = ode113(@(t,X)CW_ode(t,X,u,mu_dot),tspan,Xo,options);


%%

Xo = nan(mm,Num_Agent);
for i = 1:Num_Agent
    idx = find(time <= i*IntTime/Num_Agent, 1,'last');
    Xo(:,i) = X_nom(idx,:)';
end
% [~, X_rel] = ode113(@(t,x)CW_ode(t,x,u,mu_dot),tspan,Xo,options);


% Specify the final relative-orbital elements of the deputies in both orbits
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dela = 0;
dele1 = 1/(6*a_chief);
deli1 = 1/(2*a_chief);
delLitOmg = 0;
delBigOmg = 0;
delM = 0;
delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
deltheta = delLitOmg + delM; % Relative true latitude in rad
delCOE_1 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';

delBigOmg = -pi*1E-5;
delM = pi*1E-6;
delLitOmg = pi*1E-5;
deltheta = delLitOmg + delM; % Relative true latitude in rad
dele2 = 1/(8*a_chief);
deli2 = 2/(10*a_chief);
delq1_2 = dele2*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
delq2_2 = dele2*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
delCOE_2 = [dela, deltheta, deli2, delq1_2, delq2_2, delBigOmg]';

% Integrating the deputy final trajectories
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Xo_1 = AMap*delCOE_1;
% u = zeros(3,1);
% [time, X_nom_1] = ode113(@(t,X)CW_ode(t,X,u,mu_dot),tspan,Xo_1,options);
% 
% Xo_1 = nan(mm,Num_Agent1);
% for i = 1:Num_Agent1
%     idx = find(time <= i*IntTime/Num_Agent1, 1,'last');
%     Xo_1(:,i) = X_nom_1(idx,:)';
% end
% [~, X_rel_1] = ode113(@(t,x)CW_ode(t,x,u,mu_dot),tspan,Xo_1,options);

Xo_2 = AMap*delCOE_2;
% [time, X_nom_2] = ode113(@(t,X)CW_ode(t,X,u,mu_dot),tspan,Xo_2,options);

% Xo_2 = nan(mm,Num_Agent2);
% for i = 1:Num_Agent2
%     idx = find(time <= (2*i-1)*IntTime/Num_Agent2, 1,'last');
%     Xo_2(:,i) = X_nom_2(idx,:)';
% end
% [~, X_rel_2] = ode113(@(t,x)CW_ode(t,x,u,mu_dot),tspan,Xo_2,options);


%% Geometric Motion Planning
%===========================%
%---------------------------- AGHF parameters ----------------------------
% Boundary conditions
clc
X0 = Xo(:)/Period;
Xf = [Xo_1 Xo_2]; 
Xf = Xf(:)/Period;
% smax for steady state
smax = 9e7; % Need tpo find an analytic way to solve for this value
% # of grids in t
tgrids = 500;
% # of grids in s
sgrids = 500;
% # of grids of integration for final trajectory extraction
intgrids = 1e5;
% # of states
N = length(X0);
% # of inputs
M = 3;
% motion duration
T = .5;
% penalty value (the lambda in paper)
Penalty = 10; 

%%
PysicalTime = linspace(0,T,tgrids); 
% global variable for the Euler Lagrange function and Metric G
global get_EL get_G 

% generate functions for EL, metric G, drift term f, control matrix F
% please modify this function if system is changed
[get_EL, get_G, get_f, get_F] = generate_fh_of_model(N,M,Num_Agent); 

% opts = odeset('RelTol',2.22045e-9,'AbsTol',2.22045e-20);
int = 1;
opts = odeset('AbsTol',1e-14);

%% ------------------------- Solving the AGHF -----------------------------
% solve for trajectory, see "AGHF" function file for implementation details
tmax = smax; xpoints = tgrids; xpoints_new = intgrids; tpoints = sgrids; 
k = Penalty; X = X0; AgentNum = Num_Agent;
m = 0;
x = linspace(0,T,xpoints);                  % discretization of the curve in time
t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of s interval in log scale
tic;
disp('solving PDE...');
% Solve AGHF
sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,k,N),...
              @(x) mypdexic(x, X, Xf, T),...
              @(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t, X, Xf, N),...
              x,t);
% The solution "sol" is of form sol(t,x,i)
toc;

%% ------------------ Calculating the action integral ---------------------
tic;
disp('Computing the action integral...');
% calculate the atuated curve length for each grid of s
X_temp = zeros(N,xpoints);
dX_temp = zeros(N,xpoints);
cost(tpoints) = 0;
for j = 1:tpoints
    for kk = 1:N
    [X_temp(kk,:), dX_temp(kk,:)]=pdeval(m,x,sol(j,:,kk),x);
    end
    %for each s grid, integrate  (x_dot-f(x))'*G(x_dot-f(x)) to get curve
    %length
    for i = 1:xpoints
        X = X_temp(:,i);
        dX = dX_temp(:,i);
        f = get_f(X);
        G = get_G(X,k);
        cost(j) = cost(j) + (dX-f)'*G*(dX-f)*(T/(xpoints-1));
    end
end
toc;


tpointsCorase = 400;
xpointsCorase = 400;
Tmesh = [0 logspace(-4,log10(tmax),tpointsCorase-1)]; % discretization of time
Xmesh = linspace(0,T,xpointsCorase); % discretization of the curve
[X,Y] = meshgrid(Tmesh,Xmesh);

Agent1_u1 = interp2(t,x,Period.*sol(:,:,1),X,Y);
Agent1_u2 = interp2(t,x,Period.*sol(:,:,2),X,Y);   
Agent1_u3 = interp2(t,x,Period.*sol(:,:,3),X,Y);


%% ----------------------- control extraction ----------------------------- 
disp('Extracting the required control...');
% initialize controls and states
Crtl = struct('u', cell(1, AgentNum)); % actual control
xnew = linspace(0,T,xpoints_new); % time grids for integration

% use spline to interpolate state from steady state pde solution, p
p         = zeros(N,xpoints_new);
dpdx      = zeros(N,xpoints_new);
dx_new    = diff([xnew(3),xnew,xnew(xpoints_new-2)]);
for i = 1:N
    p(i,:)         = spline(x,sol(end,:,i),xnew);
    % Numerically computing the time derivative of the state
    dp_new    = diff([p(i,3),p(i,:),p(i,xpoints_new-2)]);
    dpdx(i,:) = (dp_new(1:end-1)./dx_new(1:end-1).*dx_new(2:end) ...
        + dp_new(2:end)./dx_new(2:end).*dx_new(1:end-1)) ...
        ./ (dx_new(2:end)+dx_new(1:end-1));
end
% 
mm = 6;
for i = 1 : xpoints_new
    % f(x(t)), note: the F_d(x) in paper
    drift = get_f(p(:,i));
    % get [Fc F], note: the F_bar matrix in paper
    Ffull = get_F(p(:,i));
    B = Ffull(:,M+1:end);
    
    % Extracting/approximate the required control
    for jj = 1:AgentNum
        idx = 1+mm*(jj-1):mm*jj;
        Crtl(jj).u(:,i) = (B.'*B)\(B.')* ( dpdx(idx,i) - drift(idx) );
    end
end
disp('Done!!!!!');


%% plotting the results
%===========================%
%------------------- Plotting the extracted control -----------------------
Int_time = 0:1:intgrids-1;
figure('Name','Extracted control','Renderer', 'painters', 'Position', [10 10 900 600])
Label = {'Agent 1', 'Agent 2', 'Agent 3', 'Agent 4'};
for i = 1:Num_Agent
    subplot(2,2,i)
    plot(Int_time,Crtl(i).u(1,:),'r','LineWidth',2);
    hold on 
    plot(Int_time,Crtl(i).u(2,:),'b','LineWidth',2);
    plot(Int_time,Crtl(i).u(3,:),'g','LineWidth',2);
    set(gca, 'XScale', 'log')
    xlabel('time (s)')
    ylabel(Label(i))
    grid on;
    
end
%%
figure('Name','Action integral across iteration','Renderer', 'painters', 'Position', [10 10 900 600])
loglog(t,cost);
title('Action integral across iteration')
xlabel('Homotopy iteration')
ylabel('curve length')
grid on;


%%
for i = 1
    % close all
    c1 = rgb('DarkGreen');
    c2 = rgb('Gold');
    c3 = rgb('Lime');
    c4 = rgb('Orange');
    c5 = rgb('DarkBlue');
    c6 = rgb('Red');
    c7 = rgb('DarkGray');
    c8 = rgb('Bisque');
    c9 = rgb('Teal');
    
    cMat1 = [c6;c5;c2;c3];
    cMat2 = [c1;c4];
    cMat = [];
    Label_1 = {'Agent1','Agent2','Agent3','Agent4'};
    Label_2 = {'Agent5','Agent6'};
    Label = {'Agent1','Agent2','Agent3','Agent4','Agent5','Agent6'};
    
    fig = figure('Name','3D view of the motion','Renderer', 'painters', 'Position', [10 10 800 600]);
    
    hold on
    for jj=1%:Num_Agent1
        idx1 = 1+mm*(jj-1):mm*jj;
        
        plot3(X_rel_1(:,idx1(1)),X_rel_1(:,idx1(2)),X_rel_1(:,idx1(3)),'-.','Color',c8);
        plot3(X_rel_1(1,idx1(1)),X_rel_1(1,idx1(2)),X_rel_1(1,idx1(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat1(jj,:),...
            'MarkerFaceColor',cMat1(jj,:)',...
            'MarkerSize',8);
        hold on
        
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
        cMat = [cMat; cMat1(jj,:)];
    end
    
    plt   = zeros(Num_Agent);
    hold on
    for jj = 1%:Num_Agent
        idx = 1+mm*(jj-1):mm*jj;
        txt = ['Agent_',num2str(jj)];
        plt(:,jj) = plot3(X_rel(1,idx(1)),X_rel(1,idx(2)),X_rel(1,idx(3)),'*',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat(jj,:),...
            'MarkerFaceColor',cMat(jj,:)',...
            'MarkerSize',8);
        
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
    end
    grid on
    h5 = plot3(0,0,0,'bo',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',15);
    h4 = plot3(X_rel(:,idx(1)),X_rel(:,idx(2)),X_rel(:,idx(3)),'--','Color',c7,'DisplayName','Initial Traj');
    h6 = plot3(X_rel_1(:,idx1(1)),X_rel_1(:,idx1(2)),X_rel_1(:,idx1(3)),'--','Color',c7,'DisplayName','Final Traj');
    view(21,23)
    axis equal
end

r = 1;
[x_sphere, y_sphere, z_sphere] = sphere(100);
h = surf(r*x_sphere, r*y_sphere, r*z_sphere);
set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
axis equal

% figure('Name','3D configuration space curve deformation','Renderer', 'painters', 'Position', [10 10 900 600])
X = Agent1_u1(1,:);
Y = Agent1_u2(1,:);
Z = Agent1_u3(1,:);
h1 = plot3(X,Y,Z,'r','LineWidth',2);
hold on
h2 = plot3(X,Y,Z,'Color',c9,'LineWidth',2);
grid ON;
title('Nonlinear Geomeric Planning - 3D configuration');
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
hL = legend([plt(end,1), h5, h2, h1],...
    {'Agent1','Chief','Initial Guess','Homotopy Iteration'},'AutoUpdate','off');
% Programatically move the Legend
newPosition = [.23 .65 0.01 0.01];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);


set(gca,'nextplot','replacechildren','visible','off')
f  = getframe;
[im,map] = rgb2ind(f.cdata,1000,'nodither'); % 65536
im(1,1,1,tpointsCorase) = 0;

for i = 1:tpointsCorase
    
    X=Agent1_u1(i,:);h1.XDataSource='X';
    Y=Agent1_u2(i,:);h1.YDataSource='Y';
    Z=Agent1_u3(i,:);h1.ZDataSource='Z';
    refreshdata(h1,'caller');
    drawnow;
    
    f = getframe;
    im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,'AGHF_CWEquation.gif','DelayTime',0,'LoopCount',inf) %g443800


%% Integrate the system's trajectory 
clear X_ode45
tolerance = 1e-16;
h = 1e-9;
X_ode45(:,1) = Xo(:);
xnew = linspace(0,T*Period,intgrids); % time grids for integration
for i = 1:intgrids-1
    for jj = 1:AgentNum
        u(:,jj) = Crtl(jj).u(:,i);
    end
    [y_f, ~, h_next] = Runge_Kutta_Fehlberg_4_5(@(t,X)CW_ode(t,X,u),X_ode45(:,i),xnew(i),h,xnew(i+1),tolerance);
    h = h_next;
    X_ode45(:,i+1) = y_f;
end


%% animation
for Penalty = 1
    % close all
    c1 = rgb('DarkGreen');
    c2 = rgb('Gold');
    c3 = rgb('Lime');
    c4 = rgb('Orange');
    c5 = rgb('DarkBlue');
    c6 = rgb('Red');
    c7 = rgb('DarkGray');
    c8 = rgb('Bisque');
    
    cMat1 = [c6;c5;c2;c3];
    cMat2 = [c1;c4];
    cMat = [];
    Label_1 = {'Agent1','Agent2','Agent3','Agent4'};
    Label_2 = {'Agent5','Agent6'};
    Label = {'Agent1','Agent2','Agent3','Agent4','Agent5','Agent6'};
    
    fig = figure('Name','3D view of the motion','Renderer', 'painters', 'Position', [10 10 700 550]);
    hold on
    for jj=1:Num_Agent1
        idx1 = 1+mm*(jj-1):mm*jj;
        
        plot3(X_rel_1(:,idx1(1)),X_rel_1(:,idx1(2)),X_rel_1(:,idx1(3)),'-.','Color',c8);
        plot3(X_rel_1(1,idx1(1)),X_rel_1(1,idx1(2)),X_rel_1(1,idx1(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat1(jj,:),...
            'MarkerFaceColor',cMat1(jj,:)',...
            'MarkerSize',8);
        hold on
        
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
        cMat = [cMat; cMat1(jj,:)];
    end
    
    
    for jj=1:Num_Agent2
        idx2 = 1+mm*(jj-1):mm*jj;
        plot3(X_rel_2(:,idx2(1)),X_rel_2(:,idx2(2)),X_rel_2(:,idx2(3)),'-.','Color',c8);
        plot3(X_rel_2(1,idx2(1)),X_rel_2(1,idx2(2)),X_rel_2(1,idx2(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat2(jj,:),...
            'MarkerFaceColor',cMat2(jj,:)',...
            'MarkerSize',8);
        cMat = [cMat; cMat2(jj,:)];
    end
    
    plt   = zeros(Num_Agent);
    hold on
    for jj = 1:Num_Agent
        idx = 1+mm*(jj-1):mm*jj;
        txt = ['Agent_',num2str(jj)];
        plt(:,jj) = plot3(X_rel(1,idx(1)),X_rel(1,idx(2)),X_rel(1,idx(3)),'*',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat(jj,:),...
            'MarkerFaceColor',cMat(jj,:)',...
            'MarkerSize',8);
        
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
    end
    grid on
    h5 = plot3(0,0,0,'bo',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',15);
    h4 = plot3(X_rel(:,idx(1)),X_rel(:,idx(2)),X_rel(:,idx(3)),'--','Color',c7,'DisplayName','Initial Traj');
    h6 = plot3(X_rel_1(:,idx1(1)),X_rel_1(:,idx1(2)),X_rel_1(:,idx1(3)),'--','Color',c7,'DisplayName','Final Traj');
    
    r = 1;
    [x_sphere y_sphere z_sphere] = sphere(100);
    h = surf(r*x_sphere, r*y_sphere, r*z_sphere);
    set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
    axis equal
    hL = legend([plt(end,1), plt(end,2), plt(end,3), plt(end,4), h5, h],...
        {'Agent1','Agent2','Agent3','Agent4','Chief','Constraints'},'AutoUpdate','off');
    % Programatically move the Legend
    newPosition = [.19 .68 0.01 0.01];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    view(31,18)
    
    for k = 1:Num_Agent
        idx = 1+mm*(k-1):mm*k;
        txt = ['Traj Agent_',num2str(k)];
        plot3(interp1(xnew,X_ode45(idx(1),:),t),interp1(xnew,X_ode45(idx(2),:),t),interp1(xnew,X_ode45(idx(3),:),t),...
            '--','Color',cMat(k,:),'LineWidth',2,'DisplayName',txt);
    end
    
end


%% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

% The followings are the PDE, initial condition and boundary condition for pdepe:

function [c,f,s] = mypdexpde(x,t,u,DuDx,k,N)    % Define PDE; right-hand-side of AGHF

global get_EL

% evaluate EL

EL = get_EL(u,DuDx,k);
% dL/dx
pLx = EL(:,1);
% dL/d(dot_x)
pLxd = EL(:,2);


f = pLxd;
s = -pLx;

c = ones(N,1);

% for details of f, s, c, please see the documentation for pdepe

end
% -------------------------------------------------------------------------
 
function u0 = mypdexic(x, X0, Xf, T)    % Initial condition for AGHF

% global iteration InitalCondition PysicalTime
% idx = find(PysicalTime <= x, 1, 'last' );
% if iteration == 1

%%  straight line connecting X0 and Xf
% u0=X0+(Xf-X0)*(x/T);
% % add some sinusoidal deviation to x1 initial guess (can be discarded)
% u0(1)=X0(1)+(Xf(1)-X0(1))*(x/T) + 0.01*sin(2*pi*x/T);

%% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;

u0 = X0*cos(freq1*x/T) ...
    + Xf*sin(freq1*x/T) ...
    - (X0+Xf)*sin(freq2*x/T);

end
% -------------------------------------------------------------------------

function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t, X0, Xf, N)    % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end
% -------------------------------------------------------------------------

% generate functions for calculating EL equation:
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!! Need to modify this function if system is changed !!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Inputs:
% N,M -- # of states and inputs
% Output:
% get_EL -- function handle for calculationg EL equation, it returns dL/dx
% and % dL/d(dot_x), it takes x and x_dot as input
% get_G -- function handle for calculating the metric G(x)
% get_f -- function handle for calculating the drift term f(x)
% get_F -- function handle for calculating the full control matrix [Fc(x) F(x)]

function [get_EL, get_G, get_f, get_F] = generate_fh_of_model(N,M,AgentNum) 
tic;
disp('generating functions for pdepe...');

% ------------------------- system setup ----------------------------------
syms Penalty
global mu_dot Period
A  = [zeros(3) eye(3);
      3*mu_dot^2 0 0 0 2*mu_dot 0;
      0 0 0 -2*mu_dot 0 0;
      0 0 -mu_dot^2 0 0 0]*Period; % Partial derivative pFd the drift vector field
N = 6;
Xaug = []; dXaug = [];
  % make the state and state deirvative simbolic
  for i = 1:AgentNum
      X_Variables{i} = [{strcat('x',num2str(i))},{strcat('y',num2str(i))},...
          {strcat('z',num2str(i))},{strcat('xdot',num2str(i))},...
          {strcat('ydot',num2str(i))},{strcat('zdot',num2str(i))}];
      dX_Variables{i} = [{strcat('dx',num2str(i))},{strcat('dy',num2str(i))},...
          {strcat('dz',num2str(i))},{strcat('dxdot',num2str(i))},...
          {strcat('dydot',num2str(i))},{strcat('dzdot',num2str(i))}];
      
      if i == 1
          AMat_Aug = A;
          % [Fc F], note: the F_bar matrix in paper
          F = sym(eye(N)); F_Aug = F;
          % penalty matrix (D matrix in the paper)
          D = diag([Penalty*ones(1,N-M) ones(1,M)]); D_Aug = D;
      else
          AMat_Aug = blkdiag(AMat_Aug,A);
          F_Aug = blkdiag(F_Aug,F);
          D_Aug = blkdiag(D_Aug,D);
      end
  end
% syms x y z xdot ydot zdot  dx dy dz dxdot dydot dzdot Penalty real
X1  = cell2sym( X_Variables{1}' ); assume(X1,'real')
X2  = cell2sym( X_Variables{2}' ); assume(X2,'real')
% X3  = cell2sym( X_Variables{3}' ); assume(X3,'real')
% X4  = cell2sym( X_Variables{4}' ); assume(X4,'real')
% X5  = cell2sym( X_Variables{5} );
% X6  = cell2sym( X_Variables{6} );

dX1 = cell2sym( dX_Variables{1}' ); assume(dX1,'real')
dX2 = cell2sym( dX_Variables{2}' ); assume(dX2,'real')
% dX3 = cell2sym( dX_Variables{3}' ); assume(dX3,'real')
% dX4 = cell2sym( dX_Variables{4}' ); assume(dX4,'real')
% dX5 = cell2sym( dX_Variables{5} );
% dX6 = cell2sym( dX_Variables{6} );

Xaug  = [X1(:);X2(:)]; %;X3(:);X4(:)]; % ;X5(:);X6(:)];
dXaug = [dX1(:);dX2(:)]; %;dX3(:);dX4(:)]; % ;dX5(:);dX6(:)];


% drift term f(x) note: Fd(x) is used in the paper
f = AMat_Aug*Xaug;
% the metric, before adding the state constraint barrier functions
H = simplify( (F_Aug.')^(-1)*D_Aug*F_Aug^(-1) );


% -------------------- state constraint barrier function ------------------
% B is a vector of barrier function for 1 state constraint, each scaler state
% constraint need 1 barrier function
B = [];
% barrier function parameters
% kb = .01; % barrier function gain
% pb = 1;   % order of barrier function

% b is one penalty term for state constraint
% b = 0; % no state contraintst
% delt12 = X1(1:3) - X2(1:3); b12 = (delt12'*delt12 - (0.05)^2)/(delt12'*delt12)^pb;
% delt13 = X1(1:3) - X3(1:3); b13 = (delt13'*delt13 - (0.05)^2)/(delt13'*delt13)^pb;
% delt14 = X1(1:3) - X4(1:3); b14 = (delt14'*delt14 - (0.05)^2)/(delt14'*delt14)^pb;
% delt23 = X2(1:3) - X3(1:3); b23 = (delt23'*delt23 - (0.05)^2)/(delt23'*delt23)^pb;
% delt24 = X2(1:3) - X4(1:3); b24 = (delt24'*delt24 - (0.05)^2)/(delt24'*delt24)^pb;
% delt34 = X3(1:3) - X4(1:3); b34 = (delt34'*delt34 - (0.05)^2)/(delt34'*delt34)^pb;
% b = b12+b13+b14+b23+b24+b34;
% b21 = kb/((0.5)^2 - norm(X1(1:3))^2)^pb;
% b22 = kb/((0.5)^2 - norm(X2(1:3))^2)^pb;
% b23 = kb/((0.5)^2 - norm(X3(1:3))^2)^pb;
% b24 = kb/((0.5)^2 - norm(X4(1:3))^2)^pb;
% b = b21 + b22 + b23 + b24;
% b = kb/((2*pi/4)^2 - z2^2)^pb; %contraint on angular velocity, |z2|<pi/2
% b = kb/(2^2 - z1^2)^pb; %contraint on translational velocity, |z1|<2
% b = kb/(z1+2)^pb + kb/(2-z1)^pb; %contraint on translational velocity |z1|<2
b =0;
B = [B b]; % attach the barrier function to B

% -------------------- Metric and curve Length ----------------------------

% the metric with state constraint barier functions
G = (sum(B)+1)*H;
% actuated curve length
L = (dXaug - f).' * G * (dXaug - f) ;

% -------------------- Function Generations -------------------------------

% taking derivatives symbolically to get the EL euqation terms
% pLx -- dL/dx
% pLxd -- dL/d(x_dot)
pLx  = sym('pLx', [N*AgentNum 1],'real');
pLxd = sym('pLxd', [N*AgentNum 1],'real');
for i=1:N*AgentNum
    pLx(i)  = diff(L,Xaug(i));
    pLxd(i) = diff(L,dXaug(i));
end

% generate functions
get_EL = matlabFunction([pLx, pLxd],'vars', {Xaug,dXaug,Penalty});
get_G  = matlabFunction(G,'vars', {Xaug,Penalty});
get_f  = matlabFunction(f,'vars',{Xaug});
get_F  = matlabFunction(F,'vars',{Xaug});

toc;
end


function [xdot] = MyOde(t,x,u)
format long 
global mu_dot Period
% Reshaping the input into an 6 by n matrix
[l,~] = size(x);
x = reshape(x,6,l/6);

% Computing the CW control and dynamic matrices
B = [zeros(3); eye(3)];
A = [zeros(3) eye(3);
    3*mu_dot^2 0 0 0 2*mu_dot 0;
    0 0 0 -2*mu_dot 0 0;
    0 0 -mu_dot^2 0 0 0]*Period;

% Computing the derivatives for the ode
xdot = A*x + B*u;
xdot = xdot(:);

end
% -------------------------------------------------------------------------
function [Grad_x,Grad_dx] = AutoDiff_Dynamics(t,x,dx,K)
global M Num_Agent
  
N = 6;
for i = 1:Num_Agent
    if i == 1
        % [Fc F], note: the F_bar matrix in paper
        F = eye(N); F_Aug = F;
        % penalty matrix (D matrix in the paper)
        D = diag([K*ones(1,N-M) ones(1,M)]); D_Aug = D;
    else
        F_Aug = blkdiag(F_Aug,F);
        D_Aug = blkdiag(D_Aug,D);
    end
end

%% Defining the drift vector field
for i = 1:Num_Agent
    idx = 1+N*(i-1):N*i;
    X0 = x(idx);
    f(idx,1) = NonDimentionalized_CWode(t, X0); % NonDimentionalized_NLode(t, X0);
end


%% Defining the metric, before adding the state constraint barrier functions
H = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);

%% ------------------- state constraint barrier function ------------------
B = [];
b = 0; % no state contraintst
B = [B b]; % attach the barrier function to B

% -------------------- Metric and curve Length ----------------------------
% the metric with state constraint barier functions
G = (sum(B)+1)*H;
% actuated curve length
L = (dx - f).' * G * (dx - f);

Grad_x  = extractdata(dlgradient(L,x));
Grad_dx = extractdata(dlgradient(L,dx));

end

function [Xdot] = NonDimentionalized_CWode(t,X)
global mu_dot Tc

% Computing the CW control and dynamic matrices
A  = [zeros(3) eye(3);
      3*mu_dot^2 0 0 0 2*mu_dot 0;
      0 0 0 -2*mu_dot 0 0;
      0 0 -mu_dot^2 0 0 0]*Tc;

% Reshaping the input into an 6 by n matrix
[l,~] = size(X);
X = reshape(X,6,l/6);
Xdot = A*X;
Xdot = Xdot(:);
end






