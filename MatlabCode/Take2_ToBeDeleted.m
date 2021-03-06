%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Hermann Kaptui Sipowa %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear
close all
start_up
format long e
clc
global mu_Earth tspan M T_normalized mu_dot ...
    Num_Agent a_chief e_chief DynamicsModel ...
    Period Rrho Tc  % iteration InitalCondition

c1  = rgb('DarkGreen');
c2  = rgb('DarkSlateGray');
c3  = rgb('Lime');
c4  = rgb('DarkOrange');
c5  = rgb('DarkBlue');
c6  = rgb('Red');
c7  = rgb('Purple');
c8  = rgb('Bisque');
c9  = rgb('Orange');
c10 = rgb('DarkGray');
c11 = rgb('Teal');
mm = 6;
addpath('/Users/hermann/Desktop/casadi-osx-matlabR2015a-v3.5.5')
import casadi.*
CasADiopts = struct('main', true, 'mex', true);

%%
%------------------- Specify the dynamics Model ------------------------  
% Chose from: 
% - CWode = Clohessy–Wiltshire, 
% - THode = Tschauner-Hempel , 
% - NLode = Full Nonlinear Model 
DynamicsModel = 'THode';
SimTime = 0.75;

%% Calculating the boundaries conditions
%=======================================%
mu_Earth   = 3.986004415E5; % Earth gravitational parameter
Num_Agent  = 2; % Number of agents in the system
Num_Agent1 = 1; % Number of agents in the first final-orbit
Num_Agent2 = 1; % Number of agents in the second final-orbit

% Specifying the chief's orbit
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a_chief      = 1.42e4;   % Semi-major axis in Km
e_chief      = 0.5;      % Eccentricity
inc_chief    = 50;       % Inclination in deg0
BigOmg_chief = 10;       % RAAN in deg
LitOmg_chief = 10;       % AOP in deg
M_chief1     = 0;        % Initial mean anomaly
OE_Chief     = [a_chief,e_chief,inc_chief,...
                BigOmg_chief,LitOmg_chief,M_chief1];


%% Integrating the chief's trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Period  = (2*pi)*sqrt(a_chief^3/mu_Earth);
IntTime = SimTime*Period;
tspan   = linspace(0,IntTime,1e4);
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-17);
mu_dot  = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
Rc  = a_chief*eye(3); Rrho = eye(3); Tc = Period; T_normalized = tspan/Tc;
inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
LitOmg_chief = deg2rad(LitOmg_chief); M_chief1 = deg2rad(M_chief1);
COE1 = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief1];
[Position_target1,Velocity_target1] = COEstoRV(COE1,mu_Earth);
X0_Chief1  = [Position_target1; Velocity_target1];

n = sqrt(mu_Earth/a_chief^3); M_chief2 = M_chief1 + n*IntTime;
COE2 = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief2];
[Position_target2,Velocity_target2] = COEstoRV(COE2,mu_Earth);
X0_Chief2  = [Position_target2; Velocity_target2];
Ttraj = linspace(0,Period,1e4)/Tc; 

%% --------------- Setting the deputy initial conditions ------------------
for i = 1
    r0 = X0_Chief1(1:3); v0 = X0_Chief1(4:6);
    COE0 = RVtoCOEs(r0,v0,mu_Earth); f1 = COE0(6);
    q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
    theta_chief1 = f1+LitOmg_chief;
    OE_chief1 = [a_chief, theta_chief1, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    AMap1 = ForwardMapping(OE_chief1, mu_Earth); % Linear mapping matrix
    
    rf = X0_Chief2(1:3); vf = X0_Chief2(4:6);
    COEf = RVtoCOEs(rf,vf,mu_Earth); f2 = COEf(6);
    theta_chief2 = f2+LitOmg_chief;
    OE_chief2 = [a_chief, theta_chief2, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    AMap2 = ForwardMapping(OE_chief2, mu_Earth); % Linear mapping matrix
    
    % Specify the final relative-orbital elements of the deputies in both orbits
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    dela = 0;
    dele1 = 1/(6*a_chief);
    deli1 = -1/(3*a_chief);
    delLitOmg = -2*pi*1E-5;
    delBigOmg = 0;
    delM = pi*1E-6;
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_1 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    
    dele1 = 1/(8*a_chief);
    deli1 = 2/(10*a_chief);
    delLitOmg = pi*1E-5;
    delBigOmg = -pi*1E-5;
    delM = pi*1E-6;
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_2 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    
    dele1 = 1/(6*a_chief);
    deli1 = 1/(2*a_chief);
    delLitOmg = 0;
    delBigOmg = 0;
    delM = 0;
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_3 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    
    dele1 = -1/(8*a_chief);
    deli1 = -1/(20*a_chief);
    delLitOmg = -5*pi*1E-7;
    delBigOmg = pi*1E-6;
    delM = 0;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    delCOE_4 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    % Integrating the deputy initial and final trajectories
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    X01 = AMap1*delCOE_1; Xi1 = X01;
    X02 = AMap1*delCOE_2; Xi2 = X02;
    Xf1 = AMap2*delCOE_3;
    Xf2 = AMap2*delCOE_4;
    
end 
% -------------------------------------------------------------------------

%% ------------------------ AutoDiff using CasADi -------------------------
for j = 1
    XChief= SX.sym('XC',6); t = SX.sym('t'); Penalty = SX.sym('Penalty');
    rho1  = SX.sym('rho1',6); rho2    = SX.sym('rho2',6);
    drho1 = SX.sym('rho1',6); drho2   = SX.sym('rho2',6);
    Xaug  = [rho1;rho2];      dXaug   = [drho1;drho2];
    M = 3; N = 6; drift = SX(N*Num_Agent,1);
    
    if exist('L','file') == 3
        delete L.c L.mexmaci64
    end
    if exist('dLdx','file') == 3
        delete dLdx.c dLdx.mexmaci64
    end
    if exist('fdrift','file') == 3
        delete fdrift.c fdrift.mexmaci64
    end
    
    for i = 1:Num_Agent
        if i == 1
            % [Fc F], note: the F_bar matrix in paper
            F = eye(N); F_Aug = F;
            % penalty matrix (D matrix in the paper)
            D = diag([Penalty*ones(1,N-M) ones(1,M)]); D_Aug = D;
        else
            F_Aug = blkdiag(F_Aug,F);
            D_Aug = blkdiag(D_Aug,D);
        end
    end
    
    % Defining the drift vector field
    for i = 1:Num_Agent
        idx = 1+N*(i-1):N*i;
        if strcmp(DynamicsModel,'THode')
            drift(idx,1) = NonDimentionalized_THodeCasADi(t,Xaug(idx),XChief);
        elseif strcmp(DynamicsModel,'CWode')
            drift(idx,1) = NonDimentionalized_CWodeCasADi(t,Xaug(idx),XChief);
        elseif strcmp(DynamicsModel,'NLode')
            drift(idx,1) = NonDimentionalized_NLodeCasADi(t,Xaug(idx),XChief);
        end
    end
    
    % Defining the metric, before adding the state constraint barrier functions
    H = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
    
    % ------------------- state constraint barrier function ------------------
    Threshold = 0.45;
    % kb = 1e-7; % barrier function gain
    % pb = 1; % order of barrier function
    % delR = norm(rho1(1:3) - rho2(1:3));
    % B = kb/(delR - Threshold)^pb;
    B = [];
    
    
    % -------------------- Metric and curve Length ----------------------------
    % the metric with state constraint barier functions
    G = (sum(B)+1)*H;
    % actuated curve length
    L = Function('L',{Xaug,dXaug,Penalty,t,XChief},...
        {(dXaug - drift).' * G * (dXaug - drift)},...
        {'X','dX','k','t','XChief'},...
        {'CurveLength'});
    dLdx = L.jacobian();
    fdrift = Function('fdrift',{t,Xaug,XChief},{drift});
    
    % Rewrite the functions in c files
    L.generate('L.c',CasADiopts);
    dLdx.generate('dLdx.c',CasADiopts);
    fdrift.generate('fdrift.c',CasADiopts);

    % Generating Mex files from the c files
    mex L.c -largeArrayDims
    mex dLdx.c -largeArrayDims
    mex fdrift.c -largeArrayDims
    
    clear L dLdx fdrift 
end
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% ##################### Geometric Motion Planning #########################
% -------------------------------------------------------------------------
%% ------- Nondimentionalized the inital condition of the problem ---------
for j = 1
    
    h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
    K_f = 1+e_chief*cos(f1);
    conts1 = mu_Earth*K_f/(h^(3/2));
    conts2 = -mu_Earth*e_chief*sin(f1)/(h^(3/2));
    conts3 = 1/conts1;
    A_map1 = [conts1*eye(3)       zeros(3);
        conts2*eye(3)  conts3*eye(3)];
    x_bar01 = A_map1*X01;
    x_bar02 = A_map1*X02;
    X0  = [x_bar01; x_bar02];
    tspan1 = 2*pi*Ttraj;
    h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
    K_f    = 1+e_chief*cos(f2);
    conts1 = mu_Earth*K_f/(h^(3/2));
    conts2 = -mu_Earth*e_chief*sin(f2)/(h^(3/2));
    conts3 = 1/conts1;
    A_map2 = [conts1*eye(3)       zeros(3);
        conts2*eye(3)  conts3*eye(3)];
    x_bar11 = A_map2*Xf1;
    x_bar21 = A_map2*Xf2;
    Xf = [x_bar11; x_bar21];
    tspan2 = linspace(f2,(f2+2*pi),length(T_normalized));
    tic
    [~, X_rel_0] = ode113(@(t,x)NonDimentionalized_THode(t,x),tspan1,X0,options);
    [~, X_rel_1] = ode113(@(t,x)NonDimentionalized_THode(t,x),tspan2,Xf,options);
    toc
    X_THpos0 = nan(length(X_rel_0),6);
    X_THpos1 = nan(length(X_rel_0),6);
    X_THpos2 = nan(length(X_rel_0),6);
    X_THpos3 = nan(length(X_rel_0),6);
    for i = 1:length(X_rel_0)
        f = tspan1(i);
        h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
        K_f = 1+e_chief*cos(f);
        conts1 = h^(3/2)/(mu_Earth*K_f);
        conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
        conts3 = 1/conts1;
        A_inv1 = [conts1*eye(3)       zeros(3);
            conts2*eye(3)  conts3*eye(3)];
        X_THpos0(i,:) = (A_inv1*X_rel_0(i,1:6).').';
        X_THpos1(i,:) = (A_inv1*X_rel_0(i,7:12).').';
        
        f = tspan2(i);
        h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
        K_f = 1+e_chief*cos(f);
        conts1 = h^(3/2)/(mu_Earth*K_f);
        conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
        conts3 = 1/conts1;
        A_inv2 = [conts1*eye(3)       zeros(3);
            conts2*eye(3)  conts3*eye(3)];
        X_THpos2(i,:) = (A_inv2*X_rel_1(i,1:6).').';
        X_THpos3(i,:) = (A_inv2*X_rel_1(i,7:12).').';
    end
    X_THpos01 = [X_THpos0 X_THpos1];
    X_THpos11 = [X_THpos2 X_THpos3];
    
    
    X0 = X0(:); Xf = Xf(:);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    for k = 1
        
        Position_target2 = X0_Chief2(1:3); Velocity_target2 = X0_Chief2(4:6);
        
        % Converting angles in radians
        LVLH1 = DCM(Position_target1,Velocity_target1);
        N_nudot1 = cross(Position_target1,Velocity_target1)/norm(Position_target1)^2;
        Position_chaser = Position_target1 + LVLH1\X01(1:3);
        Velocity_chaser = Velocity_target1 + cross(N_nudot1,(LVLH1\X01(1:3))) + LVLH1\X01(4:6);
        
        LVLH2 = DCM(Position_target2,Velocity_target2);
        N_nudot2 = cross(Position_target2,Velocity_target2)/norm(Position_target2)^2;
        Position_Desired = Position_target2 + LVLH2\Xf1(1:3);
        Velocity_Desired = Velocity_target2 + cross(N_nudot2,(LVLH2\Xf1(1:3))) + LVLH2\Xf1(4:6);
        X0_desired = [Position_Desired;Velocity_Desired]; tspan_inv = flip(tspan);
        [~, Xdesired] = ode113(@(t,x)M2BodyOde(t,x,mu_Earth),tspan_inv,X0_desired,options);
        Position_Desired = Xdesired(end,1:3)';
        Velocity_Desired = Xdesired(end,4:6)';
    end
    
end
% -------------------------------------------------------------------------

%% --------------------------- AGHF parameters ----------------------------
for ll = 1
    clc
    smax = 1e7; % Need tpo find an analytic way to solve for this value
    % # of grids in t
    tgrids = 250;
    % # of grids in s
    sgrids = 250;
    % # of grids of integration for final trajectory extraction
    intgrids = 1e5;
    % # of inputs
    M = 3;  N  = length(X0);
    % penalty value (the lambda in paper)
    Penalty = 5e9;
    % tolerance of the integrator
    % opts = odeset('AbsTol',1e-14);
    opts = odeset('RelTol',2.22045e-8,'AbsTol',2.22045e-20);
    
    % Setting the boundary condtions and the intial condition
    tmax = smax; xpoints = tgrids; xpoints_new = intgrids; tpoints = sgrids;
    m = 0; t = [0 logspace(-2,log10(tmax),tpoints-1)]; % discretization of s interval in log scale
    T    = SimTime*wrapTo2Pi(f2)-f1; % motion duration
    x    = linspace(f1,SimTime*wrapTo2Pi(f2),xpoints); % discretization of the curve in time
    xnew = linspace(f1,SimTime*wrapTo2Pi(f2),intgrids);
    
end
% -------------------------------------------------------------------------

%% ------------------------- Solving the AGHF -----------------------------
for ll = 1
    % solve for trajectory, see "AGHF" function file for implementation details
    close all
    clc
    tic;
    disp('solving PDE...');
    
    % Solve AGHF
    sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,Penalty,N),...
        @(x) mypdexic(x,X0,Xf,T),...
        @(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t,X0,Xf,N),...
        x,t);
    % The solution "sol" is of form sol(t,x,i)
    toc;
end
% -------------------------------------------------------------------------

%% LQR controller



%% Lyapunov Controller
% Integrate the trajectory in cartesian coordinates (Truth)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
for k = 1
    % Controller's gains
    K = eye(3)*3E-7;
    P = eye(3)*1E-3;
    Control_effort = nan(length(tspan),3);
    Aug_chasser_Xo = [Position_chaser; Velocity_chaser];
    Aug_target_Xo  = [Position_target1; Velocity_target1];
    Aug_Desired_Xo = [Position_Desired; Velocity_Desired];
    Aug_Xo = [Aug_chasser_Xo; Aug_target_Xo; Aug_Desired_Xo];
    tic
    [time, Aug_X] = ode113(@(t,Aug_X)ControlUnperturbed2Bodyfunc(t,Aug_X,mu_Earth,K,P),tspan,Aug_Xo,options);
    toc
    
    for j = 1:length(time)
        % Extracting relevant quantities
        Chasser_r = Aug_X(j,1:3)'; Chasser_v = Aug_X(j,4:6)';
        Target_r = Aug_X(j,7:9)'; Target_v = Aug_X(j,10:12)';
        Desired_r = Aug_X(j,13:15)'; Desired_v = Aug_X(j,16:18)';
        
        % Computing the globally asymptotically stabilizing controller
        f_desired = -mu_Earth/norm(Desired_r(1:3))^3 * Desired_r(1:3);
        f_deputy =  -mu_Earth/norm(Chasser_r(1:3))^3*Chasser_r(1:3);
        del_r = Chasser_r(1:3) - Desired_r(1:3);
        del_v = Chasser_v(1:3) - Desired_v(1:3);
        Control_effort(j,:) = f_desired - f_deputy - K*del_r - P*del_v;
        
        % Computing the necessary for trasfrom ECI to LVLH
        h_vec = cross(Target_r,Target_v); h_norm = norm(h_vec);
        N_nudot = h_vec/norm(Target_r)^2;
        TN = DCM(Target_r,Target_v); NOmega_target = cross(Target_r,Target_v)/norm(Target_r)^2;
        
        % Coverting the deputy ECI coordinates to the LVLH frame
        NR_rel = Chasser_r - Target_r; NV_rel = Chasser_v - Target_v;
        TR_rel(j,:) = TN*NR_rel; TV_rel(j,:) = TN*(NV_rel - cross(N_nudot,NR_rel));
        
        % Coverting the desired trajectory ECI coordinates to the LVLH frame
        NR_rel_desired = Desired_r - Target_r; NV_rel_desired = Desired_v - Target_v;
        TR_rel_desired(j,:) = TN*NR_rel_desired; TV_rel_desired(j,:) = TN*(NV_rel_desired - cross(N_nudot,NR_rel_desired));
    end
    
    R_chasser1 = [Aug_X(:,1) Aug_X(:,2) Aug_X(:,3)]; V_chasser1 = [Aug_X(:,4) Aug_X(:,5) Aug_X(:,6)];
    R_target1 = [Aug_X(:,7) Aug_X(:,8) Aug_X(:,9)]; V_target1 = [Aug_X(:,10) Aug_X(:,11) Aug_X(:,12)];
    
end

% Plotting the 3D realtive trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
for k = 1%:3
    close all
    % Plotting the 3D trajectories
    figure
    for jj=1
        h1 = plot3(TR_rel(:,1),TR_rel(:,2),TR_rel(:,3),'g');
        hold on
        h3 = plot3(TR_rel_desired(end,1),TR_rel_desired(end,2),TR_rel_desired(end,3),'ro',...
            'LineWidth',5,...
            'MarkerEdgeColor','r',...
            'MarkerFaceColor','r',...
            'MarkerSize',5);
        h4 = plot3(TR_rel_desired(:,1),TR_rel_desired(:,2),TR_rel_desired(:,3),'r');
        h2 = plot3(TR_rel(1,1),TR_rel(1,2),TR_rel(1,3),'go',...
            'LineWidth',2,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor','b',...
            'MarkerSize',5);
        h5 = plot3(0,0,0,'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',c2,...
            'MarkerSize',15);
        grid on
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
        legend([h5, h2, h3, h1, h4],{'Chief','Deputy IC','Deputy FC','Deputy','Reference'})
        
    end
    
end

% Plot of the Control effort
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
for k = 1
    
    figure
    YLabel={'$ u_1 ~[km/s^2]$',...
        '$u_2 ~[km/s^2]$',...
        '$ u_3 ~[km/s^2] $'};
    for i=1:3
        subplot(3,1,i)
        plot1 = plot( time/Period, Control_effort(:,i));
        ylabel(YLabel(i))
        xlabel('Period')
        grid on
    end
end


%% ------------------ Calculating the action integral ---------------------
for ll = 1
    tic;
    disp('Computing the action integral...');
    % calculate the atuated curve length for each grid of s
    X_temp  = zeros(N,xpoints);
    dX_temp = zeros(N,xpoints);
    cost    = zeros(tpoints,1);
    
    n = 6;
    F = eye(n); F_Aug = blkdiag(F,F);
    % penalty matrix (D matrix in the paper)
    D = diag([Penalty*ones(1,n-M) ones(1,M)]); D_Aug = blkdiag(D,D);
    G = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
    for j = 1:tpoints
        for kk = 1:N
            [X_temp(kk,:), dX_temp(kk,:)] = pdeval(m,x,sol(j,:,kk),x);
        end
        for i = 1:xpoints
            dX = dX_temp(:,i);
            for k = 1:Num_Agent
                idx = 1+n*(k-1):n*k;
                X = X_temp(idx,i);
                if strcmp(DynamicsModel,'THode')
                    f(idx,1) =  NonDimentionalized_THode(x(i),X);
                elseif strcmp(DynamicsModel,'CWode')
                    f(idx,1) = NonDimentionalized_CWode(x(i),X);
                elseif strcmp(DynamicsModel,'NLode')
                    f(idx,1) =  NonDimentionalized_NLode(x(i),X);
                end
            end
            cost(j) = cost(j) + (dX-f)'*G*(dX-f)*(T/(xpoints-1));
        end
    end
    toc;
end
% -------------------------------------------------------------------------
%%
for ll = 1
    dist = zeros(xpoints,1);
    if Num_Agent == 2
        disp('Computing the distance of closest approach...');
        for j = 1: length(x)
            Rho1(1,:) = spline(x,sol(j,:,1),xnew);
            Rho1(2,:) = spline(x,sol(j,:,2),xnew);
            Rho1(3,:) = spline(x,sol(j,:,3),xnew);
            Rho2(1,:) = spline(x,sol(j,:,7),xnew);
            Rho2(2,:) = spline(x,sol(j,:,8),xnew);
            Rho2(3,:) = spline(x,sol(j,:,9),xnew);
            dist(j)   = min(vecnorm(Rho1-Rho2,2,1));
            Rho1      = [];
            Rho2      = [];
        end
    end
end
% -------------------------------------------------------------------------

%% ----------------------- control extraction -----------------------------
for ll = 1
    disp('Extracting the required control...');
    % initialize controls and states
    Crtl = struct('u', cell(1, Num_Agent)); % actual control
    % use spline to interpolate state from steady state pde solution, p
    Pint   = zeros(N,xpoints_new);
    p      = zeros(N,xpoints_new);
    dpdx   = zeros(N,xpoints_new);
    dx_new = diff([xnew(3),xnew,xnew(xpoints_new-2)]);
    for i = 1:N
        Pint(i,:) = spline(x,sol(1,:,i),xnew);
        p(i,:)    = spline(x,sol(end,:,i),xnew);
        % Numerically computing the time derivative of the state
        dp_new    = diff([p(i,3),p(i,:),p(i,xpoints_new-2)]);
        dpdx(i,:) = (dp_new(1:end-1)./dx_new(1:end-1).*dx_new(2:end) ...
            + dp_new(2:end)./dx_new(2:end).*dx_new(1:end-1)) ...
            ./ (dx_new(2:end)+dx_new(1:end-1));
    end
    for i = 1 : length(xnew)
        if strcmp(DynamicsModel,'THode')
            f      = xnew(i);
            h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f    = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv  = [conts1*eye(3)       zeros(3);
                conts2*eye(3)  conts3*eye(3)];
            
            drift = NonDimentionalized_THode(xnew(i),p(:,i));
        elseif strcmp(DynamicsModel,'CWode')
            drift = NonDimentionalized_CWode(xnew(i),p(:,i));
        elseif strcmp(DynamicsModel,'NLode')
            d1 = NonDimentionalized_NLode(xnew(i),p(1:6,i));
            d2 = NonDimentionalized_NLode(xnew(i),p(7:12,i));
            drift = [d1;d2];
        end
        % get [Fc F], note: the F_bar matrix in paper
        Ffull = F;
        B = Ffull(:,M+1:end);
        % Extracting/approximate the required control
        for jj = 1:Num_Agent
            idx = 1+mm*(jj-1):mm*jj;
            Crtl(jj).u(:,i) = (B.'*B)\(B.') * ( dpdx(idx,i) - drift(idx) );
            if strcmp(DynamicsModel,'THode')
                Crtl_Aug = A_inv*[zeros(3,1); Crtl(jj).u(:,i)];
                Dim_Crtl(jj).u(:,i) = Crtl_Aug(4:6);
            elseif strcmp(DynamicsModel,'CWode')
                Dim_Crtl(jj).u(:,i) = Crtl(jj).u(:,i);
            elseif strcmp(DynamicsModel,'NLode')
                Dim_Crtl(jj).u(:,i) = (Rrho*Crtl(jj).u(:,i)).*(1/Tc);
            end
        end
    end

    if exist('u11','file') == 3
        delete u11.c u11.mexmaci64
    end
    if exist('u12','file') == 3
        delete u12.c u12.mexmaci64
    end
    if exist('u13','file') == 3
        delete u13.c u13.mexmaci64
    end
    if exist('u21','file') == 3
        delete u21.c u21.mexmaci64
    end
    if exist('u22','file') == 3
        delete u22.c u22.mexmaci64
    end
    if exist('u23','file') == 3
        delete u23.c u23.mexmaci64
    end
    
    u11 = casadi.interpolant('U','bspline',{xnew}, Dim_Crtl(1).u(1,:));
    u12 = casadi.interpolant('U','bspline',{xnew}, Dim_Crtl(1).u(2,:));
    u13 = casadi.interpolant('U','bspline',{xnew}, Dim_Crtl(1).u(3,:));
    u21 = casadi.interpolant('U','bspline',{xnew}, Dim_Crtl(2).u(1,:));
    u22 = casadi.interpolant('U','bspline',{xnew}, Dim_Crtl(2).u(2,:));
    u23 = casadi.interpolant('U','bspline',{xnew}, Dim_Crtl(2).u(3,:));
    
    u11.generate('u11.c',CasADiopts);
    u12.generate('u12.c',CasADiopts);
    u13.generate('u13.c',CasADiopts);
    u21.generate('u21.c',CasADiopts);
    u22.generate('u22.c',CasADiopts);
    u23.generate('u23.c',CasADiopts);
    mex u11.c -largeArrayDims
    mex u12.c -largeArrayDims
    mex u13.c -largeArrayDims
    mex u21.c -largeArrayDims
    mex u22.c -largeArrayDims
    mex u23.c -largeArrayDims
    clear u11 u12 u13 u21 u22 u23
end
% -------------------------------------------------------------------------

%% ------------------ Integrate the system's trajectory -------------------
for ll = 1
    disp('Integrating the resulting trajectory...');
    if strcmp(DynamicsModel,'THode')
        [~,X_ode45] = ode113(@(t,x)TH_ode(t,x),xnew,X0,options);
    elseif strcmp(DynamicsModel,'CWode')
        [~,X_ode45] = ode113(@(t,x)CW_ode(t,x),xnew,X0,options);
    elseif strcmp(DynamicsModel,'NLode')
        [~,X_ode45_1] = ode113(@(t,x)NL_ode(t,x,1),xnew,X0(1:6),options);
        [~,X_ode45_2] = ode113(@(t,x)NL_ode(t,x,2),xnew,X0(7:12),options);
        X_ode45 = [X_ode45_1 X_ode45_2];
    end
    
    if strcmp(DynamicsModel,'THode')
        for i = 1:length(xnew)
            f      = xnew(i);
            h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f    = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv  = [conts1*eye(3)       zeros(3);
                conts2*eye(3)  conts3*eye(3)];
            X_ode45_LVLH(i,:) = (blkdiag(A_inv,A_inv)*X_ode45(i,:).').';
            Pint(:,i)         = blkdiag(A_inv,A_inv)*Pint(:,i);
        end
    elseif strcmp(DynamicsModel,'CWode')
        X_ode45_LVLH = Tc*X_ode45;
        Pint         = Tc*Pint;
    elseif strcmp(DynamicsModel,'NLode')
        X_ode45_LVLH = [(blkdiag(Rrho,Rrho)*X_ode45_1.').' (blkdiag(Rrho,Rrho)*X_ode45_2.').'];
        Pint         = [blkdiag(Rrho,Rrho)*Pint(1:6,:); blkdiag(Rrho,Rrho)*Pint(7:12,:)];
    end
    
    disp('Done!!!!!');
end

%% ------------------------------------------------------------------------
% ######################## Plotting Results ###############################
% -------------------------------------------------------------------------
%% --------- Plotting the action integral across each iteratons ----------
for ll = 1
    fh2 = figure;
    loglog(t,cost,'-.','LineWidth',4,'Color',c5);
    % plt3 = loglog(t,cost05,'LineWidth',3,'Color',c6);
    % title('Action integral across iteration')
    xlabel('Homotopy iteration')
    ylabel('Action Integral')
    grid on;
    set(gca,'FontSize',30)
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 700 600];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 4320]/res;
    % save figure
    % print(fh2,'ActionIntegralCollisionAvoidance','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%%
for ll = 1
    fh2 = figure;
    hold on
    loglog(t,dist,'-.','Color',c5,'LineWidth',4);
    yl = yline(Threshold,'-.','Minimum allowable distance','LineWidth',3);
    yl.LabelHorizontalAlignment = 'left';
    yl.Color = c7;
    yl.FontSize = 23;
    yl.FontWeight = 'bold';
    % title('Action integral across iteration')
    xlabel('Homotopy iteration')
    ylabel('Distance of closest approach')
    grid on
    ylim([1e-4 1]);
    set(gca,'FontSize',30)
    
    
    % creating the zoom-in inset
    ax=axes;
    set(ax,'units','normalized','position',[.6 .18 .25 .25])
    box(ax,'on')
    loglog(t,dist2,'Color',c6,'LineWidth',4,'parent',ax)
    hold on
    yl = yline(Threshold,'-.','LineWidth',3,'parent',ax);
    yl.Color = c7;
    grid on
    set(ax,'xlim',[13000,1e7],'ylim',[.4497,0.4502])
    
    
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 700 600];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 4320]/res;
    % save figure
    % print(fh2,'DistanceOfClosestApproach','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% ------------------- Plotting the extracted control ---------------------
for ll = 1
    close all
    fh2 = figure;
    Label  = {'Agent 1', 'Agent 2', 'Agent 3', 'Agent 4'};
    for i = 1:Num_Agent
        
        if strcmp(DynamicsModel,'THode')
            Time = (Tc/(2*pi))*xnew/3600;
        else
            Time = Tc*xnew/3600;
            time1 = tspan/3600;
        end
        % set(gca,'XTick',0:pi/2:2*pi)
        % set(gca,'XTickLabel',{'0','$\pi$/2','$\pi$','3$\pi$/2','2$\pi$'})
        subplot(2,2,i)
        plt1 = plot(Time,Dim_Crtl(i).u(1,:),'r','LineWidth',3);
        hold on
        plt2 = plot(Time,Dim_Crtl(i).u(2,:),'b','LineWidth',3);
        plt3 = plot(Time,Dim_Crtl(i).u(3,:),'g','LineWidth',3);
        ylabel({'u [km/s$^{2}$]'})
        % xlabel('time (hr)')
        xlim([Time(1) Time(end)])
        grid on
        % grid minor
        title(Label{i})
        % title('Geometric Control')
        set(gca,'FontSize',30)
        
        % subplot(2,1,2)
        % plot(time1,Control_effort(:,1)','r','LineWidth',3);
        % hold on
        % plot(time1,Control_effort(:,2)','b','LineWidth',3);
        % plot(time1,Control_effort(:,3)','g','LineWidth',3);
        % ylabel('u [km/s$^{2}$]')
        % xlabel('time (hr)')
        % xlim([Time(1) Time(end)])
        % grid on
        % grid minor
        % title('Lyapunov Control')
        
        subplot(2,2,i+2)
        plt1 = plot(xnew,Crtl(i).u(1,:),'r','LineWidth',2.5);
        hold on
        plt2 = plot(xnew,Crtl(i).u(2,:),'b','LineWidth',2.5);
        plt3 = plot(xnew,Crtl(i).u(3,:),'g','LineWidth',2.5);
        ylabel({'Nondimensionalized', 'control'})
        if strcmp(DynamicsModel,'THode')
            xlabel('f [rad]')
        else
            xlabel('$\tau [-]$')
        end
        xlim([xnew(1) xnew(end)])
        grid on;
        title(Label(i))
        set(gca,'FontSize',30)
    end
    hL = legend([plt1, plt2, plt3],...
        {'u1','u2','u3'},'AutoUpdate','off');
    hL.FontSize = 30;
    newPosition = [0.22 0.28 0.1 0.1];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1  );
    
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do not show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 1000 800];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 6080 5020]/res;
    % save figure
    % print(fh2,'RequiredControlCollisionAvoidance','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% ---- 3D plot of the converged solution (Nondimensionalized units) ------
for ll = 1
    % cMat1 = [c5;c6];
    % cMat2 = [c3;c1];
    % cMat0 = [c7;c4];
    %
    % fh2 = figure;
    % hold on
    % plt0  = zeros(Num_Agent);
    % plt1   = zeros(Num_Agent);
    % plt2  = zeros(Num_Agent);
    % for jj = 1:Num_Agent
    %     idx = 1+mm*(jj-1):mm*jj;
    %
    %     plt0(:,jj) = plot3(X_ode45(:,idx(1)),X_ode45(:,idx(2)),X_ode45(:,idx(3)),...
    %         '--','Color',cMat0(jj,:),'LineWidth',2.5);
    %
    %     plt1(:,jj) = plot3(X_rel_0(1,idx(1)),X_rel_0(1,idx(2)),X_rel_0(1,idx(3)),'*',...
    %         'LineWidth',3,...
    %         'MarkerEdgeColor',cMat1(jj,:),...
    %         'MarkerFaceColor',cMat1(jj,:)',...
    %         'MarkerSize',8);
    %
    %
    %     plt2(:,jj) = plot3(X_rel_1(1,idx(1)),X_rel_1(1,idx(2)),X_rel_1(1,idx(3)),'o',...
    %         'LineWidth',3,...
    %         'MarkerEdgeColor',cMat2(jj,:),...
    %         'MarkerFaceColor',cMat2(jj,:),...
    %         'MarkerSize',8);
    %
    %     h5 = plot3(0,0,0,'bo',...
    %         'LineWidth',2,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor','k',...
    %         'MarkerSize',15);
    %     plot3(X_rel_0(:,idx(1)),X_rel_0(:,idx(2)),X_rel_0(:,idx(3)),'Color',cMat1(jj,:),'DisplayName','Initial Traj');
    %     plot3(X_rel_1(:,idx(1)),X_rel_1(:,idx(2)),X_rel_1(:,idx(3)),'Color',cMat2(jj,:));
    %
    %     grid on
    %     xlabel('$\bar{x}$')
    %     ylabel('$\bar{y}$')
    %     zlabel('$\bar{z}$')
    %     title('{\color{black} 3D Relative Trajectory (Nondimensionalized units)}','interpreter','tex')
    % end
    % % r = 1;
    % % [x_sphere, y_sphere, z_sphere] = sphere(100);
    % % h = surf(r*x_sphere, r*y_sphere, r*z_sphere);
    % % set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
    % % axis equal
    %
    % view(68,10)
    % hL = legend([h5, plt1(end,1), plt1(end,2), plt2(end,1), plt2(end,2), plt0(end,1), plt0(end,2)],...
    %     {'Chief Location','Agent1 Init Condition','Agent2 Init Condition',...
    %     'Agent1 End Condition','Agent2 End Condition','Agent1 Transfer Traj',...
    %     'Agent2 Transfer Traj'},'AutoUpdate','off');
    % % Programatically move the Legend
    % newPosition = [0.21 0.72 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    %
    % % set all units inside figure to normalized so that everything is scaling accordingly
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % % do show figure on screen
    % set(fh2, 'visible', 'on')
    % % set figure units to pixels & adjust figure size
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 900 800];
    % % define resolution figure to be saved in dpi
    % res = 500;
    % % recalculate figure size to be saved
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 5580 4320]/res;
    % % save figure
    % % print(fh2,'NondimensionalizedConvergeSolution','-dpng',sprintf('-r%d',res))
    
end
% -------------------------------------------------------------------------

%% --- 3D plot of the converged solution (in the configuration space) -----
for ll = 1
    Num_Agent = 1;
    cMat1 = [c5;c6];
    cMat2 = [c3;c1];
    cMat0 = [c4;c7];
    fh2 = figure;
    view(68,10)
    hold on
    plt0 = zeros(Num_Agent);
    plt1 = zeros(Num_Agent);
    plt2 = zeros(Num_Agent);
    
    for jj = 1:Num_Agent
        idx = 1+mm*(jj-1):mm*jj;
        plt0(:,jj) = plot3(X_ode45_LVLH(:,idx(1)),X_ode45_LVLH(:,idx(2)),X_ode45_LVLH(:,idx(3)),...
            'Color',cMat0(jj,:),'LineWidth',3);
        
        plt1(:,jj) = plot3(X_THpos01(1,idx(1)),X_THpos01(1,idx(2)),X_THpos01(1,idx(3)),'*',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat1(jj,:),...
            'MarkerFaceColor',cMat1(jj,:)',...
            'MarkerSize',12);
        
        plt2(:,jj) = plot3(X_THpos11(1,idx(1)),X_THpos11(1,idx(2)),X_THpos11(1,idx(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat2(jj,:),...
            'MarkerFaceColor',cMat2(jj,:),...
            'MarkerSize',12);
        
        h5 = plot3(0,0,0,'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','k',...
            'MarkerSize',15);
        
        plot3(X_THpos11(:,idx(1)),X_THpos11(:,idx(2)),X_THpos11(:,idx(3)),'--','Color',cMat2(jj,:),'LineWidth',3);
        plot3(X_THpos01(:,idx(1)),X_THpos01(:,idx(2)),X_THpos01(:,idx(3)),'--','Color',cMat1(jj,:),'LineWidth',3);
        
        grid on
        % grid minor
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} e = 0.5}','interpreter','tex')
    end
    
    h1 = plot3(TR_rel(:,1),TR_rel(:,2),TR_rel(:,3),'r','LineWidth',3);
    
    if Num_Agent == 2
        hL = legend([h5, plt1(end,1), plt1(end,2), plt2(end,1), plt2(end,2), plt0(end,1), plt0(end,2)],...
            {'Chief Location','Agent1 Init Condition','Agent2 Init Condition',...
            'Agent1 End Condition','Agent2 End Condition','Agent1 Transfer Traj',...
            'Agent2 Transfer Traj'},'AutoUpdate','off');
    else
        hL = legend([h5, plt1(end,1), plt2(end,1), plt0(end,1),h1],...
            {'Chief Location','Deputy Init Condition',...
            'Deputy End Condition','Geometric Ctrl Traj','Lyapunov Crtl Traj'},...
            'AutoUpdate','off');
    end
    % Programatically move the Legend
    newPosition = [0.25 0.25 0.1 0.1]; %[0.21 0.72 0.1 0.1];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    
    set(gca,'FontSize',30)
    % r = 1;
    % [x_sphere, y_sphere, z_sphere] = sphere(100);
    % h = surf(r*x_sphere, r*y_sphere, r*z_sphere);
    % set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
    % axis equal
    
    
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 1400 1000];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 8580 6320]/res;
    % save figure
    % print(fh2,'ConvergeSolution_AGHFvsLyapunov','-dpng',sprintf('-r%d',res))
    
end
% -------------------------------------------------------------------------


%% -------------------------------------------------------------------------
% ######################## Required Functions #############################
% -------------------------------------------------------------------------
%% --------------------------- PDE for pdepe ------------------------------
function [c,f,s] = mypdexpde(x,t,u,DuDx,k,N) % Define PDE; right-hand-side of AGHF

XChief = full([Chief_x(x);  Chief_y(x);  Chief_z(x);...
               Chief_dx(x); Chief_dy(x); Chief_dz(x)]);
LL = full(L(u,DuDx,k,x,XChief));
CasADiresult = full(dLdx(u,DuDx,k,x,XChief,LL))';
CasADir = [CasADiresult(1:12), CasADiresult(13:24)];

% dL/dx
pLx  = CasADir(:,1); % EL(:,1); %
% dL/d(dot_x)
pLxd = CasADir(:,2); % EL(:,2); % 

f = pLxd;
s = -pLx;
c = ones(N,1);

end
% -------------------------------------------------------------------------

%% --------------------- Initial condition for pdepe ----------------------
function u0 = mypdexic(x, X0, Xf, T)    % Initial condition for AGHF
%%  straight line connecting X0 and Xf
% u0=X0+(Xf-X0)*(x/T);
% % add some sinusoidal deviation to x1 initial guess (can be discarded)
% u0(1,:)=X0(1)+(Xf(1)-X0(1))*(x/T) + 0.01*sin(2*pi*x/T);
% u0(7,:)=X0(7)+(Xf(7)-X0(7))*(x/T) + 0.01*sin(2*pi*x/T);
% u0(1:6,:)=X0(1:6)+(Xf(1:6)-X0(1:6))*(x/T);
% % add some sinusoidal deviation to x1 initial guess (can be discarded)
% u0(1,:)=X0(1)+(Xf(1)-X0(1))*(x/T) + 0.01*sin(2*pi*x/T);
% 
%% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;

u0 = X0*cos(freq1*x/T) ...
    + Xf*sin(freq1*x/T) ...
    - (X0+Xf)*sin(freq2*x/T);

end
% -------------------------------------------------------------------------

%% --------------------- Boundary condition for pdepe ---------------------
function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t, X0, Xf, N) % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end
% -------------------------------------------------------------------------


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!   Nondimetionalized Dynamics   !!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [Xdot] = NonDimentionalized_THode(t,X)
global e_chief
k = 1+e_chief*cos(t);
A = [zeros(3) eye(3);
    3/k  0  0   0   2  0;
    0    0  0   -2  0  0;
    0    0  -1  0   0  0];
[l,~] = size(X);
X = reshape(X,6,l/6);
Xdot = A*X;
Xdot = Xdot(:);
end
% -------------------------------------------------------------------------

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!   CasADi Functions   !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [Xdot] = NonDimentionalized_THodeCasADi(t,X,XChief)
global e_chief
k = 1+e_chief*cos(t);
A = [zeros(3) eye(3);
    3/k  0  0   0   2  0;
    0    0  0   -2  0  0;
    0    0  -1  0   0  0];
[l,~] = size(X);
X = reshape(X,6,l/6);
Xdot = A*X;
Xdot = Xdot(:);
end
% -------------------------------------------------------------------------



