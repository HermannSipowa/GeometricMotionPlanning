%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hermann Kaptui Sipowa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear
close all
start_up
format long e
clc
str1 = '/Users/';
str2 = getenv('USER');
str3 = '/Desktop/casadi-osx-matlabR2015a-v3.5.5';
addpath(append(str1,str2,str3))
import casadi.*
global mu_Earth Target Chasser M ...
       Tc Num_Agent ae JD AU % iteration InitalCondition
ae        = 6378.136;
mu_Earth  = 3.986004415E5;
JD        = 2456296.25;
AU        = 149597870.7;
Chasser   = Spacecraft([10 10 1 .25 2 2 0.2 0.5]); % Deputy characteristic
Target    = Spacecraft([0 0 2 .25 1.5 .5 0.2 0.5]); % Target characteristic
Num_Agent = 1; % Number of agents in the system
SimTime   = 1.25/4;
CasADiopt = struct('main', true, 'mex', true);
M   = 3; 
N   = 13; 
c1  = rgb('DarkGreen');
c2  = rgb('Gold');
c3  = rgb('Lime');
c4  = rgb('DarkOrange');
c5  = rgb('DarkBlue');
c6  = rgb('Red');
c7  = rgb('Purple');
c8  = rgb('Bisque');
c9  = rgb('Orange');
c10 = rgb('DarkGray');
c11 = rgb('Teal');

%% Calculating the boundaries conditions
%=======================================%
% Specifying the chief's orbit
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a_chief      = 2e4; % Semi-major axis in Km
e_chief      = 0.5; % Eccentricity
inc_chief    = 50;  % Inclination in deg0
BigOmg_chief = 10;  % RAAN in deg
LitOmg_chief = 10;  % AOP in deg
M_chief1     = 0;   % Initial mean anomaly

%% Integrating the chief's trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Period  = (2*pi)*sqrt(a_chief^3/mu_Earth);
Tc      = 4*Period;
IntTime = SimTime*Tc;
Tsize   = 1e4;
tspan   = linspace(0,IntTime,Tsize);
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-30);
mu_dot  = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
Tnorm = tspan/Tc;
Ttraj1 = linspace(0,Period,Tsize); TT1 = Ttraj1/Tc;
Ttraj2 = IntTime+linspace(0,Period,Tsize); TT2 = Ttraj2/Tc;

inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
LitOmg_chief = deg2rad(LitOmg_chief); M_chief1 = deg2rad(M_chief1);
COE1 = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief1];
[Position_target1,Velocity_target1] = COEstoRV(COE1,mu_Earth);
X0_Chief1  = [Position_target1; Velocity_target1];
[~, Xnom]  = ode113(@(t,X)ChiefMotionODE(t,X,Target),tspan,X0_Chief1,options);

X0_Chief2  = Xnom(end,:)'; 
[~, Xnom0] = ode113(@(t,X)ChiefMotionODE(t,X,Target),Ttraj1,X0_Chief1,options);
[~, Xnomf] = ode113(@(t,X)ChiefMotionODE(t,X,Target),Ttraj2,X0_Chief2,options);

%% ------- Create MX files to interpolate the chief's states (Xnom) -------
for ll = 1
    if exist('Chief_x','file') == 3
        delete Chief_x.c Chief_x.mexmaci64
    end
    if exist('Chief_y','file') == 3
        delete Chief_y.c Chief_y.mexmaci64
    end
    if exist('Chief_z','file') == 3
        delete Chief_z.c Chief_z.mexmaci64
    end
    if exist('Chief_dx','file') == 3
        delete Chief_dx.c Chief_dx.mexmaci64
    end
    if exist('Chief_dy','file') == 3
        delete Chief_dy.c Chief_dy.mexmaci64
    end
    if exist('Chief_dz','file') == 3
        delete Chief_dz.c Chief_dz.mexmaci64
    end
    
    Chief_x  = casadi.interpolant('Chief_x','bspline',{Tnorm},Xnom(:,1)');
    Chief_y  = casadi.interpolant('Chief_y','bspline',{Tnorm},Xnom(:,2)');
    Chief_z  = casadi.interpolant('Chief_z','bspline',{Tnorm},Xnom(:,3)');
    Chief_dx = casadi.interpolant('Chief_dx','bspline',{Tnorm},Xnom(:,4)');
    Chief_dy = casadi.interpolant('Chief_dy','bspline',{Tnorm},Xnom(:,5)');
    Chief_dz = casadi.interpolant('Chief_dz','bspline',{Tnorm},Xnom(:,6)');
    
    Chief_x.generate('Chief_x.c',CasADiopt);
    Chief_y.generate('Chief_y.c',CasADiopt);
    Chief_z.generate('Chief_z.c',CasADiopt);
    Chief_dx.generate('Chief_dx.c',CasADiopt);
    Chief_dy.generate('Chief_dy.c',CasADiopt);
    Chief_dz.generate('Chief_dz.c',CasADiopt);
    mex Chief_x.c -largeArrayDims
    mex Chief_y.c -largeArrayDims
    mex Chief_z.c -largeArrayDims
    mex Chief_dx.c -largeArrayDims
    mex Chief_dy.c -largeArrayDims
    mex Chief_dz.c -largeArrayDims
    clear Chief_x Chief_y Chief_z Chief_dx Chief_dy Chief_dz
    
    
    
    if exist('Chief_x0','file') == 3
        delete Chief_x0.c Chief_x0.mexmaci64
    end
    if exist('Chief_y0','file') == 3
        delete Chief_y0.c Chief_y0.mexmaci64
    end
    if exist('Chief_z0','file') == 3
        delete Chief_z0.c Chief_z0.mexmaci64
    end
    if exist('Chief_dx0','file') == 3
        delete Chief_dx0.c Chief_dx0.mexmaci64
    end
    if exist('Chief_dy0','file') == 3
        delete Chief_dy0.c Chief_dy0.mexmaci64
    end
    if exist('Chief_dz0','file') == 3
        delete Chief_dz0.c Chief_dz0.mexmaci64
    end
    
    Chief_x0  = casadi.interpolant('Chief_x0','bspline',{TT1},Xnom0(:,1)');
    Chief_y0  = casadi.interpolant('Chief_y0','bspline',{TT1},Xnom0(:,2)');
    Chief_z0  = casadi.interpolant('Chief_z0','bspline',{TT1},Xnom0(:,3)');
    Chief_dx0 = casadi.interpolant('Chief_dx0','bspline',{TT1},Xnom0(:,4)');
    Chief_dy0 = casadi.interpolant('Chief_dy0','bspline',{TT1},Xnom0(:,5)');
    Chief_dz0 = casadi.interpolant('Chief_dz0','bspline',{TT1},Xnom0(:,6)');
    
    Chief_x0.generate('Chief_x0.c',CasADiopt);
    Chief_y0.generate('Chief_y0.c',CasADiopt);
    Chief_z0.generate('Chief_z0.c',CasADiopt);
    Chief_dx0.generate('Chief_dx0.c',CasADiopt);
    Chief_dy0.generate('Chief_dy0.c',CasADiopt);
    Chief_dz0.generate('Chief_dz0.c',CasADiopt);
    mex Chief_x0.c -largeArrayDims
    mex Chief_y0.c -largeArrayDims
    mex Chief_z0.c -largeArrayDims
    mex Chief_dx0.c -largeArrayDims
    mex Chief_dy0.c -largeArrayDims
    mex Chief_dz0.c -largeArrayDims
    clear Chief_x0 Chief_y0 Chief_z0 Chief_dx0 Chief_dy0 Chief_dz0
    
    
    if exist('Chief_xf','file') == 3
        delete Chief_xf.c Chief_xf.mexmaci64
    end
    if exist('Chief_yf','file') == 3
        delete Chief_yf.c Chief_yf.mexmaci64
    end
    if exist('Chief_zf','file') == 3
        delete Chief_zf.c Chief_zf.mexmaci64
    end
    if exist('Chief_dxf','file') == 3
        delete Chief_dxf.c Chief_dxf.mexmaci64
    end
    if exist('Chief_dyf','file') == 3
        delete Chief_dyf.c Chief_dyf.mexmaci64
    end
    if exist('Chief_dzf','file') == 3
        delete Chief_dzf.c Chief_dzf.mexmaci64
    end
    
    Chief_xf  = casadi.interpolant('Chief_xf','bspline',{TT2},Xnomf(:,1)');
    Chief_yf  = casadi.interpolant('Chief_yf','bspline',{TT2},Xnomf(:,2)');
    Chief_zf  = casadi.interpolant('Chief_zf','bspline',{TT2},Xnomf(:,3)');
    Chief_dxf = casadi.interpolant('Chief_dxf','bspline',{TT2},Xnomf(:,4)');
    Chief_dyf = casadi.interpolant('Chief_dyf','bspline',{TT2},Xnomf(:,5)');
    Chief_dzf = casadi.interpolant('Chief_dzf','bspline',{TT2},Xnomf(:,6)');
    
    Chief_xf.generate('Chief_xf.c',CasADiopt);
    Chief_yf.generate('Chief_yf.c',CasADiopt);
    Chief_zf.generate('Chief_zf.c',CasADiopt);
    Chief_dxf.generate('Chief_dxf.c',CasADiopt);
    Chief_dyf.generate('Chief_dyf.c',CasADiopt);
    Chief_dzf.generate('Chief_dzf.c',CasADiopt);
    mex Chief_xf.c -largeArrayDims
    mex Chief_yf.c -largeArrayDims
    mex Chief_zf.c -largeArrayDims
    mex Chief_dxf.c -largeArrayDims
    mex Chief_dyf.c -largeArrayDims
    mex Chief_dzf.c -largeArrayDims
    clear Chief_xf Chief_yf Chief_zf Chief_dxf Chief_dyf Chief_dzf
end 
% -------------------------------------------------------------------------

%% --------------- Setting the deputy initial conditions ------------------
for ll = 1
    r0 = Xnom0(1,1:3)'; v0 = Xnom0(1,4:6)';
    COE0 = RVtoCOEs(r0,v0,mu_Earth); f1 = COE0(6);
    q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
    theta_chief1 = f1+LitOmg_chief;
    OE_chief1 = [a_chief, theta_chief1, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    AMap1 = ForwardMapping(OE_chief1, mu_Earth); % Linear mapping matrix
    
    rf = Xnomf(1,1:3)'; vf = Xnomf(1,4:6)';
    COEf = RVtoCOEs(rf,vf,mu_Earth); f2 = COEf(6);
    theta_chief2 = f2+LitOmg_chief;
    OE_chief2 = [a_chief, theta_chief2, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    AMap2 = ForwardMapping(OE_chief2, mu_Earth); % Linear mapping matrix
    
    % Specify the final relative-orbital elements of the deputies in both orbits
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    dela = 0;
    dele1 = 1/(5*a_chief);
    deli1 = -1/(3*a_chief);
    delLitOmg = -2*pi*1E-8;
    delBigOmg = 0;
    delM = pi*1E-6;
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_1 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    
    dele1 = -1/(2*a_chief);
    deli1 = -1/(5*a_chief);
    delLitOmg = -5*pi*1E-7;
    delBigOmg = pi*1E-6;
    delM = 0;
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
    
    dele2 = 1/(8*a_chief);
    deli2 = 2/(4*a_chief);
    delBigOmg = -pi*1E-8;
    delM = pi*1E-6;
    delLitOmg = pi*1E-8;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    
    delq1_2 = dele2*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_2 = dele2*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    delCOE_4 = [dela, deltheta, deli2, delq1_2, delq2_2, delBigOmg]';
    
    % Integrating the deputy initial and final trajectories
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    X01 = AMap1*delCOE_1;
    X02 = AMap1*delCOE_2;
    Xf1 = AMap2*delCOE_3;
    Xf2 = AMap2*delCOE_4;
    
    qr_chasser0 = [6; 4; 9; 1];
    qr_chasser0 = qr_chasser0/norm(qr_chasser0);
    Omega_chaser0 = 5e-2*[deg2rad(-.9) deg2rad(.8) -deg2rad(.7)]'; % Angular velocity in inertial frame
    CN = Quaternion_to_DCM(qr_chasser0); Omega_chaser_B0 = CN*Omega_chaser0;
    X0 = [qr_chasser0; Omega_chaser0; X01];
    qr_chasserf = [9; 5; 1; 8];
    qr_chasserf = qr_chasserf/norm(qr_chasserf);
    Omega_chaserf = 3e-3*[deg2rad(-.3) deg2rad(.3) -deg2rad(.1)]';
    CN = Quaternion_to_DCM(qr_chasserf); Omega_chaser_Bf = CN*Omega_chaserf;
    Xf = [qr_chasserf; Omega_chaserf; Xf1];
end 
% -------------------------------------------------------------------------

%% ------- Nondimentionalized the inital condition of the problem ---------
for ll = 1
    options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-15);
    [~, Xrel0] = ode113(@(t,x)NonDimentionalized_NLode1(t,x,Target),TT1,X0(8:13),options);
    [~, Xrel1] = ode113(@(t,x)NonDimentionalized_NLode2(t,x,Target),TT2,Xf(8:13),options);
    [~, Xrel]  = ode113(@(t,x)NonDimentionalized_FlatPlateOde(t,x,Target,Chasser),Tnorm,X0,options);
end
% -------------------------------------------------------------------------
%%
for ll = 1
    % figure
    % hold on
    % plot3(Xrel0(:,1),Xrel0(:,2),Xrel0(:,3),'b','LineWidth',3)
    % plot3(Xrel1(:,1),Xrel1(:,2),Xrel1(:,3),'r','LineWidth',3)
    % plot3(Xrel(:,8),Xrel(:,9),Xrel(:,10),'g','LineWidth',3)
    % plot3(Xrel0(1,1),Xrel0(1,2),Xrel0(1,3),'*',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c5,...
    %     'MarkerFaceColor',c5,...
    %     'MarkerSize',8);
    % plot3(Xrel1(1,1),Xrel1(1,2),Xrel1(1,3),'*',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c6,...
    %     'MarkerFaceColor',c6,...
    %     'MarkerSize',8);
    % plot3(0,0,0,'bo',...
    %     'LineWidth',2,...
    %     'MarkerEdgeColor','k',...
    %     'MarkerFaceColor','k',...
    %     'MarkerSize',15);
    % view([64,26])
    % figure
    % plot(Tnorm,vecnorm(Xrel(:,8:10),2,2))
end

%% ------------------------ AutoDiff using CasADi -------------------------
for ll = 1
    import casadi.*
    t        = SX.sym('t'); 
    Xaug     = SX.sym('q1',13);  % Xaug  = [q1;omega1;rho1];
    dXaug    = SX.sym('q1',13); % dXaug = [dq1;domega1;drho1];
    XChief   = SX.sym('XChief',6);
    Penalty1 = SX.sym('Penalty1');
    Penalty2 = SX.sym('Penalty2');
    Penalty3 = SX.sym('Penalty3');
    
    Ichaser_inv = 1./(Chasser.Moment_Of_Inertia_Calculator().');
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
            F = eye(N); F(5:7,5:7) = diag(Ichaser_inv); F_Aug = F;
            % penalty matrix (D matrix in the paper)
            D = diag([Penalty1*ones(1,4), ones(1,3), Penalty2*ones(1,3), Penalty3*ones(1,3)]);
            D_Aug = D;
        else
            F_Aug = blkdiag(F_Aug,F);
            D_Aug = blkdiag(D_Aug,D);
        end
    end
    
    % Defining the metric, before adding the state constraint barrier functions
    H = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
    
    
    % % Defining the drift vector field
    % for i = 1:Num_Agent
    %     idx = 1+N*(i-1):N*i;
    %     drift(idx,1) = NonDimentionalized_FlatPlateOdeCasADi(t,Xaug(idx),XChief);
    % end
    drift = NonDimentionalized_FlatPlateOdeCasADi(t,Xaug,Target,Chasser,XChief);
    
    % ------------------- state constraint barrier function ------------------
    % Collecting relevant quantities
    r_target = XChief(1:3); v_target = XChief(4:6);
    rt_norm = norm(r_target);
    TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame
    rho = Xaug(8:10);
    rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
    r_chasser = TN\rc; % ECI position of the chasser
    [~,Costheta] = F_SPR_FlatPlate(t,r_chasser,Chasser,Xaug(1:4));
    
    Threshold = 0.05;
    kb = 1e0; % 1e-7; % barrier function gain
    pb = 1;   % order of barrier function
    B = kb/(Costheta - Threshold)^pb;
    B = [];
    
    % -------------------- Metric and curve Length ----------------------------
    % the metric with state constraint barier functions
    G = (sum(B)+1)*H;
    % actuated curve length
    L = Function('L',{Xaug,dXaug,Penalty1,Penalty2,Penalty3,t,XChief},...
        {(dXaug - drift).' * G * (dXaug - drift)},...
        {'X','dX','k1','k2','k3','t','XChief'},...
        {'CurveLength'});
    dLdx = L.jacobian();
    fdrift = Function('fdrift', {t,Xaug,XChief},...
        {drift});
    
    % Rewrite the functions in c files
    L.generate('L.c',CasADiopt);
    dLdx.generate('dLdx.c',CasADiopt);
    fdrift.generate('fdrift.c',CasADiopt);
    
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
%% --------------------------- AGHF parameters ----------------------------
for ll = 1
    clc
    smax = 5e5; % 5e4; % % Need to find an analytic way to solve for this value
    % # of grids in t
    tgrids = 1e3;
    % # of grids in s
    sgrids = 1e3;
    % # of inputs
    M = 3;  
    N = length(X0);
    % penalty value (the lambda in paper)
    Penalty1 = 1e-3; Penalty2 = 5e1; Penalty3 = 5e1;
    % Penalty1 = 1e-1; Penalty2 = 5e1; Penalty3 = 5e1;
    % tolerance of the integrator
    opts = odeset('AbsTol',1e-14);
    % opts = odeset('RelTol',2.22045e-9,'AbsTol',2.22045e-20);
    % Setting the boundary condtions and the intial condition
    tmax = smax; xpoints = tgrids; tpoints = sgrids; m = 0;
    t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of s interval in log scale
    T = SimTime; % motion duration
    x = linspace(0,T,xpoints);
end
% -------------------------------------------------------------------------

%% ------------------------- Solving the AGHF -----------------------------
for ll = 1
    % % solve for trajectory, see "AGHF" function file for implementation details
    % tic;
    % disp('solving PDE...');
    % % Solve AGHF
    % system('caffeinate -dims &');
    % sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,Penalty1,Penalty2,Penalty3,N),...
    %     @(x) mypdexic(x,X0,Xf,T),...
    %     @(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t,X0,Xf,N),...
    %     x,t); % The solution "sol" is of form sol(t,x,i)
    % system('killall caffeinate');
    % toc;
    % save sol2.mat sol -v7.3
    Sol = load('sol2.mat'); % matfile('indexYCentralized.mat');
    sol = Sol.sol;
end
% -------------------------------------------------------------------------

%% ------------------ Calculating the action integral ---------------------
for ll = 1
    tic;
    disp('Computing the action integral...');
    % calculate the atuated curve length for each grid of s
    X_temp  = zeros(N,xpoints);
    dX_temp = zeros(N,xpoints);
    cost    = zeros(tpoints,1);
    
    % F = eye(N); Ichaser_inv = 1./Chasser.Moment_Of_Inertia_Calculator().';
    % F(5:7,5:7) = diag(Ichaser_inv);
    % D = diag([Penalty1*ones(1,4) ones(1,3) Penalty2*ones(1,3) Penalty3*ones(1,3)]);
    
    % % Collecting relevant quantities
    % r_target = XChief(1:3); v_target = XChief(4:6);
    % rt_norm = norm(r_target);
    % TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame
    % rho = Xaug(8:10);
    % rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
    % r_chasser = TN\rc; % ECI position of the chasser
    % [~,Costheta] = F_SPR_FlatPlate(t,r_chasser,Chasser,Xaug(1:4));
    %
    % Threshold = 0.02;
    % kb = 1e1; % 1e-7; % barrier function gain
    % pb = 1;   % order of barrier function
    % B = kb/(Costheta - Threshold)^pb;
    % B = []; H = (F.')^(-1)*D*F^(-1); G = (sum(B)+1)*H;
    
    
    
    
    
    for j = 1:tpoints
        for kk = 1:N
            [X_temp(kk,:), dX_temp(kk,:)] = pdeval(m,x,sol(j,:,kk),x);
        end
        
        for i = 1:xpoints
            % dX = dX_temp(:,i);
            % X  = X_temp(:,i);
            % f  = NonDimentionalized_FlatPlateOde(x(i),X,Target,Chasser);
            % cost(j) = cost(j) + (dX-f).'*G*(dX-f)*(T/(xpoints-1));
            
            X_chief = full([Chief_x(x(i));  Chief_y(x(i));  Chief_z(x(i));...
                       Chief_dx(x(i)); Chief_dy(x(i)); Chief_dz(x(i))]);
            DuDx    = dX_temp(:,i);
            u       = X_temp(:,i);
            LL      = full(L(u,DuDx,Penalty1,Penalty2,Penalty3,x(i),X_chief));
            cost(j) = cost(j) + LL;
            
        end
    end
    toc;
end
% -------------------------------------------------------------------------

%% ----------------------- control extraction -----------------------------
for ll = 1
    disp('Extracting the required control...');
    % initialize controls and states
    Crtl   = struct('u', cell(1, Num_Agent)); % actual control
    % use spline to interpolate state from steady state pde solution, p
    qvec   = zeros(4,Tsize);
    uTotal = zeros(N,Tsize);
    u      = nan(N,xpoints);
    du     = nan(N,xpoints);
    Xval   = zeros(N,Tsize);
    dXdt   = zeros(N,Tsize);
    for kk = 1:N
        [Xval(kk,:), dXdt(kk,:)] = pdeval(m,x,sol(end,:,kk),Tnorm);
    end
    
    dudx = nan(N,length(x));
    DuDx  = nan(N,Tsize);
    xd = diff([x(3),x,x(N-2)]);  % <-- Corrected
    for ii = 1:13
        u = sol(end,:,ii);
        ud    = diff([u(3),u,u(N-2)]);  % <-- Corrected
        dudx(ii,:)  = (ud(1:end-1)./xd(1:end-1).*xd(2:end) ...
            + ud(2:end)./xd(2:end).*xd(1:end-1)) ...
            ./ (xd(2:end)+xd(1:end-1));
        DuDx(ii,:) = makima(x, dudx(ii,:), Tnorm);
    end

    
    
    for jj = 1:Num_Agent
        for i = 1:Tsize
            X = Xval(:,i);
            drift = NonDimentionalized_FlatPlateOde(Tnorm(i),X,Target,Chasser);
            % get [Fc F], note: the F_bar matrix in paper
            Ichaser_inv = diag(1./(Chasser.Moment_Of_Inertia_Calculator().'));
            B = [zeros(4,3); Ichaser_inv; zeros(6,3)];
            % Extracting/approximate the required control
            Crtl(jj).u(:,i) = (1/Tc).*((B.'*B)\(B.') * ( DuDx(:,i) - drift ));
            uTotal(:,i) = (1/Tc).* ( DuDx(:,i) - drift );
            
            qvec(:,i) = Xval(1:4,i);
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
        
    u11 = casadi.interpolant('U','bspline',{Tnorm}, Crtl(1).u(1,:));
    u12 = casadi.interpolant('U','bspline',{Tnorm}, Crtl(1).u(2,:));
    u13 = casadi.interpolant('U','bspline',{Tnorm}, Crtl(1).u(3,:));
    % u21 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(1,:));
    % u22 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(2,:));
    % u23 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(3,:));
    
    u11.generate('u11.c',CasADiopt);
    u12.generate('u12.c',CasADiopt);
    u13.generate('u13.c',CasADiopt);
    % u21.generate('u21.c',CasADiopts);
    % u22.generate('u22.c',CasADiopts);
    % u23.generate('u23.c',CasADiopts);
    mex u11.c -largeArrayDims
    mex u12.c -largeArrayDims
    mex u13.c -largeArrayDims
    % mex u21.c -largeArrayDims
    % mex u22.c -largeArrayDims
    % mex u23.c -largeArrayDims
    clear u11 u12 u13 % u21 u22 u23
    %========================================%
end
% -------------------------------------------------------------------------

%% 
for ll = 1
    
    Ucrtl1 = nan(3,Tsize);
    Ucrtl2 = nan(3,Tsize);
    qcheck1 = nan(4,Tsize);
    qcheck2 = nan(4,Tsize);
    ChiefHnorm1 = nan(1,Tsize);
    ChiefHnorm2 = nan(1,Tsize);
    XdotCasdi   = nan(6,Tsize);
    XdotCheck   = nan(6,Tsize);
    u0 = nan(N,Tsize);
    for i = 1:Tsize
        tt = Tnorm(i);
        X = Xval(8:13,i);
        X_chief = full([Chief_x(tt);  Chief_y(tt);  Chief_z(tt);...
            Chief_dx(tt); Chief_dy(tt); Chief_dz(tt)]);
    
    
        qcheck1(:,i) = [Xval(1,i); Xval(2,i); Xval(3,i); Xval(4,i)]; % qcheck1(:,i) = qcheck1(:,i)/norm(qcheck1(:,i));
        [Xdot, DelUsrp, U_eci_Chasser, U_eci_Target,Costf] = ControlsAccelerationCheck(tt,X,Target,Chasser,qcheck1(:,i),X_chief);
        Ucrtl1(:,i) = DelUsrp;
    
        qcheck2(:,i) = Xrel(i,1:4).'; % qcheck2(:,i) = qcheck2(:,i)/norm(qcheck2(:,i));
        X = Xrel(i,8:13).';
        [Xdot1, DelUsrp, U_eci_Chasser, U_eci_Target,Costi]  = ControlsAccelerationCheck(tt,X,Target,Chasser,qcheck2(:,i),X_chief);
        Ucrtl2(:,i) = DelUsrp;
        
        Xdot2 = fdrift(tt,Xrel(i,:).',X_chief);
        XdotCasdi(:,i) = Xdot2(8:13);
        XdotCheck(:,i) = Xdot1;
        
        
        u0(:,i) = mypdexic(tt,X0,Xf,T);
        Ci(i) = Costi;
        Cf(i) = Costf;
    
    end
    
    
    figure
    Time = Tc*Tnorm/3600;
    plot(Time,Ci,'b--')
    hold on 
    plot(Time,Cf,'r')
    grid on 
    grid minor
    
    figure
    plt1 = plot(Time,Ucrtl1(1,:),'r','LineWidth',2.5);
    hold on
    plt2 = plot(Time,Ucrtl1(2,:),'b','LineWidth',2.5);
    plt3 = plot(Time,Ucrtl1(3,:),'g','LineWidth',2.5);
    grid on
    hL1 = legend([plt1, plt2, plt3],...
            {'$u_1$','$u_2$','$u_3$'},'AutoUpdate','off');
    ylabel('SRP acceleration [N]')
    xlabel('time (hr)')
    xlim([Time(1) Time(end)])
    title('Relative SRP acceleration from the Converged trajectory')
    
    figure
    plt1 = plot(Time,Ucrtl2(1,:),'r','LineWidth',2.5);
    hold on
    plt2 = plot(Time,Ucrtl2(2,:),'b','LineWidth',2.5);
    plt3 = plot(Time,Ucrtl2(3,:),'g','LineWidth',2.5);
    grid on
    hL2 = legend([plt1, plt2, plt3],...
            {'$u_1$','$u_2$','$u_3$'},'AutoUpdate','off');
    ylabel('SRP acceleration [N]')
    xlabel('time (hr)')
    xlim([Time(1) Time(end)])
    title('Relative SRP acceleration from the Initial guess')
    
    
    figure
    plt1 = plot(Time,uTotal(11,:),'r','LineWidth',2.5);
    hold on
    plt2 = plot(Time,uTotal(12,:),'b','LineWidth',2.5);
    plt3 = plot(Time,uTotal(13,:),'g','LineWidth',2.5);
    grid on
    hL3 = legend([plt1, plt2, plt3],...
            {'$u_1$','$u_2$','$u_3$'},'AutoUpdate','off');
    ylabel('SRP acceleration PDEPE [N]')
    xlabel('time (hr)')
    xlim([Time(1) Time(end)])
    title('Initial guess')
    
    figure
    plt1 = plot(Time,uTotal(11,:)-Ucrtl1(1,:),'r','LineWidth',2.5);
    hold on
    plt2 = plot(Time,uTotal(12,:)-Ucrtl1(2,:),'b','LineWidth',2.5);
    plt3 = plot(Time,uTotal(13,:)-Ucrtl1(3,:),'g','LineWidth',2.5);
    grid on
    hL4 = legend([plt1, plt2, plt3],...
            {'$u_1$','$u_2$','$u_3$'},'AutoUpdate','off');
    ylabel('Additional acceleration needed [N]')
    xlabel('time (hr)')
    xlim([Time(1) Time(end)])
    title('Initial guess')
    
    
    % figure
    % for j = 1:6
    %     subplot(3,2,j)
    %     plot(Time,XdotCasdi(j,:),'b')
    %     hold on
    %     plot(Time,XdotCheck(j,:),'r--')
    %     grid on
    % end

    
    % figure
    % for j = 1:4
    %     subplot(2,2,j)
    %     plot(Time,qcheck1(j,:),'b--')
    %     hold on
    %     plot(Time,qcheck2(j,:),'r')
    %     grid on
    % end
    
    
    %
    % figure
    % hold on
    % plot3(Xrel0(:,1),Xrel0(:,2),Xrel0(:,3),'b','LineWidth',3)
    % plot3(Xrel1(:,1),Xrel1(:,2),Xrel1(:,3),'r','LineWidth',3)
    % plot3(u0(8,:),u0(9,:),u0(10,:),'g','LineWidth',3)
    % plot3(Xrel0(1,1),Xrel0(1,2),Xrel0(1,3),'*',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c5,...
    %     'MarkerFaceColor',c5,...
    %     'MarkerSize',8);
    % plot3(Xrel1(1,1),Xrel1(1,2),Xrel1(1,3),'*',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c6,...
    %     'MarkerFaceColor',c6,...
    %     'MarkerSize',8);
    % plot3(0,0,0,'bo',...
    %     'LineWidth',2,...
    %     'MarkerEdgeColor','k',...
    %     'MarkerFaceColor','k',...
    %     'MarkerSize',15);
    % grid on
    % view([64,26])

end

%% ------------------ Integrate the system's trajectory -------------------
for ll = 1
    disp('Integrating the resulting trajectory...');
    [~,X_ode45] = ode113(@(t,Xaug)NonDimentionalized_FlatPlateControlledOde(t,Xaug,Target,Chasser),Tnorm,X0,options);
    disp('Done!!!!!');
end

%% ------------------------------------------------------------------------
% ######################## Plotting Results ###############################
% -------------------------------------------------------------------------
%% Junk plots for checking purposes only
for ll = 1
    % figure
    % for i = 1:13
    %     subplot(4,4,i)
    %     plot(tspan,uTotal(i,:),'Color',c1)
    %     grid on
    % end
    %
    % figure
    % for i = 1:13
    %     subplot(4,4,i)
    %     plot(tspan,X_ode45(:,i)-Xrel(:,i),'r')
    %     grid on
    % end
    %
    % close all
    fh2 = figure;
    LaBeL1 = {'$q_{0}$','$q_{1}$','$q_{2}$','$q_{3}$'};
    for i = 1:4
        subplot(2,2,i)
        plt1 = plot(Time,Xval(i,:),'Color',c6);
        hold on
        plt2 = plot(Time,X_ode45(:,i).','--','Color',c11);
        
        ylabel(LaBeL1{i})
        xlabel('t [days]')
        grid on
        grid minor
        % set(gca,'FontSize',25)
        
        newPosition = [0.796 0.13 0.1 0.1];
        hL = legend([plt1, plt2],...
            {'PDEPE','ODE45'},'AutoUpdate','off', 'Position', newPosition);
    end
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 1000 600];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 3320]/res;
    % save figure
    % print(fh2,'Quaternion_Tracking_6DoF','-dpng',sprintf('-r%d',res))
    
    
    
    
    
    
    % figure
    % LaBeL2 = {'$\omega_{1}$','$\omega_{2}$','$\omega_{3}$',...
    %     'r$_{1}$','r$_{2}$','r$_{3}$','v$_{1}$','v$_{2}$','v$_{3}$'};
    % for i = 1:9
    %     j = i + 4;
    %     subplot(3,3,i)
    %     plot(tspan,Xval(j,:) - X_ode45(:,j).','Color',c5)
    %     grid on
    %     grid minor
    %     ylabel(LaBeL2{i})
    % end
    %
    % figure
    % for i = 1:13
    %     subplot(4,4,i)
    %     plot(tspan,p(i,:)-Xrel(:,i).','Color',c4)
    %     grid on
    % end
    %
    % figure
    % for i = 1:13
    %     subplot(4,4,i)
    %     plot(tspan,p(i,:),'Color',c6)
    %     grid on
    % end
end

%% --------- Plotting the action integral across each iteratons ----------
for ll = 1
    fh2 = figure;
    loglog(t,cost,'Color',c5,'LineWidth',2.5);
    title('Action integral across iteration')
    xlabel('Homotopy iteration')
    ylabel('Length of the actuated curve')
    grid on;
    % % set all units inside figure to normalized so that everything is scaling accordingly
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % % do show figure on screen
    % set(fh2, 'visible', 'on')
    % % set figure units to pixels & adjust figure size
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 700 600];
    % % define resolution figure to be saved in dpi
    % res = 500;
    % % recalculate figure size to be saved
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 5580 4320]/res;
    % % save figure
    % % print(fh2,'ActionIntegral','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% ------------------- Plotting the extracted control ---------------------
for ll = 1
    fh2 = figure;
    Label  = {'Agent 1', 'Agent 2', 'Agent 3', 'Agent 4'};
    for i = 1:Num_Agent
        Time = Tc*Tnorm/3600;
        plt1 = plot(Time,Crtl(i).u(1,:),'r','LineWidth',2.5);
        hold on
        plt2 = plot(Time,Crtl(i).u(2,:),'b','LineWidth',2.5);
        plt3 = plot(Time,Crtl(i).u(3,:),'g','LineWidth',2.5);
        ylabel(Label(i))
        xlabel('time (hr)')
        xlim([Time(1) Time(end)])
        grid on;
        title('Required Torque [N-m]')
    end
    hL = legend([plt1, plt2, plt3],...
        {'$\tau_1$','$\tau_2$','$\tau_3$'},'AutoUpdate','off');
    % hL.FontSize = 20;
    % newPosition = [0.22 0.14 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1  );
    
    % % set all units inside figure to normalized so that everything is scaling accordingly
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % % do not show figure on screen
    % set(fh2, 'visible', 'on')
    % % set figure units to pixels & adjust figure size
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 800 700];
    % % define resolution figure to be saved in dpi
    % res = 500;
    % % recalculate figure size to be saved
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 5580 4320]/res;
    % % save figure
    % % print(fh2,'RequiredControl','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% ----------------- 3D plot of the converged solution --------------------
for ll = 1
    fh2 = figure;
    LaBeL3 = {'$\delta r_{1}$ [km]','$\delta r_{2}$ [km]','$\delta r_{3}$ [km]',...
        '$\delta v_{1}$ [km/s]','$\delta v_{2}$ [km/s]','$\delta v_{3}$ [km/s]'};
    
    
    Xrel21 = squeeze(sol(end,:,:));
    for i = 1:13
        Xrel2(:,i) = makima(x,Xrel21(:,i),Tnorm);
    end
    
    for i = 1:6
        j = i+7;
        subplot(2,3,i)
        plot(Time,Xrel2(:,j) - X_ode45(:,j),'LineWidth',2.5)
        ylabel(LaBeL3{i})
        xlabel('t [days]')
        grid on
        % set(gca,'FontSize',25)
    end
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 1000 600];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 3320]/res;
    % save figure
    % print(fh2,'PositionVelocityError_6DoF','-dpng',sprintf('-r%d',res))
    
end
%%
for ll = 1
    cMat1 = [c5;c6];
    cMat2 = [c3;c1];
    cMat0 = [c7;c4];
    
    fh2 = figure;
    hold on
    plt0  = zeros(Num_Agent);
    plt1   = zeros(Num_Agent);
    plt2  = zeros(Num_Agent);
    for jj = 1:Num_Agent
        idx = 1+N*(jj-1):N*jj;
        
        h0 = plot3(Xrel2(:,idx(8)),Xrel2(:,idx(9)),Xrel2(:,idx(10)),...
            '-','Color',c4,'LineWidth',2.5);
        
        plt0(:,jj) = plot3(X_ode45(:,idx(8)),X_ode45(:,idx(9)),X_ode45(:,idx(10)),...
            '-.','Color',cMat0(jj,:),'LineWidth',2.5);
        
        plt1(:,jj) = plot3(Xrel0(1,idx(1)),Xrel0(1,idx(2)),Xrel0(1,idx(3)),'*',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat1(jj,:),...
            'MarkerFaceColor',cMat1(jj,:)',...
            'MarkerSize',8);
        
        
        plt2(:,jj) = plot3(Xrel1(1,idx(1)),Xrel1(1,idx(2)),Xrel1(1,idx(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat2(jj,:),...
            'MarkerFaceColor',cMat2(jj,:),...
            'MarkerSize',8);
        
        h5 = plot3(0,0,0,'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','k',...
            'MarkerSize',15);
        plot3(Xrel0(:,idx(1)),Xrel0(:,idx(2)),Xrel0(:,idx(3)),'Color',cMat1(jj,:),'DisplayName','Initial Traj');
        plot3(Xrel1(:,idx(1)),Xrel1(:,idx(2)),Xrel1(:,idx(3)),'Color',cMat2(jj,:));
        
        grid on
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
        view(170,10)
    end
    
    hL = legend([h5, plt1(end,1), plt2(end,1), plt0(end,1), h0],...
        {'Chief Location','Agent1 Init Condition',...
        'Agent1 End Condition','ODE45', 'AGHF'},'AutoUpdate','off');
    % % hL = legend([h5, plt1(end,1), plt1(end,2), plt2(end,1), plt2(end,2), plt0(end,1), plt0(end,2)],...
    % %     {'Chief Location','Agent1 Init Condition','Agent2 Init Condition',...
    % %     'Agent1 End Condition','Agent2 End Condition','Agent1 Transfer Traj',...
    % %     'Agent2 Transfer Traj'},'AutoUpdate','off');
    % % % Programatically move the Legend
    % hL.FontSize = 25;
    newPosition = [0.68 0.2 0.1 0.1];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    set(gca,'FontSize',25)
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 900 800];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 4320]/res;
    % save figure
    % print(fh2,'Integrated_Trajectory_6DoF','-dpng',sprintf('-r%d',res))
    
end
% -------------------------------------------------------------------------

%% ------ Animation of the homotopy trasformation during convergence ------
for ll = 1
    % cMat1 = [c5;c6];
    % cMat2 = [c3;c1];
    % cMat0 = [c7;c4];
    %
    % fh2 = figure('Renderer', 'painters', 'Position', [10 10 1400 1000]);
    % view(68,10)
    % hold on
    % plt0  = zeros(Num_Agent);
    % plt1  = zeros(Num_Agent);
    % plt2  = zeros(Num_Agent);
    % for jj = 1:Num_Agent
    %     idx = 1+N*(jj-1):N*jj;
    %
    %     h6 = plot3(X_ode45(:,idx(8)),X_ode45(:,idx(9)),X_ode45(:,idx(10)),...
    %         '--','Color',cMat0(jj,:),'LineWidth',2.5);
    %
    %     plt1(:,jj) = plot3(Xrel0(1,idx(1)),Xrel0(1,idx(2)),Xrel0(1,idx(3)),'*',...
    %         'LineWidth',3,...
    %         'MarkerEdgeColor',cMat1(jj,:),...
    %         'MarkerFaceColor',cMat1(jj,:)',...
    %         'MarkerSize',8);
    %
    %
    %     plt2(:,jj) = plot3(Xrel1(1,idx(1)),Xrel1(1,idx(2)),Xrel1(1,idx(3)),'o',...
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
    %     h4 = plot3(Xrel0(:,idx(1)),Xrel0(:,idx(2)),Xrel0(:,idx(3)),'Color',cMat1(jj,:));
    %     h3 = plot3(Xrel1(:,idx(1)),Xrel1(:,idx(2)),Xrel1(:,idx(3)),'Color',cMat2(jj,:));
    %
    %     grid on
    %     xlabel('X [km]')
    %     ylabel('Y [km]')
    %     zlabel('Z [km]')
    %     % title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
    %
    %     grid on
    %     grid minor
    % end
    %
    % % r = 1;
    % % [x_sphere, y_sphere, z_sphere] = sphere(100);
    % % h = surf(r*x_sphere, r*y_sphere, r*z_sphere);
    % % set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
    % % axis equal
    %
    %
    %
    %
    % tpointsCorase = 250;
    % xpointsCorase = 250;
    % Tmesh = [0, logspace(-4,log10(tmax),tpointsCorase-1)]; % discretization of time
    % Xmesh = linspace(0,T,xpointsCorase); % discretization of the curve
    % [X,Y] = meshgrid(Tmesh,Xmesh);
    %
    % Agent1_u1 = interp2(t,x,sol(:,:,8),X,Y);
    % Agent1_u2 = interp2(t,x,sol(:,:,9),X,Y);
    % Agent1_u3 = interp2(t,x,sol(:,:,10),X,Y);
    %
    % Xrel2 = [Agent1_u1(1,:);Agent1_u2(1,:);Agent1_u3(1,:)];
    % X = Xrel2(1,:);
    % Y = Xrel2(2,:);
    % Z = Xrel2(3,:);
    % h2 = plot3(X,Y,Z,'Color',c9,'LineWidth',2);
    % hold on
    % h1 = plot3(X,Y,Z,'r','LineWidth',2);
    % grid on;
    % % title('Nonlinear Geomeric Planning - 3D configuration');
    % xlabel('$x[km]$');
    % ylabel('$y[km]$');
    % zlabel('$z[km]$');
    % hL = legend([h5, h4, h3, h6, h2, h1],...
    %     {'Chief','Deputy Int Condition',...
    %     'Deputy End Condition','Deputy Transfer Traj',...
    %     'Initial Guess','Homotopy Iteration'},...
    %     'AutoUpdate','off');
    % newPosition = [0.2 0.7 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    % set(gca,'FontSize',30)
    %
    % set(gca,'nextplot','replacechildren','visible','on')
    %
    % fig  = getframe(fh2);
    % [im,map] = rgb2ind(fig.cdata,256,'nodither'); % 65536
    % im(1,1,1,tpointsCorase) = 0;
    % for i = 1:tpointsCorase
    %     Xrel2 = [Agent1_u1(i,:);Agent1_u2(i,:);Agent1_u3(i,:)];
    %     X = Xrel2(1,:); h1.XDataSource='X';
    %     Y = Xrel2(2,:); h1.YDataSource='Y';
    %     Z = Xrel2(3,:); h1.ZDataSource='Z';
    %     refreshdata(h1,'caller');
    %     drawnow;
    %     view(68,10)
    %
    %     % view(22,7)
    %     % if i == tpointsCorase
    %     %     % set all units inside figure to normalized so that everything is scaling accordingly
    %     %     set(findall(fh2,'Units','pixels'),'Units','normalized');
    %     %     % do show figure on screen
    %     %     set(fh2, 'visible', 'on')
    %     %     % set figure units to pixels & adjust figure size
    %     %     fh2.Units = 'pixels';
    %     %     fh2.OuterPosition = [10 10 900 800];
    %     %     % define resolution figure to be saved in dpi
    %     %     res = 500;
    %     %     % recalculate figure size to be saved
    %     %     set(fh2,'PaperPositionMode','manual')
    %     %     fh2.PaperUnits = 'inches';
    %     %     fh2.PaperPosition = [0 0 5580 4320]/res;
    %     %     % save figure
    %     %     print(fh2,'AGHF_Solutio_6DoF','-dpng',sprintf('-r%d',res))
    %     % end
    %
    %
    %
    %     f = getframe(fh2);
    %     im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    % end
    % % imwrite(im,map,'AGHF_6DoF_SRP.gif','DelayTime',0.,'LoopCount',inf) %g443800

end
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% ######################## Required Functions #############################
% -------------------------------------------------------------------------
%% --------------------------- PDE for pdepe ------------------------------
function [c,f,s] = mypdexpde(x,t,u,DuDx,k1,k2,k3,N) % Define PDE; right-hand-side of AGHF

XChief = full([Chief_x(x);  Chief_y(x);  Chief_z(x);...
               Chief_dx(x); Chief_dy(x); Chief_dz(x)]);
LL = full(L(u,DuDx,k1,k2,k3,x,XChief));
CasADiresult = full(dLdx(u,DuDx,k1,k2,k3,x,XChief,LL))';
CasADir = [CasADiresult(1:13), CasADiresult(14:26)]; % [dL/dx, dL/dxdot]

f = CasADir(:,2);  % dL/dxdot
s = -CasADir(:,1); % -dL/dx
c = ones(N,1);

end
% -------------------------------------------------------------------------

%% --------------------- Initial condition for pdepe ----------------------
function u0 = mypdexic(x,X0,Xf,T)    % Initial condition for AGHF
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

%% --------------------- Boundary condition for pdepe ---------------------
function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t,X0,Xf,N) % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end
% -------------------------------------------------------------------------

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!   Nondimetionalized Dynamics   !!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [Xdot] = NonDimentionalized_FlatPlateOde(t,X,Target,Chasser)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

% Collecting relevant quantities
X_chief = full([Chief_x(t);  Chief_y(t);  Chief_z(t);...
                Chief_dx(t); Chief_dy(t); Chief_dz(t)]);
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

q_Chasser = X(1:4)/norm(X(1:4));
Omega_Chasser_body = X(5:7);
rho = X(8:10);  rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% F_body = F_CanonBall(ti,r_chasser,Chasser); % SRP force on the Chasser
F_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);

%***************************************%
% Rotational state differential equation
% thau_Chasser         = F_body(4:6);
I_Chasser            = diag(Chasser.Moment_Of_Inertia_Calculator());
Omega_matrix_Chasser = [0, -Omega_Chasser_body(1), -Omega_Chasser_body(2), -Omega_Chasser_body(3);
                        Omega_Chasser_body(1),  0,  Omega_Chasser_body(3), -Omega_Chasser_body(2);
                        Omega_Chasser_body(2), -Omega_Chasser_body(3),  0,  Omega_Chasser_body(1);
                        Omega_Chasser_body(3), Omega_Chasser_body(2), -Omega_Chasser_body(1), 0];

qdot      = 1/2*Omega_matrix_Chasser*q_Chasser;
Omega_dot = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body)); % + thau_Chasser);

%***************************************%
% Translational state differential equation
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\F_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;


% Integrating the relative trajectory
Rho_dot = [rho_prime; DelUsrp + del_ag + Coriolis_Effect];

% Collecting the rates of change
Xdot = [qdot; Omega_dot; Rho_dot]*Tc;

end
function [Xdot] = NonDimentionalized_FlatPlateControlledOde(t,X,Target,Chasser)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

u = [u11(t) u12(t) u13(t)]';
% Collecting relevant quantities
X_chief = full([Chief_x(t);  Chief_y(t);  Chief_z(t);...
                Chief_dx(t); Chief_dy(t); Chief_dz(t)]);
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

q_Chasser = X(1:4)/norm(X(1:4));
Omega_Chasser_body = X(5:7);
rho = X(8:10);  rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% F_body = F_CanonBall(ti,r_chasser,Chasser); % SRP force on the Chasser
F_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);

%***************************************%
% Rotational state differential equation
% thau_Chasser         = F_body(4:6);
I_Chasser            = diag(Chasser.Moment_Of_Inertia_Calculator());
Omega_matrix_Chasser = [0, -Omega_Chasser_body(1), -Omega_Chasser_body(2), -Omega_Chasser_body(3);
                        Omega_Chasser_body(1),  0,  Omega_Chasser_body(3), -Omega_Chasser_body(2);
                        Omega_Chasser_body(2), -Omega_Chasser_body(3),  0,  Omega_Chasser_body(1);
                        Omega_Chasser_body(3), Omega_Chasser_body(2), -Omega_Chasser_body(1), 0];

qdot      = 1/2*Omega_matrix_Chasser*q_Chasser;
Omega_dot = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body) + u); % + thau_Chasser);

%***************************************%
% Translational state differential equation
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\F_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;

% Integrating the relative trajectory
Rho_dot = [rho_prime; DelUsrp + del_ag + Coriolis_Effect];

% Collecting the rates of change
Xdot = [qdot; Omega_dot; Rho_dot]*Tc;

end
function [Xdot] = NonDimentionalized_NLode1(t,X,Target)
global mu_Earth Tc 
format long
ti = t*Tc;
X_chief = full([Chief_x0(t);  Chief_y0(t);  Chief_z0(t);...
                Chief_dx0(t); Chief_dy0(t); Chief_dz0(t)]);
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
r_target = X_chief(1:3);  v_target = X_chief(4:6); rt_norm = (r_target.'*r_target).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho =  rho_bar; rho_prime = drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(r_target)*v_target; h_norm = norm(h_vec);  u_acc = U_eci_Target;
eh = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
u_acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot =  Xdot*Tc;

end
function [Xdot] = NonDimentionalized_NLode2(t,X,Target)
global mu_Earth Tc 
format long
ti = t*Tc;
X_chief = full([Chief_xf(t);  Chief_yf(t);  Chief_zf(t);...
                Chief_dxf(t); Chief_dyf(t); Chief_dzf(t)]);
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
r_target = X_chief(1:3);  v_target = X_chief(4:6); rt_norm = (r_target.'*r_target).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = rho_bar; rho_prime = drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(r_target)*v_target; h_norm = norm(h_vec);  u_acc = U_eci_Target;
eh = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
u_acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot = Xdot*Tc;

end
% -------------------------------------------------------------------------

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!   CasADi Functions   !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [Xdot] = NonDimentionalized_FlatPlateOdeCasADi(t,X,Target,Chasser,Xchief)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

% Collecting relevant quantities
r_target = Xchief(1:3); v_target = Xchief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

q_Chasser = X(1:4)/norm(X(1:4));
Omega_Chasser_body = X(5:7);
rho = X(8:10);  rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,Xchief,Target); % SRP force on the Target
% F_body = F_CanonBall(ti,r_chasser,Chasser); % SRP force on the Chasser
F_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);

%***************************************%
% Rotational state differential equation
% thau_Chasser         = F_body(4:6);
I_Chasser            = diag(Chasser.Moment_Of_Inertia_Calculator());
Omega_matrix_Chasser = [0, -Omega_Chasser_body(1), -Omega_Chasser_body(2), -Omega_Chasser_body(3);
                        Omega_Chasser_body(1),  0,  Omega_Chasser_body(3), -Omega_Chasser_body(2);
                        Omega_Chasser_body(2), -Omega_Chasser_body(3),  0,  Omega_Chasser_body(1);
                        Omega_Chasser_body(3), Omega_Chasser_body(2), -Omega_Chasser_body(1), 0];

qdot      = 1/2*Omega_matrix_Chasser*q_Chasser;
Omega_dot = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body)); % + thau_Chasser);

%***************************************%
% Translational state differential equation
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\F_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;


% Integrating the relative trajectory
Rho_dot = [rho_prime; DelUsrp + del_ag + Coriolis_Effect];

% Collecting the rates of change
Xdot = [qdot; Omega_dot; Rho_dot]*Tc;

end
% -------------------------------------------------------------------------

function [Xdot, DelUsrp, U_eci_Chasser, U_eci_Target,Costheta,del_ag] = ControlsAccelerationCheck(t,X,Target,Chasser,q_Chasser,X_chief)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

q_Chasser = q_Chasser/norm(q_Chasser);
% Collecting relevant quantities
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

rho = X(1:3);  rho_prime = X(4:6);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% F_body = F_CanonBall(ti,r_chasser,Chasser); % SRP force on the Chasser
[F_body,Costheta] = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);

%***************************************%
% Translational state differential equation
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\F_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;


% Integrating the relative trajectory
Rho_dot = [rho_prime; DelUsrp + del_ag + Coriolis_Effect];

% Collecting the rates of change
Xdot = Rho_dot*Tc;

end


