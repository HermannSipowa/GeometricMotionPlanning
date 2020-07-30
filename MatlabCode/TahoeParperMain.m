%% Pendulum Controlled Swing (with gravity and obstacle)
% Procedure:
% 1. Passing initial condition and boundary condition to the PDE
% 2. Solve PDE
% 3. The solution at tmax approximates the steady state solution
% 4. Extract the controls from sol(tpoints,:,:)
% 5. Simulate the motion of unicycle using extracted controls
%--------------------------------------------------------------------------
global Xrel0 Xrelf Penalty ae mu_Earth ...
       JD AU Target Chasser X_chief ...
       IntTime tspan
tmax     = 10;      % Integration time for PDE
tpoints  = 2000;   % Number of points in the time discretization
h        = 0.01;   % Step-size of integration[s]
m        = 0;      % Integrator parameter (due to cartesian coordinates)\
space    = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of the curve
Penalty  = 1E4;    % Penalty for moving in unfeasible/constrained directions
ae       = 6378.136;
mu_Earth = 3.986004415E5;
JD       = 2456296.25;
AU       = 149597870.7;
Chasser  = Spacecraft([30 30 4 116 40 0.9 0.9 3000 6000]); % Describing the Deputy sailcraft characteristic (refer Spacecraft class definition)
Target   = Spacecraft([10 6 2 40 10 0.2 0.5 6 12]);         % Describing the Chief sailcraft characteristic (refer Spacecraft class definition)

%% Target initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a1      = 1.5E+4; % Semi-major axis in Km
e1      = 0.0;    % Eccentricity
inc1    = 45;     % Inclination in deg0
BigOmg1 = 20;     % RAAN in deg
LitOmg1 = 90;     % AOP in deg
M1      = 0;      % Mean anomaly in deg

% Setting the integrator parameters
Period  = 2*pi*sqrt(a1^3/mu_Earth);
IntTime = Period;
tspan   = linspace(0,IntTime,10000);
options = odeset('RelTol',2.22045e-13,'AbsTol',2.22045e-30);

% Constructing the chief inital conditions from the chief's OEs
inc1 = deg2rad(inc1); BigOmg1 = deg2rad(BigOmg1); LitOmg1 = deg2rad(LitOmg1); M1 = deg2rad(M1);
COE_c = [a1,e1,inc1,BigOmg1,LitOmg1,M1];
[Position_target,Velocity_target]  = COEstoRV_MeanAnonaly(COE_c,mu_Earth,1);
X0_chief = [Position_target; Velocity_target];

% integrting the chief/reference trajectory
[~, X_chief] = ode113(@(t,X)ChiefMotionODE(t,X,Target),tspan,X0_chief,options);

TN = DCM(Position_target,Velocity_target); rt_norm = norm(Position_target);
h_vec = cross(Position_target,Velocity_target); h_norm = norm(h_vec);
eh = h_vec/h_norm;
U_eci_Target  = F_CanonBall(tspan(1),Position_target,Target); % SRP force on the Target
N_nudot = h_vec/rt_norm^2 + dot(U_eci_Target,eh)*Position_target/h_norm;

%% Chaser initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a2      = a1+10;        % Semi-major axis in Km
e2      = e1+1e-3;      % Eccentricity
inc2    = inc1 + 0.001; % Inclination in deg
BigOmg2 = BigOmg1;      % RAAN in deg
LitOmg2 = LitOmg1;      % AOP in deg
M2      = M1;           % True anomaly in rad

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

%% Desired trajectory initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a3      = a1;           % Semi-major axis in Km
e3      = e1 + .4/a1;   % Eccentricity
inc3    = inc1;         % Inclination in deg
BigOmg3 = BigOmg1;      % RAAN in deg
LitOmg3 = LitOmg1;      % AOP in deg
M3      = M1;           % True anomaly in rad

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
sol = pdepe(m,@TahoePDE,@TahoeIC,@TahoeBC,tspan,space);

%%  Extract Controls
n = length(Aug_chasser_Xo);
m = length(tspan);
Tau_ctrl  = zeros(2,m);
[~, Xdot] = pdeval(m,tspan,sol(end,:,n),tspan);
for i = 1:m
    Iinv = inv(diag(Chasser.Moment_Of_Inertia_Calculator()));
    Fbar = eye(n);
    Fbar(:,1) = [u(1:4);zeros(9,1)];
    Fbar(:,5:7) = [zeros(10,3);eye(3)];
    Fbar(:,11:end) = [zeros(4,3);Iinv;zeros(6,3)];
    X_aug = [X_chief; sol(end,i,:)];
    fd    = RelativeMotionODE(i,X_aug,Target,Chasser);
    Tau_ctrl(:,i) = [zeros(3,10) eye(3)]/Fbar*(Xdot(i,:)-fd);
end

%% Find x^* by integrating ODE
Statestar = zeros(m,n);
tol = 1E-11;
Statestar(1,:) = X0_chasser;
for ii = 1:(m-1)
    [Statestar_f, out, h_next] = Runge_Kutta_Fehlberg_4_5(...
        @(t,Aug_X)my2Bodyfunc(t,Aug_X,Tau_ctrl,Target,Chasser),...
        Statestar(ii,:)',tspan(ii),h,tspan(ii+1),tol);
    h = h_next;
    Statestar(ii+1,:) = Statestar_f;
end

%% Saving results
%-------------------------------------------
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
Xstar = interp1(x,Statestar',X)';
save([dataPath,'Xstar','.mat'], 'Xstar','-v7.3')
save([dataPath,data,'.mat'], 'u1','-v7.3')



