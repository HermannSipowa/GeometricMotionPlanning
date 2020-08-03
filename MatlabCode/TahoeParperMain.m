%% Pendulum Controlled Swing (with gravity and obstacle)
% Procedure:
% 1. Passing initial condition and boundary condition to the PDE
% 2. Solve PDE
% 3. The solution at tmax approximates the steady state solution
% 4. Extract the controls from sol(tpoints,:,:)
% 5. Simulate the motion of unicycle using extracted controls
%--------------------------------------------------------------------------
global Xrel0 Xrelf Penalty ae mu_Earth JD AU...
    Target Chasser X_chief IntTime tspan Progress
Progress = 0;
ae       = 6378.136;
mu_Earth = 3.986004415E5;
JD       = 2456296.25;
AU       = 149597870.7;
Chasser  = Spacecraft([100 50 4 116 40 0.9 0.9 3000 6000]); % Deputy characteristic
Target   = Spacecraft([0 0 2 40 10 0.2 0.5 6 12]); % Target characteristic
Penalty  = 1E5;  % Penalty for moving in unfeasible/constrained directions
m        = 0;    % Integrator parameter (due to cartesian coordinates)
S_Depth  = 1E3;  % Integration time for PDE
S_Int    = 1E4;  % Number of points in the time discretization
Space    = linspace(0,S_Depth,S_Int); % Discretization of the curve

%% Inital condition in Orbital Elements

% Target initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a1      = 1.5E+4; % Semi-major axis in Km
e1      = 0.1;    % Eccentricity
inc1    = 45;     % Inclination in deg0
BigOmg1 = 20;     % RAAN in deg
LitOmg1 = 45;     % AOP in deg
M1      = 0;      % Mean anomaly in deg

% Chaser initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a2      = a1+5;       % Semi-major axis in Km
e2      = e1+1e-3;    % Eccentricity
inc2    = inc1+2e-3;  % Inclination in deg
BigOmg2 = BigOmg1+.1; % RAAN in deg
LitOmg2 = LitOmg1+.1; % AOP in deg
M2      = M1;         % True anomaly in rad

% Desired trajectory initial conditions(Orbital Elements)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a3      = a1;        % Semi-major axis in Km
e3      = e1+0.4/a1; % Eccentricity
inc3    = inc1;      % Inclination in deg
BigOmg3 = BigOmg1;   % RAAN in deg
LitOmg3 = LitOmg1;   % AOP in deg
M3      = M1;        % True anomaly in rad


%% Inital condition in coratesian corditnates

% integrting the chief/reference trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Period  = 2*pi*sqrt(a1^3/mu_Earth);
IntTime = 1/2*Period;
tspan   = linspace(0,IntTime,1E4);
options = odeset('RelTol',2.22045e-13,'AbsTol',2.22045e-30);

% Constructing the chief inital conditions from the chief's OEs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
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
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
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
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
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
sol = pdepe(m,@TahoePDE,@TahoeIC,@TahoeBC,tspan,Space);
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
[~, Xstar_f] = ode113(@(t,Aug_X)ControlledRelativeMotionODE(t,Aug_X,Tau_ctrl,Target,Chasser),tspan,Xrel0);
toc

%% Computing the uncontrolled trajectory
tic
[~, XUncontrolled] = ode113(@(t,Aug_X)RelativeMotionODE(t,Aug_X,Target,Chasser),tspan,Xrel0);
toc

%% Saving results
%-------------------------------------------
dataPath = 'Outputs_Data/';
tpointsCorase = 1000;
xpointsCorase = 1000;
Tmesh = [0 logspace(-4,log10(S_Depth),tpointsCorase-1)]; % discretization of time
Xmesh = linspace(0,IntTime,xpointsCorase); % discretization of the curve
[X,Y] = meshgrid(Tmesh,Xmesh);
for i = 1:length(Xrel0)
    data = sprintf('State%d', i);
    udata = interp2(Space,tspan,sol(:,:,i),X,Y);
    save([dataPath,data,'.mat'], 'udata','-v7.3')
end

Xstar = interp1(tspan,Xstar_f,Xmesh);
save([dataPath,'Xstar','.mat'], 'Xstar','-v7.3')

XUnctrl = interp1(tspan,XUncontrolled,Xmesh);
save([dataPath,'XUnctrl','.mat'], 'XUnctrl','-v7.3')

Ctrl = interp1(tspan,Tau_ctrl',Xmesh)';
save([dataPath,'Ctrl','.mat'], 'Ctrl','-v7.3')



