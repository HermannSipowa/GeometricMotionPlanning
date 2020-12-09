clear
close all
clc
start_up
format long e
global JD AU mu_Earth tspan Xnom_chief Target Chasser t
Target  = Spacecraft([10 10 1 40 10 0.2 0.5 6 12]); 
Chasser = Spacecraft([30 30 4 116 40 0.9 0.9 3000 6000]);
JD = 2456296.25;
AU = 149597870.7;
mu_Earth     = 3.986004415E5; % Earth gravitational parameter

%% Initializing the chief's state
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a_chief      = 41798; % Semi-major axis in Km
e_chief      = 0.75;  % Eccentricity
inc_chief    = 3.22;  % Inclination in deg0
BigOmg_chief = 65;    % RAAN in deg
LitOmg_chief = 180;   % AOP in deg
M_chief      = 0;     % True anomaly in deg
OE_Chief = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief];
OE_Chief(3:end) = deg2rad(OE_Chief(3:end));
[Position_target,Velocity_target]  = COEstoRV(OE_Chief,mu_Earth);
X0_Chief = [Position_target; Velocity_target];

%% Initializing the deputy's state
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dela = 0;
dele = 1e-8;
deli = 1e-7;
delLitOmg = 0;
delBigOmg = 0;
delM = -5e-7;
for i = 1
    DelOE = [dela(i); dele(i); deli(i);...
        delLitOmg(i); delBigOmg(i); delM(i)];
    Xrel = DelOEs_to_DelX(DelOE,OE_Chief,mu_Earth);
end
sigma = [.1,.5,.3]';
Omega_chaser = deg2rad([-.5, .4, .1])'; % Angular velocity (inertial frame)


%% Integrating the chief' norminal trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Period  = 2*pi*sqrt(a_chief^3/mu_Earth);
IntTime = 1*Period;
tspan   = linspace(0,IntTime,1e5);
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-30);
mu_dot  = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
[~, Xnom_chief] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),tspan,X0_Chief,options);



%% Co,puting the Lie brackets;
t = 0;
X0 = [sigma; Omega_chaser; Xrel]; X0_array = dlarray(X0);
[f,grad] = dlfeval(@MRP_RelativeMotionODE,X0_array);


%% Ordinary Differential Equation representing the dynamics of the system
function [f,grad] = MRP_RelativeMotionODE(Xaug)
global mu_Earth Xnom_chief tspan Target Chasser t
idx = find(tspan <= t, 1, 'last' );

tilde = @(v) [0    -v(3)  v(2);
              v(3)   0   -v(1);
             -v(2)  v(1)    0];

%***************************************%
% Collecting relevant quantities
Sigma = Xaug(1:3,:); sigma = (Sigma.'*Sigma).^(1/2);
Omega_Chasser_body = Xaug(4:6,:);
rho_aug  = Xaug(7:12,:);

r_target = Xnom_chief(idx,1:3)'; v_target = Xnom_chief(idx,4:6)';
rt_norm = norm(r_target);
rho = rho_aug(1:3,:); rho_prime = rho_aug(4:6,:);
r_c = [rt_norm+rho(1) rho(2) rho(3)]'; rc_norm = (r_c.'*r_c).^(1/2);  % RIC position of the chasser
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame
r_chasser = TN.'*r_c; % ECI position of the chasser

%***************************************%
% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, V_sun, beta_chief] = F_CanonBall(t,r_target,Target); % SRP force on the Target
% U_eci_Chasser = F_CanonBall(t,r_chasser,Chasser); % SRP force on the Chasser

F_body = fSPR_FlatPlate_MRP(t, r_chasser, Chasser, Sigma); % F_CanonBall(t,r_chasser,Chasser); %  

%***************************************%
% Rotational state differential equation
thau_Chasser     = F_body(4:6,:);
I_Chasser        = diag(Chasser.Moment_Of_Inertia_Calculator());

Sigmadot = 1/4*((1-sigma^2)*eye(3)+2*tilde(Sigma)+2*(Sigma*Sigma'))*Omega_Chasser_body;
Omega_dot = inv(I_Chasser)*(-cross(Omega_Chasser_body,I_Chasser*Omega_Chasser_body) + thau_Chasser);

%***************************************%
% Translational state differential equation
CN = MRP_to_DCM(Sigma);
U_eci_Chasser = CN.'*F_body(1:3,1);

%***************************************%
% Computing the angular velocity and angular acceleration of the target frame
h_vec = cross(r_target,v_target); h_norm = norm(h_vec);  u_acc = U_eci_Target;
eh = h_vec/h_norm; h_dot = cross(r_target,u_acc); r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
V_rel = V_sun - v_target;
u_acc_dot = beta_chief*(V_rel/r_rel_norm^3 - 3*dot(Xrel,V_rel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

%***************************************%
% Computing the relative acceleration
delf     =  -mu_Earth/rc_norm^3*r_c + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelU_srp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP

%***************************************%
% Integrating the relative trajectory
Rho_dot(1:3,:) = rho_prime;
Rho_dot(4:6,:) = DelU_srp + delf - 2*tilde(Omega)*rho_prime ...
                - tilde(Omegadot)*rho - tilde(Omega)*tilde(Omega)*rho;

            
%***************************************%
% Collecting the rates of change

f = [Sigmadot; Omega_dot; Rho_dot];

n = length(Xaug);
for i = 1:n
    grad(i,:) = dlgradient(f(i),Xaug)';
end

end

function [F_body] = fSPR_FlatPlate_MRP(time,X,Spacecraft,sigma)
% ------------------------------
% [F_body] = F_SPR_FlatPlate(t, X, Spacecraft, Aug_X)
% 
% The SRP force model comes from the following paper:
% J. Van Der Ha and V. Lappas, "Long-Term Attitude Drift of Spinning 
% Spacecraft Under Solar Radiation Torques", Journal of Guidance, Control,
% and Dynamics, vol. 30, no. 5, pp. 1470-1479, 2007.
%
% [inputs]: - [time]       : simulation time (for computing the sun sun ephemeris)
%           - [X]          : Spacecraft inertial state
%           - [Spacecraft] : class defining a solar sail
%
% [outputs]: - [F_body] : F(1:3,1) = Acceleration on acting on the flat plate,
%                         F(4:6,1) = Torque acting on the flat plate 
%                         expressed in the body fixed frame of the sailcraft
% -------------------------------------------------------------------------
format long
global JD AU
tilde = @(v) [0    -v(3)  v(2);
              v(3)   0   -v(1);
             -v(2)  v(1)    0];
% Description of the sun location relative to the sailcraft body frame
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ECI_to_BodyFrame = MRP_to_DCM(sigma);
Pcrp = 1357/(299792458); % Solar pressure

JulianDay = JD+time/86400;
[XS, ~, ~] = Ephem(JulianDay,3,'EME2000');
DelR = XS-X(1:3,:); NomrDelR = (DelR.'*DelR).^(1/2);
s = DelR/NomrDelR; s = ECI_to_BodyFrame*s;

% Collecting sailcraft characteristics
rhos = Spacecraft.rhos;rhod = Spacecraft.rhod;
area = Spacecraft.L*Spacecraft.l;
Cr   = Spacecraft.Cr;
mass = Spacecraft.M+Spacecraft.m;
rs_squared = (DelR.'*DelR)/AU^2; 
P = -Cr*Pcrp*area/(mass*rs_squared);
n = [0 0 1]'; Costheta = n.'*s;

% Computing the external force exerted on the sailcraft (in km/s^2)
if Costheta<0
    a_srp = zeros(3,1);
else
    a_srp = -P*Costheta*((1-rhos)*s+(rhod+2*rhos*Costheta)*n)/1000;
end
F = mass*a_srp;
% Computing the torque acting on the sailcraft
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rb     = Spacecraft.rhos*[0 0 1]';   % lever arm wrt the spacecraft center of mass 
thau_B = tilde(rb)*F;                % the solar radiation torque

F_body = [a_srp; thau_B];
end

function [BN] = MRP_to_DCM(x)

tilde = @(v) [0    -v(3)  v(2);
              v(3)   0   -v(1);
             -v(2)  v(1)    0];
BN = eye(3)+(8*tilde(x)^2-4*(1-x'*x)*tilde(x))/(1+x'*x)^2;

end


















