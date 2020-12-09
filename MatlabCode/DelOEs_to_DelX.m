function [DelX] = DelOEs_to_DelX(DelOE,OE_Chief,mu_Earth)


a_chief      = OE_Chief(1); % Semi-major axis in Km 
e_chief      = OE_Chief(2); % Eccentricity
inc_chief    = OE_Chief(3); % Inclination in deg
BigOmg_chief = OE_Chief(4); % RAAN in deg
LitOmg_chief = OE_Chief(5); % AOP in deg
M_chief      = OE_Chief(6); % True anomaly in deg

% Computing the nearly non-singular orbital elements
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
LitOmg_chief = deg2rad(LitOmg_chief); M_chief = deg2rad(M_chief);
q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
f_chief = M_chief; theta_chief = f_chief+LitOmg_chief;

OE_chief = [a_chief, theta_chief, inc_chief, q1_chief, q2_chief, BigOmg_chief];
AMap = ForwardMapping(OE_chief, mu_Earth); % Linear mapping matrix


%% Specify the initial relative-orbital element of the deputies
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dela = DelOE(1);
dele = DelOE(2);
deli = DelOE(3);
delLitOmg = DelOE(4);
delBigOmg = DelOE(5);
delM = DelOE(6);

delq1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
delq2 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
deltheta = delLitOmg + delM; % Relative true latitude in rad
delCOE = [dela, deltheta, deli, delq1, delq2, delBigOmg]';

DelX = AMap*delCOE;


end