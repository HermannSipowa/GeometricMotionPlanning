function [F_body] = F_SPR_FlatPlate(time,X,Spacecraft,q)
% -------------------------------------------------------------------------
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

% Description of the sun location relative to the sailcraft body frame
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ECI_to_BodyFrame = Quaternion_to_DCM(q);
Pcrp = 1357/(299792458); % Solar pressure

JulianDay = JD+time/86400;
[XS, ~, ~] = Ephem(JulianDay,3,'EME2000');
s = (XS-X(1:3,:))/norm((XS-X(1:3,:))); s = ECI_to_BodyFrame*s;

% Collecting sailcraft characteristics
rhos = Spacecraft.rhos;rhod = Spacecraft.rhod;
area = Spacecraft.L*Spacecraft.l;
Cr   = Spacecraft.Cr;
mass = Spacecraft.M+Spacecraft.m;
rs_squared = norm((XS-X(1:3,:))/AU)^2;
P = -Cr*Pcrp*area/(mass*rs_squared);
n = [0 0 1]'; Costheta = dot(n,s);

% Computing the external force exerted on the sailcraft (in km/s^2)
if Costheta<0
    a_srp = zeros(3,1);
else
    a_srp = -P*dot(n,s)*((1-rhos)*s+(rhod+2*rhos*dot(n,s))*n)/1000;
end
F = mass*a_srp;
% Computing the torque acting on the sailcraft
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rb     = Spacecraft.rhos*[0 0 1]';   % lever arm wrt the spacecraft center of mass 
thau_B = cross(rb,F);                % the solar radiation torque

F_body = [a_srp; thau_B];
end