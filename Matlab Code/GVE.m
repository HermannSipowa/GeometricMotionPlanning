function [OEdot] = GVE(t,COE,Spacecraft)
%--------------------------------------------------------------------------
% [OEdot] = GVE(t,OE,nu,u)
%
% This function impliments Gauss' Variational Equation
%
% [inputs]: - [t]  : time vector (integration purposes)
%           - [OE] : Keplerian (osculating) orbital elements
%           - [nu] : Body's gravitational constant
%           - [u]  : Pertubating acceleration (expressed in LVLH frame)
%           - [f]  : true anomaly
%
% [outputs]: -[OEdot] : rate of change of Keplerian (osculating) orbital elements (last element is mean anomaly)
%--------------------------------------------------------------------------
format long
global mu_Earth

%% Computing the acceleration acting on the spacecraft
[r, ~]  = COEstoRV(COE,mu_Earth);
F_body = 0*F_CanonBall(t,r,Spacecraft); % F_SPR_FlatPlate(t, r_chasser, Chasser, q_Chasser); % 

[ON] = OE_DCM(COE);
u = ON*F_body;

%% Computing rate of change of the true anomaly
a = COE(1); e = COE(2); f = COE(6);
h = sqrt(mu_Earth*a*(1-e^2)); 
r = (a*(1-e^2))/(1+e*cos(f));
nudot = h/r^2;

%% Gauss' Variational Equation
B = BMatrix(COE);
OEdot = B*u;
OEdot(end) = OEdot(end)+nudot;

end