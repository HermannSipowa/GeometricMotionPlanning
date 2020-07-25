function [acc,Xrel,Vs,beta] = F_CanonBall(t,X,Spacecraft)
%--------------------------------------------------------------------------
% [F] = F_CanonBall(t,X,Spacecraft)
%
% The SRP force model comes from the following paper:
% Farres, Ariadna, Cassandra Webster, and David Folta. 
% "High fidelity modeling of SRP and its effect on the relative motion
%  of Starshade and WFIRST." In 2018 Space Flight Mechanics Meeting,
%  p. 2227. 2018.
%
% [inputs]: - [X]          : 6 by 1 vector of cartesian states (position and velocity)
%           - [t]          : simulation time (for computing the sun sun ephemeris)
%           - [Spacecraft] : class defining the spacecraft's charateristics
%
% [outputs]: - [F] : 3 by 1 vector representing the CanonBall SRP force
%--------------------------------------------------------------------------
format long
global JD

Cr   = Spacecraft.Cr;
area = Spacecraft.L*Spacecraft.l;
mass = Spacecraft.M+Spacecraft.m;

Pcrp = 1357/(299792458); % Solar pressure
AU = 149597870.7; % Distance between Earth and Sun (in Km)
gamma = (Cr*area)/mass;

JulianDay = JD+t/86400;
[XS, Vs, ~] = Ephem(JulianDay,3,'EME2000');
Vs = -Vs;
XS = -XS;
Xrel = XS - X(1:3,:);

beta = -Pcrp * AU^2 * gamma /1000;

acc =  beta * Xrel/ norm(Xrel)^3;

end