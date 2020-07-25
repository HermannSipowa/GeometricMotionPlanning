function [ON] = OE_DCM(COE)
%--------------------------------------------------------------------------
% [ON] = OE_DCM(COE)
%
% This function computes the direct cosine matrix between the inertial
% frame to the LVLH frame
% 
% [inputs]:  - [COE] : Keplerian Orbital elements (COE = [a e inc BigOmg LitOmg nu])
%
% [outputs]: - [ON]  : DCM from ECI to LVLH
%--------------------------------------------------------------------------
inc = COE(3); W = COE(4); theta = COE(5)+COE(6);

R1=[cos(W) sin(W) 0; -sin(W) cos(W) 0; 0 0 1];
R2=[1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)];
R3=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];

ON = R3*R2*R1;

end