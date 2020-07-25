function [B] = BMatrix(COE)
%--------------------------------------------------------------------------
% [B] = BMatrx(OE, nu, f)
%
% This function computes Gauss' Variational Equation matrix
%
% [inputs]: - [OE] : Keplerian (osculating) orbital elements
%
% [outputs]: -[B] : Gauss' Variational Equation matrix
%--------------------------------------------------------------------------
format long
global mu_Earth

a = COE(1); e = COE(2); i = COE(3); omega = COE(5); f = COE(6);

B(1,1) = 2.*a.^2.*e.*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*sin(f);
B(1,2) = 2.*a.^2.*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*(1+e.*cos(f));
B(1,3) = 0;

B(2,1) = a.*(1+(-1).*e.^2).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*sin(f);
B(2,2) = (a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*(a.*e.*(1+(-1).*e.^2).*(1+e.*cos( ...
  f)).^(-1)+cos(f).*(a.*(1+(-1).*e.^2)+a.*(1+(-1).*e.^2).*(1+e.*cos( ...
  f)).^(-1)));
B(2,3) = 0;

B(3,1) = 0;
B(3,2) = 0;
B(3,3) = a.*(1+(-1).*e.^2).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*(1+e.*cos(f)) ...
  .^(-1).*cos(f+omega);

B(4,1) = 0;
B(4,2) = 0;
B(4,3) = a.*(1+(-1).*e.^2).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*(1+e.*cos(f)) ...
  .^(-1).*csc(i).*sin(f+omega);
B(5,1) = (-1).*a.*e.^(-1).*(1+(-1).*e.^2).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2) ...
  .*cos(f);

B(5,2) = e.^(-1).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*(a.*(1+(-1).*e.^2)+a.*( ...
  1+(-1).*e.^2).*(1+e.*cos(f)).^(-1)).*sin(f);
B(5,3) = (-1).*a.*(1+(-1).*e.^2).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*(1+e.* ...
  cos(f)).^(-1).*cot(i).*sin(f+omega);

B(6,1) = a.*e.^(-1).*(1+(-1).*e.^2).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*cos(f);
B(6,2) = e.^(-1).*(a.*(1+(-1).*e.^2).*mu_Earth).^(-1/2).*((-1).*a.*(1+(-1).*e.^2) ...
  +(-1).*a.*(1+(-1).*e.^2).*(1+e.*cos(f)).^(-1)).*sin(f);
B(6,3) = 0;

end