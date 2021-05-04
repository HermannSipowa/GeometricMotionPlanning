function Sigmadot=myMRPsODE(t,Sigma,Omega)

sigma=norm(Sigma);
sigma_tilde=[0 -Sigma(3) Sigma(2); Sigma(3) 0 -Sigma(1); -Sigma(2) Sigma(1) 0];
Sigmadot=1/4*((1-sigma^2)*eye(3)+2*sigma_tilde+2*(Sigma*Sigma'))*Omega;
end