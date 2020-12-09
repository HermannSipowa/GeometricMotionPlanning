function [Xdot] = NonDimentionalized_2BP(t, X, mu, R, Tc)
format long
X = X(1:6,:); 
r_brar = X(1:3,:); v_bar = X(4:6,:); 
Xdot(1:3,:) = v_bar;
Xdot(4:6,:) = -mu*r_brar*Tc^2/(r_brar'*(R'*R)*r_brar)^(3/2);
end