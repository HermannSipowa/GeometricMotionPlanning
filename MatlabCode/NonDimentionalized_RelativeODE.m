function [Xdot] = NonDimentionalized_RelativeODE(t, Xaug, mu, Rc, Rrho, Tc)
format long

%% Collecting relevant quantities
X_chief = Xaug(1:6,:); 
X_rho  = Xaug(7:12,:);
rho_bar = X_rho(1:3,:); drho_bar = X_rho(4:6,:); 
r_brar = X_chief(1:3,:); v_bar = X_chief(4:6,:); 

%% Redimensionalizing the variables in the ODE
r_target = Rc*r_brar; v_target = 1/Tc*Rc*v_bar; rt_norm = norm(r_target);
r_c = [rt_norm; 0; 0] + rho_bar; rc_norm = norm(r_c);
rho = Rrho*rho_bar; rho_prime = 1/Tc*Rrho*drho_bar;

%% Calculating the respective accelerations
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = cross(r_target,v_target);
Omega = TN*( h_vec/rt_norm^2);  
Omegadot = TN*( -2*dot(r_target,v_target)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu/rc_norm^3*r_c + mu/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*cross(Omega,rho_prime) ...
                  - cross(Omega,cross(Omega,rho)) ...
                  - cross(Omegadot,rho);

%% ODE of the chief motion
Xdot(1:3,:) = v_bar;
Xdot(4:6,:) = Rc\(-mu*r_target/rt_norm^3)*Tc^2;
%% Nondimetional relative ODE of the deputy
Xdot(7:9,:) = drho_bar;
Xdot(10:12,:) = Rrho\(del_ag + Coriolis_Effect)*Tc^2;

end