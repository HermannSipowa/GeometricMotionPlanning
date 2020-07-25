function Xaug_dot = UnperturbedRelativeMotionODE(t, Xaug,mu_Earth)
format long
% global mu_Earth

%% Collecting relevant quantities
Target_States = Xaug(1:6,:); 
rho_aug  = Xaug(7:12,:);

r_target = Target_States(1:3,:); v_target = Target_States(4:6,:);
rt_norm = norm(r_target);
rho = rho_aug(1:3,:); rho_prime = rho_aug(4:6,:);
r_c = [rt_norm+rho(1) rho(2) rho(3)]'; rc_norm = norm(r_c); % RIC position of the chasser
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

%% Computing the angular velocity and angular acceleration of the target frame
h_vec = cross(r_target,v_target);
Omega = TN*( h_vec/rt_norm^2);  
Omegadot = TN*( -2*dot(r_target,v_target)*h_vec/rt_norm^4 );

%% Computing the relative acceleration
delf     =  -mu_Earth/rc_norm^3*r_c + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal

%% Integrating the the target trajectory
Xdot(1:3,:) = v_target;
Xdot(4:6,:) = -mu_Earth/rt_norm^3*r_target;

%% Integrating the relative trajectory
Rho_dot(1:3,:) = rho_prime;
Rho_dot(4:6,:) =  delf - 2*cross(Omega,rho_prime) ...
                - cross(Omegadot,rho) - cross(Omega,cross(Omega,rho));

%% Collecting the rates of change
Xaug_dot = [Xdot; Rho_dot];

end