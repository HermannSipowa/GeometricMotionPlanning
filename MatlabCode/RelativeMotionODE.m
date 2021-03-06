function Xaug_dot = RelativeMotionODE(t, Xaug, Target, Chasser)
format long
global mu_Earth Tc
ti = t/Tc;
X_chief = full([Chief_x(t);  Chief_y(t);  Chief_z(t);...
    Chief_dx(t); Chief_dy(t); Chief_dz(t)]);
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

%% Collecting relevant quantities
q_Chasser = Xaug(1:4,1)/norm(Xaug(1:4,:));
Omega_Chasser_body = Xaug(5:7,:);
rho_aug  = Xaug(8:13,:);

r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
rho = rho_aug(1:3,:); rho_prime = rho_aug(4:6,:);
r_c = [rt_norm+rho(1) rho(2) rho(3)]'; rc_norm = (r_c.'*r_c).^(1/2);  % RIC position of the chasser
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame
r_chasser = TN\r_c; % ECI position of the chasser


%% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, V_sun, beta_chief] = F_CanonBall(ti,r_target,Target); % SRP force on the Target
% U_eci_Chasser = F_CanonBall(t,r_chasser,Chasser); % SRP force on the Chasser

F_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser); % F_CanonBall(t,r_chasser,Chasser); %  

%***************************************%
% Rotational state differential equation
thau_Chasser     = F_body(4:6,:);
I_Chasser        = diag(Chasser.Moment_Of_Inertia_Calculator());
Omega_matrix_Chasser = [0 -Omega_Chasser_body(1) -Omega_Chasser_body(2) -Omega_Chasser_body(3);
                        Omega_Chasser_body(1) 0 Omega_Chasser_body(3) -Omega_Chasser_body(2);
                        Omega_Chasser_body(2) -Omega_Chasser_body(3) 0 Omega_Chasser_body(1);
                        Omega_Chasser_body(3) Omega_Chasser_body(2) -Omega_Chasser_body(1) 0];

qdot(1:4,:) = 1/2*Omega_matrix_Chasser*q_Chasser;
Omega_dot(1:3,:) = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body) + thau_Chasser);

%***************************************%
% Translational state differential equation
CN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = CN\F_body(1:3,1);

%% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(r_target)*v_target; h_norm = norm(h_vec);  u_acc = U_eci_Target;
eh = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
V_rel = V_sun - v_target;
u_acc_dot = beta_chief*(V_rel/r_rel_norm^3 - 3*dot(Xrel,V_rel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );


%% Computing the relative acceleration
delf     =  -mu_Earth/rc_norm^3*r_c + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelU_srp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP

%% Integrating the relative trajectory
Rho_dot(1:3,:) = rho_prime;
Rho_dot(4:6,:) = DelU_srp + delf - 2*tilde(Omega)*rho_prime ...
                - tilde(Omegadot)*rho - tilde(Omega)*(tilde(Omega)*rho);

%% Collecting the rates of change
Xaug_dot = [qdot; Omega_dot; Rho_dot]*Tc;

end