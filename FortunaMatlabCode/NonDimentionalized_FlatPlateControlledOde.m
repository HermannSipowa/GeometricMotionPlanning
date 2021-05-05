function [Xdot] = NonDimentionalized_FlatPlateControlledOde(t,X,Target,Chasser)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

u = [u11(t) u12(t) u13(t)]';
% Collecting relevant quantities
X_chief = full([Chief_x(t);  Chief_y(t);  Chief_z(t);...
                Chief_dx(t); Chief_dy(t); Chief_dz(t)]);
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

q_Chasser = X(1:4)/norm(X(1:4));
Omega_Chasser_body = X(5:7);
rho = X(8:10);  rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% F_body = F_CanonBall(ti,r_chasser,Chasser); % SRP force on the Chasser
% F_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);

%***************************************%
% Rotational state differential equation
% thau_Chasser         = F_body(4:6);
I_Chasser            = diag(Chasser.Moment_Of_Inertia_Calculator());
Omega_matrix_Chasser = [0, -Omega_Chasser_body(1), -Omega_Chasser_body(2), -Omega_Chasser_body(3);
                        Omega_Chasser_body(1),  0,  Omega_Chasser_body(3), -Omega_Chasser_body(2);
                        Omega_Chasser_body(2), -Omega_Chasser_body(3),  0,  Omega_Chasser_body(1);
                        Omega_Chasser_body(3), Omega_Chasser_body(2), -Omega_Chasser_body(1), 0];

qdot      = 1/2*Omega_matrix_Chasser*q_Chasser;
Omega_dot = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body) + u);

%***************************************%
% Translational state differential equation
q_Chasser = [q1(t),q2(t),q3(t),q4(t)].'; q_Chasser = q_Chasser/norm(q_Chasser);
F_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\F_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;

% Integrating the relative trajectory
Rho_dot = [rho_prime; DelUsrp + del_ag + Coriolis_Effect];

% Collecting the rates of change
Xdot = [qdot; Omega_dot; Rho_dot]*Tc;

end