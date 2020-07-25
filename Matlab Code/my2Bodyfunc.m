function [Aug_Xdot] = my2Bodyfunc(t,Aug_X,Target,Chasser)
format long
global mu_Earth index acc2body time2body

Aug_X_Chasser = Aug_X(1:13); Aug_X_Target = Aug_X(14:26,1);
X_Chasser = Aug_X_Chasser(8:13,1); X_Target = Aug_X_Target(8:13,1);

q_Chasser = Aug_X_Chasser(1:4,1)/norm(Aug_X_Chasser(1:4,1));
Omega_Chasser_body = Aug_X_Chasser(5:7,1);

q_Target = Aug_X_Target(1:4,1)/norm(Aug_X_Target(1:4,1));
Omega_Target_body = Aug_X_Target(5:7,1);


%% Computing the dynamic of the CanonBall
%****************************************%
% Rotational state differential equation
U_eci_Target        = F_CanonBall(t,X_Target,Target); % F_SPR_FlatPlate(t, X_Target, Target, q_Target); % 
thau_Target         = zeros(3,1); % F_target(4:6,1); % thau_Target = zeros(3,1);
I_Target            = diag(Target.Moment_Of_Inertia_Calculator());
Omega_matrix_Target = [0 -Omega_Target_body(1) -Omega_Target_body(2) -Omega_Target_body(3);
    Omega_Target_body(1) 0 Omega_Target_body(3) -Omega_Target_body(2);
    Omega_Target_body(2) -Omega_Target_body(3) 0 Omega_Target_body(1);
    Omega_Target_body(3) Omega_Target_body(2) -Omega_Target_body(1) 0];

Aug_Xdot_Target(1:4,1) = 1/2*Omega_matrix_Target*q_Target;
Aug_Xdot_Target(5:7,1) = I_Target\(-cross(Omega_Target_body,I_Target*Omega_Target_body) + thau_Target);

%***************************************%
% Translational state differential equation
r_Target = norm(X_Target(1:3,:));

%***************************************%
% Integrating the velocity vector to get the possition at each time step
Aug_Xdot_Target(8,1)  = X_Target(4);
Aug_Xdot_Target(9,1)  = X_Target(5);
Aug_Xdot_Target(10,1) = X_Target(6);

%***************************************%
% Integrating the acceleration to the velocity at each time step
Aug_Xdot_Target(11,1) = -mu_Earth/(r_Target)^3*X_Target(1)+U_eci_Target(1);
Aug_Xdot_Target(12,1) = -mu_Earth/(r_Target)^3*X_Target(2)+U_eci_Target(2);
Aug_Xdot_Target(13,1) = -mu_Earth/(r_Target)^3*X_Target(3)+U_eci_Target(3);


%% Computing the dynamic of the Flat plate
%****************************************%
% Computing the force acting on the system
F_body = F_SPR_FlatPlate(t, X_Chasser, Chasser, q_Chasser); % F_CanonBall(t,X_Chasser,Chasser); %  

%***************************************%
% Rotational state differential equation
thau_Chasser  = F_body(4:6,:); 
I_Chasser     = diag(Chasser.Moment_Of_Inertia_Calculator());


Omega_matrix_Chasser = [0 -Omega_Chasser_body(1) -Omega_Chasser_body(2) -Omega_Chasser_body(3);
    Omega_Chasser_body(1) 0 Omega_Chasser_body(3) -Omega_Chasser_body(2);
    Omega_Chasser_body(2) -Omega_Chasser_body(3) 0 Omega_Chasser_body(1);
    Omega_Chasser_body(3) Omega_Chasser_body(2) -Omega_Chasser_body(1) 0];

Aug_Xdot_Chasser(1:4,1) = 1/2*Omega_matrix_Chasser*q_Chasser;
Aug_Xdot_Chasser(5:7,1) = I_Chasser\(-cross(Omega_Chasser_body,I_Chasser*Omega_Chasser_body) + thau_Chasser);

%***************************************%
% Translational state differential equation
CN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = CN\F_body(1:3,1);
acc2body(index,:) = [U_eci_Chasser; thau_Chasser];
time2body(index) = t;
index = index+1;
% U_eci_Chasser = F_CanonBall(t,X_Chasser,Chasser);
r_Chasser = norm(X_Chasser(1:3,:));

%***************************************%
% Integrating the velocity vector to get the possition at each time step
Aug_Xdot_Chasser(8,1)  = X_Chasser(4);
Aug_Xdot_Chasser(9,1)  = X_Chasser(5);
Aug_Xdot_Chasser(10,1) = X_Chasser(6);

%***************************************%
% Integrating the acceleration to the velocity at each time step
Aug_Xdot_Chasser(11,1) = -mu_Earth/(r_Chasser)^3*X_Chasser(1)+U_eci_Chasser(1);
Aug_Xdot_Chasser(12,1) = -mu_Earth/(r_Chasser)^3*X_Chasser(2)+U_eci_Chasser(2);
Aug_Xdot_Chasser(13,1) = -mu_Earth/(r_Chasser)^3*X_Chasser(3)+U_eci_Chasser(3);


%% Collecting the rates of change
Aug_Xdot = [Aug_Xdot_Chasser; Aug_Xdot_Target];

end