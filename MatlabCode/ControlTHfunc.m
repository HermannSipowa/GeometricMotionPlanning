function [Aug_Xdot] = ControlTHfunc(t,Aug_X,K,P)
format long
global e_chief
X_Chasser = Aug_X(1:6,1); X_Desired = Aug_X(7:12,1);

%% Computing the globally asymptotically stabilizing controller
k = 1+e_chief*cos(t);
A = [zeros(3)         eye(3);
     3/k  0    0    0   2  0;
     0    0    0   -2   0  0;
     0    0   -1    0   0  0];
B = [zeros(3); eye(3)];
f_deputy  = A*X_Chasser;
f_desired = A*X_Desired;
del_r = X_Chasser(1:3) - X_Desired(1:3);
del_v = X_Chasser(4:6) - X_Desired(4:6);
u = f_desired(4:6) - f_deputy(4:6) - K*del_r - P*del_v;

%% Integrating the deputy's trajectory
Aug_Xdot_Chasser = f_deputy + B*u;

%% Integrating the reference's trajectory
Aug_Xdot_Desired = f_desired;

%% Collecting the rates of change
Aug_Xdot = [Aug_Xdot_Chasser; Aug_Xdot_Desired];

end