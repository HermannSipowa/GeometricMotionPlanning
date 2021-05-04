function [Aug_Xdot] = ControlUnperturbed2Bodyfunc(t,Aug_X,mu_Earth,K,P)
format long
X_Chasser = Aug_X(1:6,1); X_Target = Aug_X(7:12,1); X_Desired = Aug_X(13:18,1);

%% Computing the globally asymptotically stabilizing controller
f_desired = -mu_Earth/norm(X_Desired(1:3,:))^3 * X_Desired(1:3,1);
f_deputy =  -mu_Earth/norm(X_Chasser(1:3,:))^3 * X_Chasser(1:3,1);
del_r = X_Chasser(1:3,:) - X_Desired(1:3,:);
del_v = X_Chasser(4:6,1) - X_Desired(4:6,1);
u = f_desired - f_deputy - K*del_r - P*del_v;


%% Integrating the deputy's trajectory
Aug_Xdot_Chasser(1:3,1) = X_Chasser(4:6,1);
Aug_Xdot_Chasser(4:6,1) = -mu_Earth/norm(X_Chasser(1:3,:))^3 * X_Chasser(1:3,1) + u;

%% Integrating the chief's trajectory
Aug_Xdot_Target(1:3,1) = X_Target(4:6,1);
Aug_Xdot_Target(4:6,1) = -mu_Earth/norm(X_Target(1:3,:))^3 * X_Target(1:3,1);

%% Integrating the reference's trajectory
Aug_Xdot_Desired(1:3,1) = X_Desired(4:6,1);
Aug_Xdot_Desired(4:6,1) = -mu_Earth/norm(X_Desired(1:3,:))^3 * X_Desired(1:3,1);

%% Collecting the rates of change
Aug_Xdot = [Aug_Xdot_Chasser; Aug_Xdot_Target; Aug_Xdot_Desired];

end