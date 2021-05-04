function [Aug_Xdot] = Unperturbed2Bodyfunc(t,Aug_X,mu_Earth)
format long
X_Chasser = Aug_X(1:6,1); X_Target = Aug_X(7:12,1);

%% Integrating the deputy's trajectory
Aug_Xdot_Chasser(1:3,1) = X_Chasser(4:6,1);
Aug_Xdot_Chasser(4:6,1) = -mu_Earth/norm(X_Chasser(1:3,:))^3*X_Chasser(1:3,1);

%% Integrating the chief's trajectory
Aug_Xdot_Target(1:3,1)  = X_Target(4:6,1);
Aug_Xdot_Target(4:6,1) = -mu_Earth/norm(X_Target(1:3,:))^3 * X_Target(1:3,1);

%% Collecting the rates of change
Aug_Xdot = [Aug_Xdot_Chasser; Aug_Xdot_Target];

end