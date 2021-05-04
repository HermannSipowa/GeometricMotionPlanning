function [xdot] = MyPendulumControlContraints_Dynamics(t,X,U)
global Mass g
r        = X(1);
theta    = X(2);
rdot     = X(3);
thetadot = X(4);
fd = [rdot; thetadot; r*thetadot^2 + g*cos(theta)/Mass;
      -2*rdot*thetadot/r + g*sin(theta)/(Mass*r)];
f1 = [0; 0; 1/Mass; 0];
f2 = [0; 0; 0; 1/(r*Mass)];

xdot = fd + f1*U(1) + f2*U(2);

end