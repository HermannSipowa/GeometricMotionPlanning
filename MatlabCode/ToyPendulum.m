function xdot = ToyPendulum(t,X)
r = X(1);
rdot = X(3);
thetadot = X(4);
xdot = [rdot, thetadot, r*thetadot^2, -(rdot^2 + thetadot^2)/r]';
end 