function u0 = MyPendulumControlContraints_ic(x) % Initial condition for GHF
global X0 Xf xmax
u0 = [X0(1) + (Xf(1)-X0(1))*sin(5*pi*x/(2*xmax))
      X0(2) + (Xf(2)-X0(2))*sin(pi*x/(2*xmax));
      X0(3)*cos(pi*x/(2*xmax)) + Xf(3)*sin(pi*x/(2*xmax));
      X0(4)*cos(pi*x/(2*xmax)) + Xf(4)*sin(pi*x/(2*xmax));
      sin(x*pi/xmax);
      sin(x*pi/xmax)];

end
% --------------------------------------------------------------------------