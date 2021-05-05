function u0 = mypdexic(x,X0,Xf,T)    % Initial condition for AGHF
%%  semi-straight line connecting X0 and Xf
% u0 = X0+(Xf-X0)*(x/T);
% u0(8:10) = X0(8:10)+(Xf(8:10)-X0(8:10))*(x/T) + 0.1*sin(2*pi*x/T);

%% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;
u0 = X0*cos(freq1*x/T) ...
    + Xf*sin(freq1*x/T) ...
    - (X0+Xf)*sin(freq2*x/T);
u0(1:4) = u0(1:4)/norm(u0(1:4)); % imposing the quaternion constraint
end