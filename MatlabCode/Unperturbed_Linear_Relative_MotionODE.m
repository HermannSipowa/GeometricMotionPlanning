function Xaug_dot = Unperturbed_Linear_Relative_MotionODE(t,Xaug,OE,delCOE,Mu_body)
format long

%% Collecting relevant quantities
theta = Xaug(1); delCOE(2) = Xaug(2);
delq1 = delCOE(4); delq2 = delCOE(5);


a = OE(1); q1 = OE(4); q2 = OE(5);

% Defining canstants
p = a*(1 - q1^2 - q2^2);
h = sqrt(Mu_body*p);
R = p/(1 + q1*cos(theta) + q2*sin(theta));
Vr = h/p * (q1*sin(theta) - q2*cos(theta));
Vt = h/p * (1 + q1*cos(theta) + q2*sin(theta));


% Create the first row for the forward mapping [A]
A = zeros(1,6);

% Nonzero the first row of the matrix elements of [A]
A(1,1) = R/a;
A(1,2) = Vr/Vt*R;
A(1,4) = -R/p * ( 2*a*q1 + R*cos(theta) );
A(1,5) = -R/p * ( 2*a*q2 + R*sin(theta) );


%% Integrating the chief true anomaly
Xdot(1) = h/R^2;

%% Integrating the relative tue latitude
delr = A*delCOE;
delp = p/a*delCOE(1) - 2*a*(q1*delq1 + q2*delq2);

Rho_dot(1) = h/R^2 * ( delp/(2*p) - 2*delr/R );

%% Collecting the rates of change
Xaug_dot = [Xdot; Rho_dot];

end