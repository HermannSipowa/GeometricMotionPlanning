function [delXf] = TH_Approximation(Delt,mu,X0,Xf,delX0)

OE0 = RVtoCOEs(X0(1:3),X0(4:6),mu); f0 = OE0(6);
OEf = RVtoCOEs(Xf(1:3),Xf(4:6),mu); f  = OEf(6);
A    = Afoward(f0);
Ainv = Areverse(f);
STMinv_f0 = STM(f0,Delt);
STM_f     = STMinv(f);
Phi_LVLH  = Ainv*STM_f*STMinv_f0*A;

delXf = Phi_LVLH*delX0;

end


function PhiMat = STM(f,t)

rho_f  = 1+e*cos(f);
phi1 = rho_f*sin(f);
phi1_prime = rho_f*cos(f)-e*(sin(f))^2;
S_phi1 = -cos(f)-1/2*e*(cos(f))^2;

%% Elliptic and hyperbolic orbits
if e==1
    Jf = 1/4*tan(f/2) - 1/20*tan(f)^5;
    phi2 = 2*phi1*Jf - cos(f)/rho_f;
    phi3 = -2*phi1*Jf - (cos(f))^2/rho_f - (cos(f))^2;
    phi2_prime = 2*phi1_prime*Jf + sin(f)*cos(f)/rho_f^2 - sin(f)/rho_f;
    S_phi2 = -rho_f^2 * Jf;
    phi3_prime = 2*( phi1_prime*S_phi2 - phi2_prime*S_phi1 );
    S_2phi3p1 = 2*( rho_f^2*Jf - sin(f) ) - sin(f)*cos(f);
else
    Df = sin(f)*(2+e*cos(f))/rho_f^2;
    Lf = sqrt(mu/p^3)*t;
    phi2 = e*phi1/(1-e^2) * (Df - 3*e*Lf) - cos(f)/rho_f;
    phi3 = -phi1/(1-e^2)  * (Df - 3*e*Lf) -  (cos(f))^2/rho_f - (cos(f))^2;
    phi2_prime = e*phi1_prime/(1-e^2) * (Df - 3*e*Lf) ...
                 + e*sin(f)*cos(f)/rho_f^2 + sin(f)/rho_f;
    
    S_phi2 = -rho_f^2/(2*(1-e^2)) * (Df - 3*e*Lf);
    phi3_prime = 2*(phi1_prime*S_phi2 - phi2_prime*S_phi1);
    S_2phi3p1  = ( e*sin(f)*(2+e*cos(f)) - 3*rho_f^2*Lf )/(1-e^2);
end  

PhiMat = [phi1         phi2         phi3         0   0        0;
       -2*S_phi1    -2*S_phi2    -S_2phi3p1   1   0        0;
       0            0            0            0   cos(f)   sin(f);
       phi1_prime   phi2_prime   phi3_prime   0   0        0;
       -2*phi1      -2*phi2      (-2*phi3-1)  0   0        0;
       0            0            0            0   -sin(f)  cos(f)];


end

function PhiMat = STMinv(f)

rho_f  = 1+e*cos(f);
phi1 = rho_f*sin(f);
phi1_prime = rho_f*cos(f)-e*(sin(f))^2;
S_phi1 = -cos(f)-1/2*e*(cos(f))^2;

%% Elliptic and hyperbolic orbits
if e==1
    Jf = 1/4*tan(f/2) - 1/20*tan(f)^5;
    phi2 = 2*phi1*Jf - cos(f)/rho_f;
    phi3 = -2*phi1*Jf - (cos(f))^2/rho_f - (cos(f))^2;
    phi2_prime = 2*phi1_prime*Jf + sin(f)*cos(f)/rho_f^2 - sin(f)/rho_f;
    S_phi2 = -rho_f^2 * Jf;
    phi3_prime = 2*( phi1_prime*S_phi2 - phi2_prime*S_phi1 );
    S_2phi3p1 = 2*( rho_f^2*Jf - sin(f) ) - sin(f)*cos(f);
else
    Df = sin(f)*(2+e*cos(f))/rho_f^2;
    Lf = sqrt(mu/p^3)*t;
    phi2 = e*phi1/(1-e^2) * (Df - 3*e*Lf) - cos(f)/rho_f;
    phi3 = -phi1/(1-e^2)  * (Df - 3*e*Lf) -  (cos(f))^2/rho_f - (cos(f))^2;
    phi2_prime = e*phi1_prime/(1-e^2) * (Df - 3*e*Lf) ...
                 + e*sin(f)*cos(f)/rho_f^2 + sin(f)/rho_f;
    
    S_phi2 = -rho_f^2/(2*(1-e^2)) * (Df - 3*e*Lf);
    phi3_prime = 2*(phi1_prime*S_phi2 - phi2_prime*S_phi1);
    S_2phi3p1  = ( e*sin(f)*(2+e*cos(f)) - 3*rho_f^2*Lf )/(1-e^2);
end  

PhiMat = [(4*S_phi2+phi2_prime)      0  0       -phi2  2*S_phi2    0;
          (-4*S_phi1-phi1_prime)     0  0       phi1   -2*S_phi1   0;
          -2                         0  0       0      1 0;
          (-2*S_2phi3p1-phi3_prime)  1  0       phi3   -S_2phi3p1  0;
          0                          0  cos(f)  0      0           -sin(f);
          0                          0  sin(f)  0      0           cos(f)];

end


function AMat = Afoward(f)

rho_f  = 1+e*cos(f);
const1 = ( mu*rho_f ) / h^(3/2);
const2 = -mu*e*sin(f) / h^(3/2);

AMat = [const1*eye(3)  zeros(3);
        const2*eye(3)  (1/const1)*eye(3)];

end


function AMat = Areverse(f)

rho_f  = 1+e*cos(f);
const1 = ( mu*rho_f ) / h^(3/2);
const2 = -mu*e*sin(f) / h^(3/2);

AMat = [(1/const1)*eye(3)  zeros(3);
        -const2*eye(3)     const1*eye(3)];

end












