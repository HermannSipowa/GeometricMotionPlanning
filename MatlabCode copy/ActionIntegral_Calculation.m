global Mass g Penalty tMin tMax rMin rMax ...
    tClearance rClearance z1 z2 rr Alpha
for ii = 1:length(Tmesh)
    % [p1, dp1] = pdeval(m,x,sol(ii,:,1),x);
    % [p2, dp2] = pdeval(m,x,sol(ii,:,2),x);
    % [p3, dp3] = pdeval(m,x,sol(ii,:,3),x);
    % [p4, dp4] = pdeval(m,x,sol(ii,:,4),x);
    
    [p1, dp1] = pdeval(m,Xmesh,u1(ii,:),Xmesh);
    [p2, dp2] = pdeval(m,Xmesh,u2(ii,:),Xmesh);
    [p3, dp3] = pdeval(m,Xmesh,u3(ii,:),Xmesh);
    [p4, dp4] = pdeval(m,Xmesh,u4(ii,:),Xmesh);
    
    for jj = 1:length(Xmesh)
        r         = p1(jj);
        theta     = p2(jj);
        rdot      = p3(jj);
        thetadot  = p4(jj);
        
        
        %% Definint the Metric
        Tmin = tMin + tClearance; % Minimum angle
        Tmax = tMax - tClearance; % Minimum angle
        Rmin = rMin + rClearance; % Minimum radius
        Rmax = rMax - rClearance; % Minimum radius
        R    = rr + rClearance;
        tip  = [r.*cos(theta); r.*sin(theta)];
        obs1 = [z1(1).*cos(z1(2)); z1(1).*sin(z1(2))];
        obs2 = [z2(1).*cos(z2(2)); z2(1).*sin(z2(2))];
        l5   = norm(tip-obs1);
        l6   = norm(tip-obs2);
        
        
        if l5>=R
            lambda5  = 0;
            plambda5 = zeros(4,1);
        else
            lambda5  = Alpha*((l5^2-R^2)/l5^2)^2;
            plambda5 = Alpha*4*R^2*(l5^2-R^2)/l5^6*[r-z1(1)*cos(theta-z1(2));
                r*z1(1)*sin(theta-z1(2)); 0; 0];
        end
        
        if l6>=R
            lambda6  = 0;
            plambda6 = zeros(4,1);
        else
            lambda6  = Alpha*((l6^2-R^2)/l6^2)^2;
            plambda6 = Alpha*4*R^2*(l6^2-R^2)/l6^6*[r-z2(1)*cos(theta-z2(2));
                r*z2(1)*sin(theta-z2(2)); 0; 0];
        end
        
        if theta>=Tmin
            lambda1  = 0;
            plambda1 = zeros(4,1);
        else
            l1       = theta-Tmin;
            lambda1  = Alpha*((l1^2-tClearance^2)/l1^2)^2;
            plambda1 = Alpha*[0; 2*tClearance^2/(theta-Tmin)^3; 0; 0];
        end
        
        if theta<=Tmax
            lambda2  = 0;
            plambda2 = zeros(4,1);
        else
            l2       = theta-Tmax;
            lambda2  = Alpha*((l2^2-tClearance^2)/l2^2)^2;
            plambda2 = Alpha*[0; 2*tClearance^2/(theta-Tmax)^3; 0; 0];
        end
        
        if r>=Rmin
            lambda3  = 0;
            plambda3 = zeros(4,1);
        else
            l3       = r-Rmin;
            lambda3  = Alpha*((l3^2-rClearance^2)/l3^2)^2;
            plambda3 = Alpha*[2*rClearance^2/(r-Rmin)^3; 0; 0; 0];
        end
        
        if r<=Rmax
            lambda4  = 0;
            plambda4 = zeros(4,1);
        else
            l4       = r-Rmax;
            lambda4  = Alpha*((l4^2-rClearance^2)/l4^2)^2;
            plambda4 = Alpha*[2*rClearance^2/(r-Rmax)^3; 0; 0; 0];
        end
        
        lambda  = 1 +lambda1+lambda2+lambda3+lambda4+lambda5+lambda6;
        plambda = plambda1+plambda2+plambda3+plambda4+plambda5+plambda6;
        
        
        
        H = [Penalty 0 0 0;
            0 Penalty 0 0;
            0 0 Mass^2  0;
            0 0 0 (Mass*r)^2]; % Riemannian metric
        G = lambda*H;
        
        %--------------------------------------------------------------------------
        
        Xdot      = [dp1(jj) dp2(jj) dp3(jj) dp4(jj)]';
        fd        = [rdot; thetadot; r*thetadot^2 + g*cos(theta)/Mass;
                     -2*rdot*thetadot/r + g*sin(theta)/(Mass*r)];
        Lagrangian(jj) = 1/2*(Xdot - fd)'*G*(Xdot - fd);
    end
    ActionInt(ii) =  trapz(Xmesh,Lagrangian);
    Lagrangian = [];
end

%% Plotting the action integral vs iteration 
close all
loglog(Tmesh,ActionInt)
grid on
xlabel('Iteration space (s)')
ylabel('Action itegral')






