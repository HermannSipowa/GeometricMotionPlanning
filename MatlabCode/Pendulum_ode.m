dt=.001; % step size
t = 1:dt:5; % Calculates upto y(3)
y = zeros(3,length(t));
y(:,1) = [0 0 0]'; % initial condition
Sigma=[0 0 0]';
Omega=[1 0.5 -0.7]';
F_xy =@(t,Sigma) myMRPsODE(t,Sigma,Omega); % change the function as you desire

for i=1:(length(t)-1) % calculation loop
    k_1 = F_xy(t(i),y(:,i));
    k_2 = F_xy(t(i)+0.5*dt,y(:,i)+0.5*dt*k_1);
    k_3 = F_xy((t(i)+0.5*dt),(y(:,i)+0.5*dt*k_2));
    k_4 = F_xy((t(i)+dt),(y(:,i)+k_3*dt));
    
    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;  % main equation
    
    if norm(y(:,i+1))>1
        y(:,i+1)=-y(:,i+1)/norm(y(:,i+1))^2;
    end
end

dataPath = 'Outputs_Data/';
save([dataPath,'yData','.mat'], 'y','-v7.3')

