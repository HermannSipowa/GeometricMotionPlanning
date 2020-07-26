dt=.001; % step size
t = 1:dt:5; % Calculates upto y(3)
y = zeros(3,length(t)); 
y(:,1) = [0 0 0]'; % initial condition
Sigma=[0 0 0]';
Omega=[1 0.5 -0.7]';
F_xy =@(t,Sigma) myODE(t,Sigma,Omega); % change the function as you desire

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
% plot(t,y,'.')
% xlabel('Time (s)')
% ylabel('SRP')
% title('Intrgation of SRP')

function Sigmadot=myODE(t,Sigma,Omega)

sigma=norm(Sigma);
sigma_tilde=[0 -Sigma(3) Sigma(2); Sigma(3) 0 -Sigma(1); -Sigma(2) Sigma(1) 0];
Sigmadot=1/4*((1-sigma^2)*eye(3)+2*sigma_tilde+2*(Sigma*Sigma'))*Omega;
end



% %% Example 1: 2D Surface in 3D Space
% %-------------------------------------------
% close all
% clear all
% clc
% start_up
% 
% c1 = rgb('RosyBrown'); c2 = rgb('Black'); c3 = rgb('Lime');
% c4 = rgb('Tomato'); c5 = rgb('DarkBlue'); c6 = rgb('LightSlateGray');
% 
% Smap  = 10;
% Straj = 15;
% basis = 20;
% xgrid = linspace(-Smap, Smap, 100);
% ygrid = linspace(-Smap, Smap, 100);
% [xi, yi] = meshgrid(xgrid, ygrid);
% 
% a = .01; b = -.01; c = .03; d = -.02; e = .04; f = 0;
% zi = a.*xi.^2 + b.*yi.^2 + c.*xi.*yi + d.*xi + e.*yi + f;
% zfunc = @(x,y) a.*x.^2 + b.*y.^2 + c.*x.*y + d.*x + e.*y + f;
% zgrid = zfunc(xgrid,ygrid);
% 
% %% Solve BVP
% %-------------------------------------------
% X0 = [-8,-7,0,0]';
% n = 4;
% Smap    = 3;
% opts    = bvpset('RelTol',2.53447e-07,'AbsTol',2.22045e-20,'Stats','on'); % 'FJacobian',@jac,
% xmesh   = linspace(-Smap, Smap, 10);
% solinit = bvpinit(xmesh, X0);
% sol5c   = bvp5c(@Pendulum_ode, @Pendulum_bcfcn, solinit, opts);
% 
% 
% 
% %% Ploting the time history of the system's states
% %-------------------------------------------
% ex_func = @(x,y) [1, 0, 2.*a.*x + c.*y + d] / norm([1, 0, 2.*a.*x + c.*y + d]);
% ey_func = @(x,y) [0, 1, 2.*b.*y + c.*x + e] / norm([0, 1, 2.*b.*y + c.*x + e]);
% 
% 
% X = X0(1);
% Y = X0(2);
% Z = zfunc(X,Y);
% ex = ex_func(X,Y);
% ey = ey_func(X,Y);
% 
% surf(xi,yi,zi, 'EdgeAlpha', 0)
% hold on
% 
% plt1 = scatter3(sol5c.y(1,1),sol5c.y(2,1),zfunc(sol5c.y(1,1),sol5c.y(2,1)),80);
% plt1.MarkerEdgeColor = c2;
% plt1.MarkerFaceColor = c2;
% plt1.Marker = 'o';
% 
% plt2 = scatter3(sol5c.y(1,end),sol5c.y(2,end),zfunc(sol5c.y(1,end),sol5c.y(2,end)),80);
% plt2.MarkerEdgeColor = c3;
% plt2.MarkerFaceColor = c3;
% plt2.Marker = 'o';
% 
% plt3 = plot3(sol5c.y(1,:),sol5c.y(2,:),zfunc(sol5c.y(1,:),sol5c.y(2,:)),'c');
% 
% plt4 = quiver3(X,Y,Z,ex(1),ex(2),ex(3),2);
% plt4.LineWidth = 1;
% plt4.Color = 'm';
% hold on
% plt5 = quiver3(X,Y,Z,ey(1),ey(2),ey(3),2);
% plt5.LineWidth = 1;
% plt5.Color = 'k';
% 
% r_unit = [X,Y,Z]/norm([X,Y,Z]);
% 
% colormap copper
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% set(gca,'view',[58 23])
% axis equal
% legend([plt1,plt2,plt3,plt4,plt5],{'Start','Goal','traj','$\hat{\textbf{\textit{u}}}$','$\hat{\textbf{\textit{v}}}$'})
% 
% %% Ploting the time history of the system's states
% %-------------------------------------------
% figure
% YLabel = {'$r(t)$','$\theta(t)$','$\dot{r}(t)$','$\dot{\theta}(t)$'};
% for i = 1:n
% subplot(2,2,i)
% plot(sol5c.x,sol5c.y(i,:),'r')
% xlabel('t')
% ylabel(YLabel(i))
% grid on
% end
% sgt = sgtitle('Example Toy Problem');
% sgt.FontSize = 20;
% 
% 
% %-------------------------------------------
% %-------------------------------------------
% function xdot = Pendulum_ode(t,X) % equation to solve
% 
% u = X(1); v = X(2); udot = X(3); vdot = X(4);
% a = .01; b = -.01; c = .03; d = -.02; e = .04;
% denominator = 1 + d^2 + e^2 + 2*c*u*(e + 2*(a + b)*v) + 2*d*(2*a*u + c*v)...
%                 + c^2*(u^2 + v^2) + 4*(a^2*u^2 + b*v*(e + b*v));
%             
% xdot(1) = udot;
% xdot(2) = vdot;
% xdot(3) = 2*(d + 2*a*u + c*v)*(a*udot^2 + vdot*(c*udot + b*vdot)) / denominator;
% xdot(4) = 2*(e + c*u + 2*b*v)*(a*udot^2 + vdot*(c*udot + b*vdot)) / denominator;
% 
% end
% %-------------------------------------------
% function res = Pendulum_bcfcn(ya,yb) % boundary conditions
% x0 = -8; y0 = -7; x_End = 6; y_End = 7;
% Start = [ x0; y0]; End = [x_End; y_End]; 
% 
% res = [ya(1:2) - Start
%        yb(1:2) - End];
% end
% %-------------------------------------------