close all
clear all
clc
x_i = 1/(3*pi);
x_f = 1;
opts = bvpset('FJacobian',@jac,'RelTol',2.53447e-07,'AbsTol',2.22045e-20,'Stats','on');
xmesh = linspace(x_i, x_f, 10);
solinit = bvpinit(xmesh, [1; 1]);
sol4c = bvp4c(@bvpfcn, @bcfcn, solinit, opts);
sol5c = bvp5c(@bvpfcn, @bcfcn, solinit, opts);

xplot = linspace(x_i, x_f, 500);
yplot = [sin(1./xplot); -cos(1./xplot)./xplot.^2];

subplot(2,1,1)
plot(xplot,yplot(1,:),'k',sol4c.x,sol4c.y(1,:),'r*',sol5c.x,sol5c.y(1,:),'bo')
% legend('True','BVP4C','BVP5C')
xlabel('x')
ylabel('Position')

subplot(2,1,2)
plot(xplot,yplot(2,:),'k',sol4c.x,sol4c.y(2,:),'r*',sol5c.x,sol5c.y(2,:),'bo')
% legend('True','BVP4C','BVP5C')
xlabel('x')
ylabel('Velocity')
sgt = sgtitle('Comparison of BVP Solvers with Crude Error Tolerance');
sgt.FontSize = 20;

%-------------------------------------------
function dydx = bvpfcn(x,y) % equation to solve
dydx = [y(2)
       -2*y(2)/x - y(1)/x^4];
end
%-------------------------------------------
function res = bcfcn(ya,yb) % boundary conditions
res = [ya(1)
       yb(1)-sin(1)];
end
%-------------------------------------------
function dfdy = jac(x,y) % analytical jacobian for f
dfdy = [0      1
       -1/x^4 -2/x];
end
%-------------------------------------------



% %% Toy Problem: Pendulum
% options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-20);
% tvec = linspace(0,100,1000);
% r0 = 10; theta0 = 3*pi/2; rdot0 = 0.02; thetadot0 = 0.01;
% X0 = [r0, theta0, rdot0, thetadot0];
% [~, Xstate] = ode113(@(t,X)ToyPendulum(t,X),tvec,X0,options);
% 
% figure
% Label = {'$r(t)$','$\theta(t)$','$\dot{r}(t)$','$\dot{\theta}(t)$'};
% plt   = zeros(4);
% for k = 1:4
%     subplot(2,2,k)
%     plt(k) = plot(tvec,Xstate(:,k));
%     grid on
%     xlabel('t')
%     ylabel(Label(k))
% end
