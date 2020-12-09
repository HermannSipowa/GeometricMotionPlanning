%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hermann Kaptui Sipowa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear all
close all
start_up
format long e
clc
global mu_dot tspan Period

%% Calculating the boundaries conditions
%=======================================%

mm = 6; % Size of the state space
mu_Earth   = 3.986004415E5; % Earth gravitational parameter
Num_Agent  = 6; % Number of agents in the system
Num_Agent1 = 4; % Number of agents in the first final-orbit
Num_Agent2 = 2; % Number of agents in the second final-orbit

% Specifying the chief's orbit
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a_chief      = 1.5E+4;   % Semi-major axis in Km
e_chief      = 0.0;      % Eccentricity
inc_chief    = 50;       % Inclination in deg0
BigOmg_chief = 10;       % RAAN in deg
LitOmg_chief = 10;       % AOP in deg
M_chief      = 0;        % True anomaly in deg


% Computing the nearly non-singular orbital elements
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
LitOmg_chief = deg2rad(LitOmg_chief); M_chief = deg2rad(M_chief);
COE = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief];
[Position_target,Velocity_target]  = COEstoRV(COE,mu_Earth); f_chief = M_chief;
q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
theta_chief = f_chief+LitOmg_chief;

% Specify the initial relative-orbital element of the deputies
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dela = 0;
dele = 1/(2*a_chief);
deli = 1/(4*a_chief);
delLitOmg = 0;
delBigOmg = 0;
delM = 0;

delq1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
delq2 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
deltheta = delLitOmg + delM; % Relative true latitude in rad
delCOE = [dela, deltheta, deli, delq1, delq2, delBigOmg]';


% Integrating the deputies' initial trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Period = 2*pi*sqrt(a_chief^3/mu_Earth);
IntTime = 1*Period;
tspan  = linspace(0,IntTime,10000);
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-30);
mu_dot = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
OE_chief = [a_chief, theta_chief, inc_chief, q1_chief, q2_chief, BigOmg_chief];
AMap = ForwardMapping(OE_chief, mu_Earth); % Linear mapping matrix
Xo = AMap*delCOE;
u = zeros(3,1);
[time, X_nom] = ode113(@(t,X)CW_ode(t,X,u,mu_dot),tspan,Xo,options);

% tspan_normaized = tspan/Period;
% Xnomalized = Xo/Period;
% [~, X_nomalized] = ode113(@(t,X)NormalizeedCW_ode(t,X,u,mu_dot,Period),tspan_normaized,Xnomalized,options);
% X_rescaled = Period.*X_nomalized;
% 
% plot3(X_nom(:,1),X_nom(:,2),X_nom(:,3),'b','LineWidth',2);
% hold on
% plot3(X_rescaled(:,1),X_rescaled(:,2),X_rescaled(:,3),'r--','LineWidth',2);
% grid on

%%

Xo = nan(mm,Num_Agent);
for i = 1:Num_Agent
    idx = find(time <= i*IntTime/Num_Agent, 1,'last');
    Xo(:,i) = X_nom(idx,:)';
    % if i == 1
    %     X0 = X_nom(idx,:)';
    % elseif i == 4
    %     Xf = X_nom(idx,:)';
    % end
end
[~, X_rel] = ode113(@(t,x)CW_ode(t,x,u,mu_dot),tspan,Xo,options);


% Specify the final relative-orbital elements of the deputies in both orbits
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dela = 0;
dele1 = 1/(6*a_chief);
deli1 = 1/(2*a_chief);
delLitOmg = 0;
delBigOmg = 0;
delM = 0;
delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
deltheta = delLitOmg + delM; % Relative true latitude in rad
delCOE_1 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';

delBigOmg = -pi*1E-5;
delM = pi*1E-6;
delLitOmg = pi*1E-5;
deltheta = delLitOmg + delM; % Relative true latitude in rad
dele2 = 1/(8*a_chief);
deli2 = 2/(10*a_chief);
delq1_2 = dele2*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
delq2_2 = dele2*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
delCOE_2 = [dela, deltheta, deli2, delq1_2, delq2_2, delBigOmg]';


% Integrating the deputy final trajectories
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Xo_1 = AMap*delCOE_1;
u = zeros(3,1);
[time, X_nom_1] = ode113(@(t,X)CW_ode(t,X,u,mu_dot),tspan,Xo_1,options);

Xo_1 = nan(mm,Num_Agent1);
for i = 1:Num_Agent1
    idx = find(time <= i*IntTime/Num_Agent1, 1,'last');
    Xo_1(:,i) = X_nom_1(idx,:)';
end
[~, X_rel_1] = ode113(@(t,x)CW_ode(t,x,u,mu_dot),tspan,Xo_1,options);

Xo_2 = AMap*delCOE_2;
[time, X_nom_2] = ode113(@(t,X)CW_ode(t,X,u,mu_dot),tspan,Xo_2,options);

Xo_2 = nan(mm,Num_Agent2);
for i = 1:Num_Agent2
    idx = find(time <= (2*i-1)*IntTime/(2*Num_Agent2), 1,'last');
    Xo_2(:,i) = X_nom_2(idx,:)';
end
[~, X_rel_2] = ode113(@(t,x)CW_ode(t,x,u,mu_dot),tspan,Xo_2,options);





%% Geometric Motion Planning
%===========================%
%---------------------------- AGHF parameters ----------------------------
% Boundary conditions
clc
jj = 2; idx = 1+mm*(jj-1):mm*jj; X0 = X_rel(1,idx)'/Period; 
jj = 1; idx = 1+mm*(jj-1):mm*jj; Xf = X_rel_1(1,idx)'/Period;
% smax for steady state
smax = 9e3;
% # of grids in t
tgrids = 500;
% # of grids in s
sgrids = 500;
% # of grids of integration for final trajectory extraction
intgrids = 1e5;
% # of states
N = length(X0);
% # of inputs
M = 3;
% motion duration
T = .5;
% penalty value (the lambda in paper)
Penalty = 2; 

% solve for trajectory, see "AGHF" function file for implementation details
[p, dp, xfull, xdrift, bufull, budrift, sol, s, t, tint, cost,X_ode45,p_initial] = AGHF(smax, tgrids, intgrids, sgrids, Penalty, X0, Xf, T, N, M);


%% plotting
%===========================%
%------------------- Plotting the extracted control -----------------------
idx = 0:1:intgrids-1;
figure('Name','Extracted control')
plot(idx,budrift);
set(gca, 'XScale', 'log')
title('actuated curve length vs s')
xlabel('s')
ylabel('curve length')
grid on;


figure('Name','actuated curve length vs s')
loglog(s,cost);
title('actuated curve length vs s')
xlabel('s')
ylabel('curve length')
grid on;

%%
u1 = Period.*sol(:,:,1);    
u2 = Period.*sol(:,:,2);    
u3 = Period.*sol(:,:,3);    
u4 = Period.*sol(:,:,4);    
u5 = Period.*sol(:,:,5);   

figure('Name','3D configuration space curve deformation')
X = u1(1,:);
Y = u2(1,:);
Z = u3(1,:);
h1 = plot3(X,Y,Z,'k','LineWidth',2);
hold on
plot3(X,Y,Z,'r','LineWidth',2);
% legend('$x-y-\theta$ trajectory','initial guess')
% axis([min(u1(end,:))-0.1,max(u1(end,:))+0.1,min(u2(end,:))-0.1,max(u2(end,:))+0.1,min(u3(end,:))-0.1,max(u3(end,:))+0.1]);
grid ON;
title('3D configuration space curve deformation');
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
%pause;

for i=1:sgrids
    X=u1(i,:);h1.XDataSource='X';
    Y=u2(i,:);h1.YDataSource='Y';
    Z=u3(i,:);h1.ZDataSource='Z';
    refreshdata(h1,'caller');
    drawnow;
    pause(0.05);
end

%% animation
%%
clear X_ode45
tolerance = 1e-16;
h = 1e-9;
X_ode45(:,1) = X0test;
xnew = linspace(0,T*Period,intgrids); % time grids for integration
for i = 1:intgrids-1
    u = budrift(:,i);
    [y_f, ~, h_next] = Runge_Kutta_Fehlberg_4_5(@(t,X)CW_ode(t,X,u,mu_dot),X_ode45(:,i),xnew(i),h,xnew(i+1),tolerance);
    h = h_next;
    X_ode45(:,i+1) = y_f;
end


%% animation
for Penalty = 1
c1 = rgb('Aqua');
c2 = rgb('Gold');
c3 = rgb('Lime');
c4 = rgb('LightBlue');
c5 = rgb('DarkBlue');
c6 = rgb('Red');

cMat1 = [c6;c5;c2;c3];
cMat2 = [c1;c4];
Label_1 = {'Agent1','Agent2','Agent3','Agent4'};
Label_2 = {'Agent5','Agent6'};
Label = {'Agent1','Agent2','Agent3','Agent4','Agent5','Agent6'};

fig = figure('Name','3D view of the motion');
hold on
for jj=1%:Num_Agent1
    idx = 1+mm*(jj-1):mm*jj;

    plot3(X_rel_1(:,idx(1)),X_rel_1(:,idx(2)),X_rel_1(:,idx(3)),'m-.');
    plot3(X_rel_1(1,idx(1)),X_rel_1(1,idx(2)),X_rel_1(1,idx(3)),'o',...
        'LineWidth',3,...
        'MarkerEdgeColor',cMat1(jj,:),...
        'MarkerFaceColor',cMat1(jj,:)',...
        'MarkerSize',8);
    hold on

    xlabel('X [km]')
    ylabel('Y [km]')
    zlabel('Z [km]')
    title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
    if jj == 1
        Xf = X_rel_1(1,idx)';
    end
    % elseif jj == 2
    %     Xf = [Xf;X_rel_1(1,idx1)'];
    % end
end


%     for jj=1%:Num_Agent2
%         idx2 = 1+mm*(jj-1):mm*jj;
%         plot3(X_rel_2(:,idx2(1)),X_rel_2(:,idx2(2)),X_rel_2(:,idx2(3)),'m-.');
%         plot3(X_rel_2(1,idx2(1)),X_rel_2(1,idx2(2)),X_rel_2(1,idx2(3)),'o',...
%             'LineWidth',3,...
%             'MarkerEdgeColor',cMat2(jj,:),...
%             'MarkerFaceColor',cMat2(jj,:)',...
%             'MarkerSize',8);
%     end
%     


c1 = rgb('Aqua'); c2 = rgb('Gold'); c3 = rgb('Lime');
c4 = rgb('LightBlue'); c5 = rgb('DarkBlue'); c6 = rgb('Red');
c7 = rgb('Tomato');
cMat = [c6;c5;c2;c3;c1;c4];
plt   = zeros(Num_Agent);
%     figure
hold on
for jj=2%:Num_Agent
    idx = 1+mm*(jj-1):mm*jj;

    % h4 = plot3(X_rel(:,idx(1)),X_rel(:,idx(2)),X_rel(:,idx(3)),'k-.','DisplayName','Chief');
    plt(jj) = plot3(X_rel(1,idx(1)),X_rel(1,idx(2)),X_rel(1,idx(3)),'*',...
        'LineWidth',3,...
        'MarkerEdgeColor',cMat(jj,:),...
        'MarkerFaceColor',cMat(jj,:)',...
        'MarkerSize',8);

    xlabel('X [km]')
    ylabel('Y [km]')
    zlabel('Z [km]')
    title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
    % hL = legend(plt(1:jj),Label(1:jj));
    % axis equal
    % if jj == 5
    %     X0 = X_rel(1,idx)';
    % elseif jj == 6
    %     X0 = [X0; X_rel(1,idx)'];
    % end
    if jj == 2
        X0 = X_rel(1,idx)';
    end
end
grid on
h5 = plot3(0,0,0,'bo',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',15,'DisplayName','Chief');
h4 = plot3(X_rel(:,idx(1)),X_rel(:,idx(2)),X_rel(:,idx(3)),'k','DisplayName','Initial Traj');
h6 = plot3(X_rel_1(:,idx(1)),X_rel_1(:,idx(2)),X_rel_1(:,idx(3)),'m','DisplayName','Final Traj');

% Programatically move the Legend
newPosition = [.2 .7 0.05 0.05];
newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
view(31,18)
    
plot3(p_initial(1,:),p_initial(2,:),p_initial(3,:),'c-.','LineWidth',2);
% plot3(interp1(tint,xdrift(1,:),t),interp1(tint,xdrift(2,:),t),interp1(tint,xdrift(3,:),t),':c','LineWidth',2);
plot3(interp1(tint,X_ode45(1,:),t),interp1(tint,X_ode45(2,:),t),interp1(tint,X_ode45(3,:),t),'--r','LineWidth',2);
% legend('AGHF solution','initial Condition', 'integrated path',  'integrated ODE45', 'AutoUpdate','off')

end

% %%
% figure('Name','xy-plane trajectory deformation')
% 
% X=u1(1,:);
% Y=u2(1,:);
% h1=plot(X,Y,'k','LineWidth',2);
% hold on
% plot(X,Y,'r','LineWidth',2);
% legend('x-y trajectory','initial guess')
% % axis([min(u1(end,:))-0.1,max(u1(end,:))+0.1,min(u2(end,:))-0.1,max(u2(end,:))+0.1]);
% grid ON;
% title('xy-plane trajectory and unicycle deformation');
% xlabel('x');
% ylabel('y');
% %pause;
% 
% for i=1:sgrids
%     
%     X=u1(i,:);h1.XDataSource='X';
%     Y=u2(i,:);h1.YDataSource='Y';
%     
%     refreshdata(h1,'caller');
%     drawnow;
%     pause(0.05);
%     
% end
% 
% %%
% % plot actual trajectory
% Xs=interp1(tint,xdrift(1,:),t);
% Ys=interp1(tint,xdrift(2,:),t);
% Zs=interp1(tint,xdrift(3,:),t);
% Z1s=interp1(tint,xdrift(4,:),t);
% Z2s=interp1(tint,xdrift(5,:),t);
% 
% w = 0.2;
% h = 0.1;
% Rot_s = [cos(Zs(1)) -sin(Zs(1));
%        sin(Zs(1)) cos(Zs(1))];
% p1s = [Xs(1),Ys(1)]' + Rot_s*[w/2;h/2];
% p2s = [Xs(1),Ys(1)]' + Rot_s*[-w/2;h/2];
% p3s = [Xs(1),Ys(1)]' + Rot_s*[-w/2;-h/2];
% p4s = [Xs(1),Ys(1)]' + Rot_s*[w/2;-h/2];
% h_box_patch=patch([p1s(1) p2s(1) p3s(1) p4s(1)], [p1s(2) p2s(2) p3s(2) p4s(2)],[0.4,0.4,1]);
% h_v=quiver(Xs(1),Ys(1),Z1s(1)*cos(Zs(1)),Z1s(1)*sin(Zs(1)),'r');
% qfactor = 0.05;
% qhead = 10;
% qwidth = 1.5;
% set(h_v,'AutoScale','off','AutoScaleFactor',qfactor,'MaxHeadSize',qhead,'linewidth',qwidth)
% 
% axis equal
% axis([-1.1,1.1,-1.6,0.6]);
% grid ON;
% xlabel('x');
% ylabel('y');
% %pause;
% 
% ind_plot = 0:10:tgrids;
% ind_plot(1) = 1;
% ind_count = 1;
% for i=1:tgrids
%     delete(h_box_patch);
%     delete(h_v);
%     Rot_s = [cos(Zs(i)) -sin(Zs(i));
%         sin(Zs(i)) cos(Zs(i))];
%     p1s = [Xs(i),Ys(i)]' + Rot_s*[w/2;h/2];
%     p2s = [Xs(i),Ys(i)]' + Rot_s*[-w/2;h/2];
%     p3s = [Xs(i),Ys(i)]' + Rot_s*[-w/2;-h/2];
%     p4s = [Xs(i),Ys(i)]' + Rot_s*[w/2;-h/2];
%     pxs = [p1s(1) p2s(1) p3s(1) p4s(1) p1s(1)];
%     pys = [p1s(2) p2s(2) p3s(2) p4s(2) p1s(2)];
%     
%     h_box_patch=patch(pxs(1:4),pys(1:4),[0.4,0.4,1]);
%     h_v=quiver(Xs(i),Ys(i),Z1s(i)*cos(Zs(i)),Z1s(i)*sin(Zs(i)),'r');
%     set(h_v,'AutoScaleFactor',qfactor,'MaxHeadSize',qhead,'linewidth',qwidth)
%     
%     if i == ind_plot(ind_count)
%         patch(pxs(1:4),pys(1:4),[0.4,0.4,1],'EdgeAlpha', 0.4+(i/tgrids)*0.6, 'FaceAlpha', 0.4+(i/tgrids)*0.6);
%         U = Z1s(i)*cos(Zs(i));
%         V = Z1s(i)*sin(Zs(i));
%         if (abs(U)>0.0001)&&(abs(V)>0.0001)
%             quiver(Xs(i),Ys(i),U,V,'r','AutoScaleFactor',qfactor,'MaxHeadSize',qhead,'linewidth',qwidth);
%         end
%         ind_count = ind_count + 1;
%     end
% 
%     drawnow;
%     pause(0.005)
% end
% 
% delete(h_box_patch);
% delete(h_v);
% 





