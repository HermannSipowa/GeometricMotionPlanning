for i = 1:length(tspan)
    IC(i,:) = TahoeIC(tspan(i));
end
%%
close all
start_up
c1 = rgb('Black'); c2 = rgb('Cyan'); c3 = rgb('Goldenrod');
c4 = rgb('RosyBrown'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');

%~~~~~~~~~~~~~~~~~~~~~~~~~~% sol(end,:,i)
% Plotting the 3D homotopies
figure('Renderer', 'painters', 'Position', [15 15 900 600])
h1 = plot3(IC(:,8),IC(:,9),IC(:,10),'r--');
hold on
h2 = plot3(IC(1,8),IC(1,9),IC(1,10),'o','color', c2,...
    'LineWidth',2,...
    'MarkerEdgeColor',c2,...
    'MarkerFaceColor',c2,...
    'MarkerSize',5);
h3 = plot3(IC(end,8),IC(end,9),IC(end,10),'go',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',c1,...
    'MarkerSize',5);
h5 = plot3(0,0,0,'bo',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',c5,...
    'MarkerSize',15);
h = legend([h1, h5, h2, h3],{'Relative','chief','$\delta X_0$','$\delta X_f$'});
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')

figure
plot(tspan,IC(:,8))
ylabel('Xrel')
figure
plot(tspan,IC(:,9))
ylabel('Yrel')
figure
plot(tspan,IC(:,10))
ylabel('Zrel')

%%
    close all
    clc
    start_up
    c1 = rgb('Black'); c2 = rgb('Cyan'); c3 = rgb('Goldenrod');
    c4 = rgb('RosyBrown'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~% sol(end,:,i)
    % Plotting the 3D homotopies
    figure('Renderer', 'painters', 'Position', [15 15 900 600])
    h1 = plot3(sol(end,:,8),sol(end,:,9),sol(end,:,10),'r');
    hold on
    h2 = plot3(sol(1,1,8),sol(1,1,9),sol(1,1,10),'o','color', c2,...
        'LineWidth',2,...
        'MarkerEdgeColor',c2,...
        'MarkerFaceColor',c2,...
        'MarkerSize',5);
    h3 = plot3(sol(1,end,8),sol(1,end,9),sol(1,end,10),'go',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',c1,...
        'MarkerSize',5);
    h5 = plot3(0,0,0,'bo',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',c5,...
        'MarkerSize',15);
    %arrow3D([.5e4,.5e4,.5e4] ,8e3*XSUN_LVLH,c3)
    h = legend([h1, h5, h2, h3],{'Relative','chief','$\delta X_0$','$\delta X_f$'});
    xlabel('X [km]')
    ylabel('Y [km]')
    zlabel('Z [km]')
    
    
    
    %% Visualizing the results
%-------------------------------------------
 for k = 1
    close all
    clc
    start_up
    c1 = rgb('Black'); c2 = rgb('Cyan'); c3 = rgb('Goldenrod');
    c4 = rgb('Lime'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
    c7 = rgb('Crimson'); c8 = rgb('DarkSlateGray');
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~% sol(end,:,i)
    % Plotting the 3D homotopies
    figure('Renderer', 'painters', 'Position', [15 15 900 600])
    h1 = plot3(sol(1,:,8),sol(1,:,9),sol(1,:,10),'r');
    hold on
    h2 = plot3(sol(1,1,8),sol(1,1,9),sol(1,1,10),'o','color', c2,...
        'LineWidth',2,...
        'MarkerEdgeColor',c2,...
        'MarkerFaceColor',c2,...
        'MarkerSize',5);
    h3 = plot3(sol(1,end,8),sol(1,end,9),sol(1,end,10),'go',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',c1,...
        'MarkerSize',5);
    h5 = plot3(0,0,0,'bo',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',c8,...
        'MarkerSize',15);
    %arrow3D([.5e4,.5e4,.5e4] ,8e3*XSUN_LVLH,c3)
    h = legend([h1, h5, h2, h3],{'Relative','chief','$\delta X_0$','$\delta X_f$'});
    xlabel('X [km]')
    ylabel('Y [km]')
    zlabel('Z [km]')
    for i=1:10
        idx= 1000*i;
        X = sol(idx,:,8)  ;h1.XDataSource='X';
        Y = sol(idx,:,9)  ;h1.YDataSource='Y';
        Z = sol(idx,:,10) ;h1.ZDataSource='Z';
        
        refreshdata(h1,'caller');
        drawnow;
        
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %% Plotting the 3D trajectory
%     close all
    for ii =1
        figure('Renderer', 'painters', 'Position', [85 85 1100 700])
        for kk = 1:4
            subplot(2,2,kk)
            for jj=1
                hold on
                h0 = plot3(XUncontrolled(:,8),XUncontrolled(:,9),XUncontrolled(:,10),'color', c8);
                h1 = plot3(Xstar_f(:,8),Xstar_f(:,9),Xstar_f(:,10),'color', c7);
                h2 = plot3(Xrel0(8),Xrel0(9),Xrel0(10),'o',...
                    'LineWidth',2,...
                    'MarkerEdgeColor',c6,...
                    'MarkerFaceColor',c6,...
                    'MarkerSize',5);
                h3 = plot3(Xrelf(8),Xrelf(9),Xrelf(10),'o',...
                    'LineWidth',2,...
                    'MarkerEdgeColor',c4,...
                    'MarkerFaceColor',c4,...
                    'MarkerSize',5);
                h5 = plot3(0,0,0,'bo',...
                    'LineWidth',2,...
                    'MarkerEdgeColor',c1,...
                    'MarkerFaceColor',c1,...
                    'MarkerSize',15);
                h = legend([h0,h1, h5, h2, h3],{'Uncontrolled','Controlled','chief','$\delta X_0$','$\delta X_f$'});
                rect = [0.48, 0.46, 0.08, 0.08];
                set(h, 'Position', rect)
                box on
                grid on
                xlabel('X [km]')
                ylabel('Y [km]')
                zlabel('Z [km]')
            end
            
            if kk == 1
                view(70,10)
            elseif kk==2
                view(0,90)  % XY
            elseif kk==3
                view(0,0)   % XZ
            else
                view(90,0)  % YZ
            end
            
        end
        sgt = sgtitle('Relative Trajectory','interpreter','tex');
        sgt.FontSize = 30;
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %% Plotting the difference between ECI Cartesian coordinates of the deputy
    for jj=1
        
        DQ1 = [Xstar_f(:,8) Xstar_f(:,11) Xstar_f(:,9) Xstar_f(:,12) Xstar_f(:,10) Xstar_f(:,13)];
        figure('Renderer', 'painters', 'Position', [15 15 900 600])
        YLabel={'$\delta x [km]$',...
            '$\delta \dot{x} [km/s]$',...
            '$\delta y [km]$',...
            '$\delta \dot{y} [km/s]$',...
            '$\delta z [km]$',...
            '$\delta \dot{z} [km/s]$'};
        for i=1:6
            subplot(3,2,i)
            plot1 = plot(tspan/86400,DQ1(:,i),'color', c4);
            ylabel(YLabel(i))
            xlabel('days')
            grid on
        end
        sgt = sgtitle('Relative Position','interpreter','tex');
        sgt.FontSize = 30;
    end
    
    
    for jj=1
        figure('Renderer', 'painters', 'Position', [15 15 900 600])
        YLabel={'$\tau_{_{\!1}}$','$\tau_{_{\!2}}$','$\tau_{_{\!3}}$'};
        for i=1:3
            subplot(3,1,i)
            plot(tspan/86400,Tau_ctrl(i,:),'color', c6)
            ylabel(YLabel(i))
            xlabel('days')
            grid on
            % legend('Classical Dynamics', 'Dual Quaternions')
        end
        sgt = sgtitle('Acceleration and Torque on the Flat-Plate','interpreter','tex');
        sgt.FontSize = 30;
    end
    
    
    for jj=1
        Q1 = Xstar_f(:,1:4);
        figure('Renderer', 'painters', 'Position', [15 15 900 600])
        YLabel={'$q_{_{\!1}}$','$q_{_{\!2}}$','$q_{_{\!3}}$','$q_{_{\!4}}$'};
        for i=1:4
            subplot(2,2,i)
            plot(tspan/86400,Q1(:,i),'color', c5) % -Q2(:,i)
            ylabel(YLabel(i))
            xlabel('days')
            grid on
            % legend('Classical Dynamics', 'Dual Quaternions')
        end
        sgt = sgtitle('Deputy Attitude (Quaternions)','interpreter','tex');
        sgt.FontSize = 30;
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Plotting the angular velocy
    for jj=1
        Q2 = Xstar_f(:,5:7);
        figure('Renderer', 'painters', 'Position', [15 15 900 600])
        YLabel={'$\omega_{_{\!1}}$','$\omega_{_{\!2}}$','$\omega_{_{\!3}}$'};
        for i=1:3
            subplot(3,1,i)
            % plot(time/86400,omega_chasser1(:,i)-omega_target1(:,i))
            plot(tspan/86400,Q2(:,i),'color', c5)
            ylabel(YLabel(i))
            xlabel('days')
            grid on
            % legend('Classical Dynamics', 'Dual Quaternions')
        end
        sgt = sgtitle('Deputy Angular Velocity','interpreter','tex');
        sgt.FontSize = 30;
    end
    
end % Plots for relative motions (Chaser and Target) with classical dynamics




%% Check for the normalization
COE = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief];
[Position_target,Velocity_target]  = COEstoRV(COE,mu_Earth);

Alpha = a_chief*eye(3); Tc = Period;
X0_normilized = [Alpha\Position_target; Alpha\Velocity_target*Tc]; T_normalized = tspan/Tc;
[~, X_nomalized] = ode113(@(t,X)NonDimentionalized_2BP(t, X, mu_Earth, Alpha, Tc),T_normalized,X0_normilized,options);
Xcheck_Pos = (Alpha*X_nomalized(:,1:3).').'; Xcheck_Vel = (1/Tc)*(Alpha*X_nomalized(:,4:6).').';

plot3(X_nom(:,1),X_nom(:,2),X_nom(:,3),'r')
hold on
plot3(Xcheck_Pos(:,1),Xcheck_Pos(:,2),Xcheck_Pos(:,3),'b-.')



clc
Rc = a_chief*eye(3); Rrho = eye(3); Tc = Period; T_normalized = tspan/Tc;
% Integrating the relative motion in the chief's LVLH frame
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
for i = 1
    
    X_aug0 = [Position_target; Velocity_target; Xo_1(:,1)];
    [~, X_rel] = ode113(@(t,X_aug)UnperturbedRelativeMotionODE(t, X_aug,mu_Earth),tspan,X_aug0,options);
    
    X0_normilized = [Rc\Position_target; Rc\Velocity_target*Tc; Rrho\Xo_1(1:3,1); Rrho\Xo_1(4:6,1)*Tc];
    [~, Xnom_rel] = ode113(@(t,Xaug)NonDimentionalized_RelativeODE(t, Xaug, mu_Earth, Rc, Rrho, Tc),T_normalized,X0_normilized,options);
    Xchief_Pos = (Rc*Xnom_rel(:,1:3).').'; Xchief_Vel = (1/Tc)*(Rc*Xnom_rel(:,4:6).').';
    Xrel_Pos = (Rrho*Xnom_rel(:,7:9).').'; Xrel_Vel = (1/Tc)*(Rrho*Xnom_rel(:,10:12).').';
end

% Plotting error between the two approaches
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
for k = 1
    plot3(X_rel(:,1),X_rel(:,2),X_rel(:,3),'r')
    hold on
    plot3(Xchief_Pos(:,1),Xchief_Pos(:,2),Xchief_Pos(:,3),'b-.')
    
    figure
    plot3(X_rel(:,7),X_rel(:,8),X_rel(:,9),'r')
    hold on
    plot3(Xrel_Pos(:,1),Xrel_Pos(:,2),Xrel_Pos(:,3),'b-.')
    %~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Plotting the difference between ECI Cartesian coordinates of the deputy
    for jj=1
        
        DQ1 = [X_rel(:,7) X_rel(:,10) X_rel(:,8) X_rel(:,11) X_rel(:,9) X_rel(:,12)];
        DQ2 = [Xrel_Pos(:,1) Xrel_Vel(:,1) Xrel_Pos(:,2) Xrel_Vel(:,2) Xrel_Pos(:,3) Xrel_Vel(:,3)];
        figure
        YLabel={'$\delta x_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!True}}} - \delta x_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!\!\!Approx}}} ~[km]$',...
            '$\delta \dot{x}_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!True}}} - \delta \dot{x}_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!\!\!Approx}}} ~[km]$',...
            '$\delta y_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!True}}} - \delta y_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!\!\!Approx}}} ~[km]$',...
            '$\delta \dot{y}_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!True}}} - \delta \dot{y}_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!\!\!Approx}}} ~[km]$',...
            '$\delta z_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!True}}} - \delta z_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!\!\!Approx}}} ~[km]$',...
            '$\delta \dot{z}_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!True}}} - \delta \dot{z}_{_{_{\!\!\!\!\!\!\!\!\!\!\!\!\!\!Approx}}} ~[km]$'};
        for i=1:6
            subplot(3,2,i)
            plot1 = plot( time/Period, DQ1(:,i) - DQ2(:,i) );
            ylabel(YLabel(i))
            xlabel('Period')
            grid on
        end
        % suptitle('Trajectories comparaison (Difference between truth and approximation)')
    end
    
    
end


%%
for i = 1
% X0 = [X0_Chief; Xrel]; X0_array = dlarray(X0);
% [fval,gradval] = dlfeval(@RelativeMotionODE,X0_array);
% function [f,grad] = RelativeMotionODE(Xaug)
% global mu_Earth
% 
% tilde = @(v) [0    -v(3)  v(2);
%               v(3)   0   -v(1);
%              -v(2)  v(1)    0];
% 
% % Collecting relevant quantities
% Target_States = Xaug(1:6,:);
% rho_aug       = Xaug(7:12,:);
% 
% r_target = Target_States(1:3,:); v_target = Target_States(4:6,:);
% rt_norm  = (r_target.'*r_target).^(1/2);
% rho = rho_aug(1:3,:); rho_prime = rho_aug(4:6,:);
% r_c = [rt_norm+rho(1) rho(2) rho(3)]';
% rc_norm = (r_c.'*r_c).^(1/2); % RIC position of the chasser
% TN  = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame
% 
% % Computing the angular velocity and angular acceleration of the target frame
% h_vec = tilde(r_target)*v_target;
% Omega = TN*( h_vec/rt_norm^2);
% Omegadot = TN*( -2*dot(r_target,v_target)*h_vec/rt_norm^4 );
% 
% % Computing the relative acceleration
% delf = -mu_Earth/rc_norm^3*r_c + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
% 
% % Integrating the the target trajectory
% f(1:3,:) = v_target;
% f(4:6,:) = -mu_Earth/rt_norm^3*r_target;
% 
% 
% % Integrating the relative trajectory
% f(7:9,:)   = rho_prime;
% f(10:12,:) = delf - 2*tilde(Omega)*rho_prime ...
%     - tilde(Omega)*(tilde(Omega)*rho) ...
%     - tilde(Omegadot)*rho;
% 
% % Numerically calculating gradient of f w.r.t the state X
% n = length(Xaug);
% grad = nan(n,n);
% for i = 1:n
%     grad(i,:) = dlgradient(f(i),Xaug)';
% end
% 
% end
end

%% Specify the initial relative-orbital element of the deputies
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% dela = zeros(4,1);
% dele = 3e-6*ones(4,1);
% deli = 1e-5*[1;-1;1;-1];
% delLitOmg = 5e-5*[-3, 4, 9, -2];
% delBigOmg = 1e-5*[1, 2, 3, 8];
% delM = 0*5e-7*[-1;1;0;0];
% mm = size(Xnom,2);
% Xint = zeros(mm,Num_Agent);
% for i = 1:Num_Agent
%    DelOE     = [dela(i); dele(i); deli(i);...
%                 delLitOmg(i); delBigOmg(i); delM(i)];
%    Xrel      = DelOEs_to_DelX(DelOE,OE_Chief,mu_Earth);
%    Xint(:,i) = [Rrho\Xrel(1:3,1); Rrho\Xrel(4:6,1)*Tc];
% end
% X01 = Xint(:,1); X02 = Xint(:,2);



%%

for i = 1:length(x)
    t_i = x(i);
    Xint(:,i) = mypdexic(t_i, X0, Xf, T);
end
% use spline to interpolate state from steady state pde solution
dpdx  = zeros(N,xpoints);
dx    = diff([x(3),x,x(xpoints-2)]);
for i = 1:N
    % Numerically computing the time derivative of the state
    dp        = diff([Xint(i,3),Xint(i,:),Xint(i,xpoints-2)]);
    dpdx(i,:) = (dp(1:end-1)./dx(1:end-1).*dx(2:end) ...
        + dp(2:end)./dx(2:end).*dx(1:end-1)) ...
        ./ (dx(2:end)+dx(1:end-1));
end
for i= 1:length(x)
    u = Xint(:,i);
    dudx = dpdx(:,i);
    EL = get_EL(u,dudx,Penalty);
    state  = dlarray(u);
    dstate = dlarray(dudx);
    [pLx,pLxd] = dlfeval(@AutoDiff_Dynamics,x,state,dstate,Penalty);
    
    
    Autodiff_pLx(:,i)  = pLx;
    Autodiff_pLxd(:,i) = pLxd;
    
    Analytic_pLx(:,i)  = EL(:,1);
    Analytic_pLxd(:,i) = EL(:,2);
    
end

for i = 1:12
    subplot(4,3,i)
    plot(x,Xint(i,:))
    hold on
    plot(x(1),X0(i),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',c5,...
            'MarkerFaceColor',c5,...
            'MarkerSize',8);
    plot(x(end),Xf(i),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',c6,...
            'MarkerFaceColor',c6,...
            'MarkerSize',8);
    xlabel('$\tau$')
    grid on
end

figure
for i = 1:12
    subplot(4,3,i)
    plot(x,dpdx(i,:))
    xlabel('$\tau$')
    grid on
end

figure
for i = 1:12
    subplot(4,3,i)
    plot(x,Autodiff_pLx(i,:),'b')
    hold on
    plot(x,Analytic_pLx(i,:),'r--')
    xlabel('$\tau$')
    ylabel('pLx')
    grid on
end

figure
for i = 1:12
    subplot(4,3,i)
    plot(x,Autodiff_pLxd(i,:),'b')
    hold on
    plot(x,Analytic_pLxd(i,:),'r--')
    xlabel('$\tau$')
    ylabel('pLxd')
    grid on
end




%% ----------------------- control extraction -----------------------------
disp('Extracting the required control...');
% initialize controls and states 
Crtl = struct('u', cell(1, Num_Agent)); % actual control
Dim_Crtl   = struct('u', cell(1, Num_Agent)); % actual control
xnew = linspace(0,T,xpoints_new); % time grids for integration
scale = Tc/(2*pi);
Xnom_inter = interp1(tspan,Xnom,xnew*scale);


% use spline to interpolate state from steady state pde solution, p
p      = zeros(N,xpoints_new);
dpdx   = zeros(N,xpoints_new);
dx_new = diff([xnew(3),xnew,xnew(xpoints_new-2)]);
for i = 1:N
    p(i,:)    = spline(x,sol(end,:,i),xnew);
    % Numerically computing the time derivative of the state
    dp_new    = diff([p(i,3),p(i,:),p(i,xpoints_new-2)]);
    dpdx(i,:) = (dp_new(1:end-1)./dx_new(1:end-1).*dx_new(2:end) ...
        + dp_new(2:end)./dx_new(2:end).*dx_new(1:end-1)) ...
        ./ (dx_new(2:end)+dx_new(1:end-1));
end
%
mm = 6;
for i = 1 : xpoints_new
    % f(x(t)), note: the F_d(x) in paper
    if strcmp(DynamicsModel,'THode')
        drift = NonDimentionalized_THode(xnew(i),p(:,i));
    elseif strcmp(DynamicsModel,'CWode')
        drift = NonDimentionalized_THode(xnew(i),p(:,i));
    elseif strcmp(DynamicsModel,'NLode')
        drift = NonDimentionalized_THode(xnew(i),p(:,i));
    end
    % get [Fc F], note: the F_bar matrix in paper
    Ffull = F;
    B = Ffull(:,M+1:end);
    
    OE = RVtoCOEs(Xnom_inter(i,1:3),Xnom_inter(i,4:6),mu_Earth); f = OE(6);
    h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
    K_f = 1+e_chief*cos(f);
    conts1 = h^(3/2)/(mu_Earth*K_f);
    conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
    conts3 = 1/conts1;
    A_inv = [conts1*eye(3)       zeros(3);
        conts2*eye(3)  conts3*eye(3)];
    
    % Extracting/approximate the required control
    for jj = 1:Num_Agent
        idx = 1+mm*(jj-1):mm*jj;
        Crtl(jj).u(:,i) = (B.'*B)\(B.')* ( dpdx(idx,i) - drift(idx) );
        
        Crtl_Aug = A_inv*[zeros(3,1); Crtl(jj).u(:,i)];
        Dim_Crtl(jj).u(:,i) = Crtl_Aug(4:6);
    end
end

%% plotting the results
%===========================%
%------------------- Plotting the extracted control -----------------------
if strcmp(DynamicsModel,'THode')
    xnew = linspace(0,T,intgrids);
else
    xnew = linspace(0,T*Period,intgrids);
end
figure('Name','Extracted control','Renderer', 'painters', 'Position', [10 10 900 600])
Label  = {'Agent 1', 'Agent 2', 'Agent 3', 'Agent 4'};
Label2 = {'Agent 1', 'Agent 2', 'Agent 3', 'Agent 4'};
for i = 1:Num_Agent
    subplot(2,2,i)
    plot(xnew,Crtl(i).u(1,:),'r','LineWidth',2);
    hold on
    plot(xnew,Crtl(i).u(2,:),'b','LineWidth',2);
    plot(xnew,Crtl(i).u(3,:),'g','LineWidth',2);
    ylabel(Label(i))
    % xlabel('f (rad)')
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'0','$\pi$/2','$\pi$','3$\pi$/2','2$\pi$'})
    grid on;
    
    
    subplot(2,2,i+2)
    plot(scale*xnew,Dim_Crtl(i).u(1,:),'r','LineWidth',2);
    hold on
    plot(scale*xnew,Dim_Crtl(i).u(2,:),'b','LineWidth',2);
    plot(scale*xnew,Dim_Crtl(i).u(3,:),'g','LineWidth',2);
    ylabel(Label2(i))
    xlim([scale*xnew(1) scale*xnew(end)])
    xlabel('time (s)')
    grid on;
    
    
    % if strcmp(DynamicsModel,'THode')
    %     xlabel('True Anomaly (rad)')
    % else
    %     set(gca, 'XScale', 'log')
    %     xlabel('time (s)')
    % end
    
    
end

%%
for ll = 1
    if strcmp(DynamicsModel,'THode')
        cMat1 = [c3;c5;c6;c1];
        cMat2 = [c1;c4];
        cMat = [c5;c6];
        Label_1 = {'Agent1','Agent2','Agent3','Agent4'};
        Label_2 = {'Agent5','Agent6'};
        % Label = {'Agent1','Agent2','Agent3','Agent4','Agent5','Agent6'};
        
        fh2 = figure;
        hold on
        for jj=1:Num_Agent1
            idx1 = 1+mm*(jj-1):mm*jj;
            
            plot3(X_THpos2(:,idx1(1)),X_THpos2(:,idx1(2)),X_THpos2(:,idx1(3)),'Color',cMat1(jj,:));
            plot3(X_THpos2(1,idx1(1)),X_THpos2(1,idx1(2)),X_THpos2(1,idx1(3)),'o',...
                'LineWidth',3,...
                'MarkerEdgeColor',cMat1(jj,:),...
                'MarkerFaceColor',cMat1(jj,:)',...
                'MarkerSize',8);
            hold on
            
            xlabel('X [km]')
            ylabel('Y [km]')
            zlabel('Z [km]')
            title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
            % cMat = [cMat; cMat1(jj,:)];
        end
        
        
        for jj = 1:Num_Agent2
            idx2 = 1+mm*(jj-1):mm*jj;
            plot3(X_THpos3(:,idx2(1)),X_THpos3(:,idx2(2)),X_THpos3(:,idx2(3)),'Color',cMat2(jj,:));
            plot3(X_THpos3(1,idx2(1)),X_THpos3(1,idx2(2)),X_THpos3(1,idx2(3)),'o',...
                'LineWidth',3,...
                'MarkerEdgeColor',cMat2(jj,:),...
                'MarkerFaceColor',cMat2(jj,:)',...
                'MarkerSize',8);
            % cMat = [cMat; cMat2(jj,:)];
        end
        
        plt   = zeros(Num_Agent);
        hold on
        for jj = 1:Num_Agent
            idx = 1+mm*(jj-1):mm*jj;
            txt = ['Agent_',num2str(jj)];
            plt(:,jj) = plot3(X_THpos01(1,idx(1)),X_THpos01(1,idx(2)),X_THpos01(1,idx(3)),'*',...
                'LineWidth',3,...
                'MarkerEdgeColor',cMat(jj,:),...
                'MarkerFaceColor',cMat(jj,:)',...
                'MarkerSize',8);
            
            xlabel('X [km]')
            ylabel('Y [km]')
            zlabel('Z [km]')
            title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
            
            grid on
            h5 = plot3(0,0,0,'bo',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',15);
            h4 = plot3(X_THpos01(:,idx(1)),X_THpos01(:,idx(2)),X_THpos01(:,idx(3)),'Color',cMat(jj,:),'DisplayName','Initial Traj');
            
        end
        % r = 1;
        % [x_sphere, y_sphere, z_sphere] = sphere(100);
        % h = surf(r*x_sphere, r*y_sphere, r*z_sphere);
        % set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
        % axis equal
        
        n1 = 2;
        n2 = 400;
        
        hL = legend([plt(end,1), plt(end,2), h5,],...
            {'Agent1','Agent2','Chief'},'AutoUpdate','off');
        % Programatically move the Legend
        newPosition = [.19 .68 0.01 0.01];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
        view(31,18)
        
        for k = 1:Num_Agent
            idx = 1+mm*(k-1):mm*k;
            txt = ['Traj Agent_',num2str(k)];
            plot3(X_ode45_LVLH(idx(1),n1:n2),X_ode45_LVLH(idx(2),n1:n2),X_ode45_LVLH(idx(3),n1:n2),...
                '--','Color',cMat(k,:),'LineWidth',2,'DisplayName',txt);
        end
        view(27,12)
    end
end

%%
for j = 1
    
    % Crt = zeros(3,1);
    % [~, X_rel_0] = ode113(@(t,x)CW_ode(t,x,Crt,mu_dot),tspan,X01,options);
    % [~, X_nom_1] = ode113(@(t,X)CW_ode(t,X,Crt,mu_dot),tspan,Xo,options);
    % Xo_1 = nan(mm,Num_Agent1);
    % for i = 1:Num_Agent1
    %     idx = find(time <= i*IntTime/(2*Num_Agent1), 1,'last');
    %     Xrel = [Rrho\X_nom_1(idx,1:3)'; (Rrho\X_nom_1(idx,4:6)')*Tc];
    %     Xo_1(:,i) = Xrel;
    % end
    % [~, X_rel_1] = ode113(@(t,x)CW_ode(t,x,Crt,mu_dot),tspan,Xo_1,options);
    %
    % [~, X_nom_2] = ode113(@(t,X)CW_ode(t,X,Crt,mu_dot),tspan,Xo,options);
    % Xo_2 = nan(mm,Num_Agent2);
    % for i = 1:Num_Agent2
    %     idx = find(time <= (2*i-1)*IntTime/(2*Num_Agent2), 1,'last');
    %     Xrel = [Rrho\X_nom_2(idx,1:3)'; (Rrho\X_nom_2(idx,4:6)')*Tc];
    %     Xo_2(:,i) = Xrel;
    % end
    % [~, X_rel_2] = ode113(@(t,x)CW_ode(t,x,Crt,mu_dot),tspan,Xo_2,options);
    
    
    

 

    % % %%
    % % close all
    % fh2 = figure;
    % plt0 = plot3(0,0,0,'bo',...
    %     'LineWidth',2,...
    %     'MarkerEdgeColor','k',...
    %     'MarkerFaceColor','k',...
    %     'MarkerSize',15);
    % hold on
    % %plot3(X_CW(:,1),X_CW(:,2),X_CW(:,3),'m')
    % plot3(X_THpos0(:,1),X_THpos0(:,2),X_THpos0(:,3),'Color',c5,'DisplayName','Agnet 1 Init Traj')
    % plot3(X_THpos1(:,1),X_THpos1(:,2),X_THpos1(:,3),'Color',c6,'DisplayName','Agnet 2 Init Traj')
    % plot3(X_THpos2(:,1),X_THpos2(:,2),X_THpos2(:,3),'Color',c3,'DisplayName','Agnet 1 Final Traj')
    % plot3(X_THpos3(:,1),X_THpos3(:,2),X_THpos3(:,3),'Color',c1,'DisplayName','Agnet 2 Final Traj')
    % % plot3(X_NL2(:,1),X_NL2(:,2),X_NL2(:,3),'m-.')
    % plt1 = plot3(X01(1),X01(2),X01(3),'o',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c5,...
    %     'MarkerFaceColor',c5,...
    %     'MarkerSize',8);
    % plt2 = plot3(X02(1),X02(2),X02(3),'o',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c6,...
    %     'MarkerFaceColor',c6,...
    %     'MarkerSize',8);
    % plt3 = plot3(Xf1(1),Xf1(2),Xf1(3),'o',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c3,...
    %     'MarkerFaceColor',c3,...
    %     'MarkerSize',8);
    % plt4 = plot3(Xf2(1),Xf2(2),Xf2(3),'o',...
    %     'LineWidth',3,...
    %     'MarkerEdgeColor',c1,...
    %     'MarkerFaceColor',c1,...
    %     'MarkerSize',8);
    % grid on
    % view(68,10)
    % hL = legend([plt0, plt1, plt2, plt3,plt4],...
    %     {'Chief Location','Agent1 Init Cond','Agent2 Init Cond','Agent1 End Cond','Agent2 End Cond'},'AutoUpdate','off');
    % hL.FontSize = 20;
    % newPosition = [0.21 0.72 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1  );
    % % set all units inside figure to normalized so that everything is scaling accordingly
    %  set(findall(fh2,'Units','pixels'),'Units','normalized');
    %   % do not show figure on screen
    %  set(fh2, 'visible', 'on')
    %  % set figure units to pixels & adjust figure size
    %  fh2.Units = 'pixels';
    %  fh2.OuterPosition = [10 10 900 800];
    %  % define resolution figure to be saved in dpi
    %  res = 500;
    %  % recalculate figure size to be saved
    %  set(fh2,'PaperPositionMode','manual')
    %  fh2.PaperUnits = 'inches';
    %  fh2.PaperPosition = [0 0 5580 4320]/res;
    %  % save figure
    %  print(fh2,'InitialConditions','-dpng',sprintf('-r%d',res))
    
    % %%
    % figure
    % for i=1:3
    % subplot(3,1,i)
    % plot(T_normalized,X_CW(:,i),'m')
    % % plot(T_normalized,X_THpos2(:,i),'k--')
    % hold on
    % plot(T_normalized,X_THpos(:,i),'r--')
    % plot(T_normalized,X_NL(:,i),'b-.')
    % % plot(T_normalized,X_NL2(:,i),'m-.')
    % end
end


%% --------------------------- AGHF parameters ----------------------------
for ll = 1
    smax = 9e7; % Need tpo find an analytic way to solve for this value
    % # of grids in t
    tgrids = 250;
    % # of grids in s
    sgrids = 250;
    % # of grids of integration for final trajectory extraction
    intgrids = 1e4;
    % # of inputs
    M = 3;
    % motion duration
    T = SimTime;
    % penalty value (the lambda in paper)
    Penalty = 2e4;
    % tolerance of the integrator
    % opts = odeset('AbsTol',1e-14);
    opts = odeset('RelTol',2.22045e-8,'AbsTol',2.22045e-20);
    
    % Setting the boundary condtions and the intial condition
    tmax = smax; xpoints = tgrids; xpoints_new = intgrids; tpoints = sgrids;
    m = 0; t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of s interval in log scale
    
    if strcmp(DynamicsModel,'THode')
        %############## Initial boundary conditions #################
        f = tspan1(1);
        h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
        K_f = 1+e_chief*cos(f);
        conts1 = mu_Earth*K_f/(h^(3/2));
        conts2 = -mu_Earth*e_chief*sin(f)/(h^(3/2));
        conts3 = 1/conts1;
        A_map1 = [conts1*eye(3)       zeros(3);
            conts2*eye(3)  conts3*eye(3)];
        X01 = A_map1*X01;
        X02 = A_map1*X02;
        
        %############## Final boundary conditions #################
        f = tspan2(1);
        h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
        K_f = 1+e_chief*cos(f);
        conts1 = mu_Earth*K_f/(h^(3/2));
        conts2 = -mu_Earth*e_chief*sin(f)/(h^(3/2));
        conts3 = 1/conts1;
        A_map2 = [conts1*eye(3)       zeros(3);
            conts2*eye(3)  conts3*eye(3)];
        Xf1 = A_map2*Xf1;
        Xf2 = A_map2*Xf2;
        
        %############## Final boundary conditions #################
        T    = PNum*wrapTo2Pi(f2)-f1;
        x    = linspace(f1,PNum*wrapTo2Pi(f2),xpoints); % discretization of the curve in time
        xnew = linspace(f1,PNum*wrapTo2Pi(f2),intgrids);
        
    elseif strcmp(DynamicsModel,'CWode')
        X01  = X01/Tc;
        X02  = X02/Tc;
        Xf1  = Xf1/Tc;
        Xf2  = Xf2/Tc;
        x    = linspace(0,T,xpoints); % discretization of the curve in time
        xnew = linspace(0,T*Period,intgrids);
        
    elseif strcmp(DynamicsModel,'NLode')
        X01  = [Rrho\X01(1:3); (Rrho\X01(4:6))*Tc];
        X02  = [Rrho\X02(1:3); (Rrho\X02(4:6))*Tc];
        Xf1  = [Rrho\Xf1(1:3); (Rrho\Xf1(4:6))*Tc];
        Xf2  = [Rrho\Xf2(1:3); (Rrho\Xf2(4:6))*Tc];
        x    = linspace(0,T,xpoints);
        xnew = linspace(0,T*Period,intgrids);
        
    end
    
    X0 = [X01; X02]; X0 = X0int; %X0(:);
    Xf = [Xf1; Xf2]; Xf = X0end; %Xf(:);
    N  = length(X0);
end


% % Specify the initial relative-orbital element of the deputies
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% dela = 0;
% dele = 1/(8*a_chief);
% deli = 1/(3*a_chief);
% delLitOmg = 0;
% delBigOmg = 0;
% delM = 0;
% 
% delq1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
% delq2 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
% deltheta = delLitOmg + delM; % Relative true latitude in rad
% delCOE = [dela, deltheta, deli, delq1, delq2, delBigOmg]';
% % Defining the deputies' initial trajectory
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% OE_chief = [a_chief, theta_chief, inc_chief, q1_chief, q2_chief, BigOmg_chief];
% AMap = ForwardMapping(OE_chief, mu_Earth); % Linear mapping matrix
% Xo = AMap*delCOE;
% u = zeros(3,1);
% [time, X_nom] = ode113(@(t,X)CW_ode(t,X,u,mu_dot),tspan,Xo,options);
% Xo = nan(mm,Num_Agent);
% for i = 1:Num_Agent
%     idx = find(time <= i*IntTime/Num_Agent, 1,'last');
%     Xo(:,i) = X_nom(idx,:)';
% end



% generate functions for EL, metric G, drift term f, control matrix F
% please modify this function if system is changed
[get_EL, get_G, get_f, get_F] = generate_fh_of_model(Num_Agent); 


dX0 = rand(12,1); k = 200*rand(1); x = rand(1);
XChief = full([Chief_x(x);  Chief_y(x);  Chief_z(x);...
               Chief_dx(x); Chief_dy(x); Chief_dz(x)]);
LL = full(L(X0,dX0,k,x,XChief));
CasADiresult = full(dLdx(X0,dX0,k,x,XChief,LL))';
CasADir = [CasADiresult(1:12), CasADiresult(13:24)];



drift_CasADi = full(fdrift(x,X0,XChief));


state  = dlarray(X0);
dstate = dlarray(dX0);
[pLx2,pLxd2] = dlfeval(@AutoDiff_Dynamics,state,dstate,k,x);

test = CasADir - [pLx2,pLxd2]





%%
X_rel_0 = []; X_rel_1 = [];
for i = 1:2
    X0 = [X01,X02];
    Xf = [Xf1,Xf2];
    [~, X_rel_01] = ode113(@(t,x)NonDimentionalized_NLode(t,x),T_normalized,X0(:,i),options);
    [~, X_rel_11] = ode113(@(t,x)NonDimentionalized_NLode2(t,x),T_normalized,Xf(:,i),options);
    X_rel_0   = [X_rel_0, X_rel_01];  X_rel_1 = [X_rel_1, X_rel_11];
end
cMat1 = [c5;c6];
cMat2 = [c3;c1];
cMat0 = [c7;c4];
for jj = 1:Num_Agent
        idx = 1+mm*(jj-1):mm*jj;
        plot3(X_rel_0(:,idx(1)),X_rel_0(:,idx(2)),X_rel_0(:,idx(3)),'Color',cMat1(jj,:),'DisplayName','Initial Traj');
        hold on
        plot3(X_rel_1(:,idx(1)),X_rel_1(:,idx(2)),X_rel_1(:,idx(3)),'Color',cMat2(jj,:));
         plot3(X_rel_0(1,idx(1)),X_rel_0(1,idx(2)),X_rel_0(1,idx(3)),'*',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat1(jj,:),...
            'MarkerFaceColor',cMat1(jj,:)',...
            'MarkerSize',8);
        
        
        plot3(X_rel_1(1,idx(1)),X_rel_1(1,idx(2)),X_rel_1(1,idx(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat2(jj,:),...
            'MarkerFaceColor',cMat2(jj,:),...
            'MarkerSize',8);
        
        plot3(0,0,0,'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','k',...
            'MarkerSize',15);
        grid on
        
end

%%
for ll = 1
    disp('Extracting the required control...');
    % initialize controls and states
    Crtl = struct('u', cell(1, Num_Agent)); % actual control
    % use spline to interpolate state from steady state pde solution, p
    Pint   = zeros(N,xpoints_new);
    p      = zeros(N,xpoints_new);
    dpdx   = zeros(N,xpoints_new);
    dx_new = diff([xnew(3),xnew,xnew(xpoints_new-2)]);
    for i = 1:N
        Pint(i,:) = spline(x,sol(1,:,i),xnew);
        p(i,:)    = spline(x,sol(end,:,i),xnew);
        % Numerically computing the time derivative of the state
        dp_new    = diff([p(i,3),p(i,:),p(i,xpoints_new-2)]);
        dpdx(i,:) = (dp_new(1:end-1)./dx_new(1:end-1).*dx_new(2:end) ...
            + dp_new(2:end)./dx_new(2:end).*dx_new(1:end-1)) ...
            ./ (dx_new(2:end)+dx_new(1:end-1));
    end
    for i = 1 : length(xnew)
        if strcmp(DynamicsModel,'THode')
            f      = xnew(i);
            h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f    = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv  = [conts1*eye(3)       zeros(3);
                conts2*eye(3)  conts3*eye(3)];
            
            drift = NonDimentionalized_THode(xnew(i),p(:,i));
        elseif strcmp(DynamicsModel,'CWode')
            drift = NonDimentionalized_CWode(xnew(i),p(:,i));
        elseif strcmp(DynamicsModel,'NLode')
            d1 = NonDimentionalized_NLode(xnew(i),p(1:6,i));
            d2 = NonDimentionalized_NLode(xnew(i),p(7:12,i));
            drift = [d1;d2];
        end
        % get [Fc F], note: the F_bar matrix in paper
        Ffull = F;
        B = Ffull(:,M+1:end);
        % Extracting/approximate the required control
        for jj = 1:Num_Agent
            idx = 1+mm*(jj-1):mm*jj;
            Crtl(jj).u(:,i) = (B.'*B)\(B.') * ( dpdx(idx,i) - drift(idx) );
            if strcmp(DynamicsModel,'THode')
                Crtl_Aug = A_inv*[zeros(3,1); Crtl(jj).u(:,i)];
                Dim_Crtl(jj).u(:,i) = Crtl_Aug(4:6);
            elseif strcmp(DynamicsModel,'CWode')
                Dim_Crtl(jj).u(:,i) = Tc*Crtl(jj).u(:,i);
            elseif strcmp(DynamicsModel,'NLode')
                Dim_Crtl(jj).u(:,i) = Rrho*Crtl(jj).u(:,i);
            end
        end
    end
    %%
    if exist('u11','file') == 3
        delete u11.c u11.mexmaci64
    end
    if exist('u12','file') == 3
        delete u12.c u12.mexmaci64
    end
    if exist('u13','file') == 3
        delete u13.c u13.mexmaci64
    end
    if exist('u21','file') == 3
        delete u21.c u21.mexmaci64
    end
    if exist('u22','file') == 3
        delete u22.c u22.mexmaci64
    end
    if exist('u23','file') == 3
        delete u23.c u23.mexmaci64
    end
    
    u11 = casadi.interpolant('U','bspline',{xnew}, Crtl(1).u(1,:));
    u12 = casadi.interpolant('U','bspline',{xnew}, Crtl(1).u(2,:));
    u13 = casadi.interpolant('U','bspline',{xnew}, Crtl(1).u(3,:));
    u21 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(1,:));
    u22 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(2,:));
    u23 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(3,:));
    
    u11.generate('u11.c',CasADiopts);
    u12.generate('u12.c',CasADiopts);
    u13.generate('u13.c',CasADiopts);
    u21.generate('u21.c',CasADiopts);
    u22.generate('u22.c',CasADiopts);
    u23.generate('u23.c',CasADiopts);
    mex u11.c -largeArrayDims
    mex u12.c -largeArrayDims
    mex u13.c -largeArrayDims
    mex u21.c -largeArrayDims
    mex u22.c -largeArrayDims
    mex u23.c -largeArrayDims
    clear u11 u12 u13 u21 u22 u23
end
% -------------------------------------------------------------------------
for i = 1 : length(xnew)
    u1test = [u11(i) u12(i) u13(i)]';
end
plot(xnew,u1test)

%% ------------------ Integrate the system's trajectory -------------------
for ll = 1
    disp('Integrating the resulting trajectory...');
    tolerance = 1e-17;
    % dt = 1e-16; dt1 = 1e-16; dt2 = 1e-16;
    % X_ode45(:,1) = X0;
    % X_ode45_LVLH(:,1) = [Xi1;Xi2];
    % Pint(:,1)         = [Xi1;Xi2];
    
    if strcmp(DynamicsModel,'THode')
        [~,X_ode45] = ode113(@(t,x)TH_ode(t,x),xnew,X0,options);
    elseif strcmp(DynamicsModel,'CWode')
        [~,X_ode45] = ode113(@(t,x)CW_ode(t,x),xnew,X0,options);
    elseif strcmp(DynamicsModel,'NLode')
        [~,X_ode45_1] = ode113(@(t,x)NL_ode(t,x,1),xnew,Xi1,options);
        [~,X_ode45_2] = ode113(@(t,x)NL_ode(t,x,2),xnew,Xi2,options);
        X_ode45 = [X_ode45_1 X_ode45_2];
    end
    
    % Crt = nan(M,Num_Agent);
    
    % for jj = 1:Num_Agent
    %     U1 = [u11(i) u12(i) u13(i)];
    %     U2 = [u21(i) u22(i) u23(i)];
    %     Crt(:,jj) = [U1 U2]; %Crtl(jj).u(:,i);
    % end
    % if strcmp(DynamicsModel,'THode')
    %     [y_f, ~, dt] = Runge_Kutta_Fehlberg_4_5(@(t,X)TH_ode(t,X),X_ode45(:,i),xnew(i),dt,xnew(i+1),tolerance);
    %     X_ode45(:,i+1) = y_f;
    % elseif strcmp(DynamicsModel,'CWode')
    %     [y_f, ~, dt] = Runge_Kutta_Fehlberg_4_5(@(t,X)CW_ode(t,X),X_ode45(:,i),xnew(i),dt,xnew(i+1),tolerance);
    %     X_ode45(:,i+1) = y_f;
    % elseif strcmp(DynamicsModel,'NLode')
    %     [y_f1, ~, dt1] = Runge_Kutta_Fehlberg_4_5(@(t,X)NL_ode(t,X),X_ode45(1:6,i),xnew(i),dt1,xnew(i+1),tolerance);
    %     [y_f2, ~, dt2] = Runge_Kutta_Fehlberg_4_5(@(t,X)NL_ode(t,X),X_ode45(7:12,i),xnew(i),dt2,xnew(i+1),tolerance);
    %     X_ode45(:,i+1) = [y_f1;y_f2];
    %
    % end
    
    % Crt = nan(M,Num_Agent);
    
    if strcmp(DynamicsModel,'THode')
        for i = 1:length(xnew)-1
            f      = xnew(i+1);
            h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f    = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv  = [conts1*eye(3)       zeros(3);
                conts2*eye(3)  conts3*eye(3)];
            X_ode45_LVLH(:,i+1) = (blkdiag(A_inv,A_inv)*X_ode45_1(i+1,:).').';
            Pint(:,i+1)         = blkdiag(A_inv,A_inv)*Pint(:,i+1);
        end
    elseif strcmp(DynamicsModel,'CWode')
        X_ode45_LVLH = Tc*X_ode45;
        Pint         = Tc*Pint;
    elseif strcmp(DynamicsModel,'NLode')
        X_ode45_LVLH = [(blkdiag(Rrho,Rrho)*X_ode45_1.').'; (blkdiag(Rrho,Rrho)*X_ode45_2.').'];
        Pint         = [blkdiag(Rrho,Rrho)*Pint(1:6,:); blkdiag(Rrho,Rrho)*Pint(7:12,:)];
    end
    
    disp('Done!!!!!');
end
% -------------------------------------------------------------------------

%%
% X0 = [X01 X02];
% Xf = [Xf1 Xf2];
% tic
% options = odeset('RelTol',2.22045e-9,'AbsTol',2.22045e-30);
% X_THpos01 = []; X_THpos11 = []; X_rel_0 = []; X_rel_1 = [];
% for i = 1:Num_Agent
%     X0integral = [X01,X02];
%     Xfintegral = [Xf1,Xf2];
%     [~, X_rel_01] = ode113(@(t,x)NonDimentionalized_FlatPlateOde(t,x),T_normalized,X0(:,i),options);
%     [~, X_rel_11] = ode113(@(t,x)NonDimentionalized_FlatPlateOde2(t,x),T_normalized,Xf(:,i),options);
%     X_rel_0 = [X_rel_0, X_rel_01];
%     X_rel_1 = [X_rel_1, X_rel_11];
% end
% toc
% X0 = X0(:); Xf = Xf(:);




%%
function [Xdot] = NonDimentionalized_NLode1(t,X,Target)
global mu_Earth Rrho Tc 
format long
ti = t*Tc;
X_chief = full([Chief_x0(t);  Chief_y0(t);  Chief_z0(t);...
                Chief_dx0(t); Chief_dy0(t); Chief_dz0(t)]);
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
r_target = X_chief(1:3);  v_target = X_chief(4:6); rt_norm = (r_target.'*r_target).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = Rrho*rho_bar; rho_prime = Rrho*drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, V_sun, beta_chief] = F_CanonBall(ti,r_target,Target); % SRP force on the Target
% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(r_target)*v_target; h_norm = norm(h_vec);  u_acc = U_eci_Target;
eh = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
V_rel = V_sun - v_target;
u_acc_dot = beta_chief*(V_rel/r_rel_norm^3 - 3*dot(Xrel,V_rel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot = (blkdiag(Rrho,Rrho)\Xdot)*Tc;

end












