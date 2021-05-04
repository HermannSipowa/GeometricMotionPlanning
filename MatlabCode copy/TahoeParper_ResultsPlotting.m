%% Importing the data
dataPath = 'Outputs_Data/';
for i = 1
% data = sprintf('State%d', i);
% LoadedData = load([dataPath,data,'.mat']);
% solMat(:,:,i) = ;
end
LoadedData = load([dataPath,'State1.mat']);  solMat(:,:,1)  = LoadedData.State1;
LoadedData = load([dataPath,'State2.mat']);  solMat(:,:,2)  = LoadedData.State2;
LoadedData = load([dataPath,'State3.mat']);  solMat(:,:,3)  = LoadedData.State3;
LoadedData = load([dataPath,'State4.mat']);  solMat(:,:,4)  = LoadedData.State4;
LoadedData = load([dataPath,'State5.mat']);  solMat(:,:,5)  = LoadedData.State5;
LoadedData = load([dataPath,'State6.mat']);  solMat(:,:,6)  = LoadedData.State6;
LoadedData = load([dataPath,'State7.mat']);  solMat(:,:,7)  = LoadedData.State7;
LoadedData = load([dataPath,'State8.mat']);  solMat(:,:,8)  = LoadedData.State8;
LoadedData = load([dataPath,'State9.mat']);  solMat(:,:,9)  = LoadedData.State9;
LoadedData = load([dataPath,'State10.mat']); solMat(:,:,10) = LoadedData.State10;
LoadedData = load([dataPath,'State11.mat']); solMat(:,:,11) = LoadedData.State11;
LoadedData = load([dataPath,'State12.mat']); solMat(:,:,12) = LoadedData.State12;
LoadedData = load([dataPath,'State13.mat']); solMat(:,:,13) = LoadedData.State13;
XUn = load([dataPath,'XUnctrl','.mat']); XUncontrolled = XUn.XUnctrl;
Xcontrol = load([dataPath,'Xstar','.mat']); Xstar_f = Xcontrol.Xstar;
% tData = load([dataPath,'tspanData','.mat']); tspanData = tData.tspanData;
Ctrl = load([dataPath,'Ctrl','.mat']); Tau_ctrl = Ctrl.Ctrl;

%% Visualizing the results
%-------------------------------------------
 for k = 1
    close all
    clc
    start_up
    c1 = rgb('Black'); c2 = rgb('Cyan'); c3 = rgb('Goldenrod');
    c4 = rgb('Lime'); c5 = rgb('DarkBlue'); c6 = rgb('DarkTurquoise');
    c7 = rgb('Crimson'); c8 = rgb('DarkSlateGray');
    Size_array = floor(size(solMat(:,:,1),1)/10);
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~% sol(end,:,i)
    % Plotting the 3D homotopies
    figure('Renderer', 'painters', 'Position', [15 15 900 600])
    h1 = plot3(solMat(:,8),solMat(1,:,9),solMat(1,:,10),'r');
    hold on
    h2 = plot3(solMat(1,1,8),solMat(1,1,9),solMat(1,1,10),'o','color', c2,...
        'LineWidth',2,...
        'MarkerEdgeColor',c2,...
        'MarkerFaceColor',c2,...
        'MarkerSize',5);
    h3 = plot3(solMat(1,end,8),solMat(1,end,9),solMat(1,end,10),'go',...
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
    grid on
    for i=1:10
        idx= Size_array*i;
        X = solMat(idx,:,8)  ;h1.XDataSource='X';
        Y = solMat(idx,:,9)  ;h1.YDataSource='Y';
        Z = solMat(idx,:,10) ;h1.ZDataSource='Z';
        
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
            plot1 = plot(tspanData/86400,DQ1(:,i),'color', c4);
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
            plot(tspanData/86400,Tau_ctrl(i,:),'color', c6)
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
            plot(tspanData/86400,Q1(:,i),'color', c5) % -Q2(:,i)
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
            plot(tspanData/86400,Q2(:,i),'color', c5)
            ylabel(YLabel(i))
            xlabel('days')
            grid on
            % legend('Classical Dynamics', 'Dual Quaternions')
        end
        sgt = sgtitle('Deputy Angular Velocity','interpreter','tex');
        sgt.FontSize = 30;
    end
    
end % Plots for relative motions (Chaser and Target) with classical dynamics
