%%%%%%%%%%%%%%%%%%%%%%% Shailesh Garg %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 2-DOF DO System : GP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Deterministic and Stochastic Excitation %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

f = figure; set(gcf,'Position',[1000,400,900,350]);

%%%%% Saved Workspace From File 'saved_dp_two_dof.m'
load('workspace_saved_dps_2dof.mat')

t_upper = 10000; time = 0:50:t_upper;

points = 150; y = K1_est((1:points),1); x = time(1,(1:points))';

mdl1 = fitrgp(x,y,'Basis','pureQuadratic');

[a,~,b] = predict(mdl1,time');

figure(f); subplot(1,2,1); plot(time,K1,'r'); hold on;
patch([time';flipud(time')],[b(:,1);flipud(b(:,2))],'c','FaceAlpha',0.1);
plot(time,a,'b'); xline(time(points));

legend({'True','95% CI','GP Estimated','Last Training Data'},'FontSize',14,'FontWeight','bold')
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold','FontName','Times');
set(gcf,'Position',[1000,200,900,400]);
set(get(gca,'xlabel'),'String','Time (days)','FontSize',20,'FontWeight','bold','FontName','Times','Interpreter','tex');
set(get(gca,'ylabel'),'String','k_1','FontSize',20,'FontWeight','bold','FontName','Times','Interpreter','tex');
set(gcf,'color','w'); box on;

y = K2_est((1:points),1); x = time(1,(1:points))';

mdl2 = fitrgp(x,y,'Basis','pureQuadratic');

[a,~,b] = predict(mdl2,time');

subplot(1,2,2); plot(time,K2,'r'); hold on;
patch([time';flipud(time')],[b(:,1);flipud(b(:,2))],'c','FaceAlpha',0.1);
plot(time,a,'b'); xline(time(points));

set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold','FontName','Times');
set(gcf,'Position',[1000,200,900,400]);
set(get(gca,'xlabel'),'String','Time (days)','FontSize',20,'FontWeight','bold','FontName','Times','Interpreter','tex');
set(get(gca,'ylabel'),'String','k_2','FontSize',20,'FontWeight','bold','FontName','Times','Interpreter','tex');
set(gcf,'color','w'); box on;
