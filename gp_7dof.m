%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 7-DOF DVP System : GP %%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Deterministic and Stochastic Excitation %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

f = figure('Name','GP Plots'); set(gcf,'Position',[1000,200,900,500]);

%%%%% Saved Workspace From File 'saved_dp_seven_dof.m'
load('workspace_saved_dps_7dof.mat')

t_upper = 10000; time = 0:50:t_upper;

for s = 1:6
    points = 150; y = param_est(s,(1:points))';
    x = time(1,(1:points))'; t_str = 0:50:t_upper;
    
    mdl = fitrgp(x,y,'basis','pureQuadratic',...
        'KernelFunction','squaredexponential');
    [a,~,b] = predict(mdl,t_str');
    
    if s==1
        figure(f); subplot(2,3,s); hold on; plot(time,param_tr(s,(1:length(param_tr))),'r');
        patch([t_str';flipud(t_str')],[b(:,1);flipud(b(:,2))],'c','FaceAlpha',0.1);
        plot(t_str,a,'b'); xline(time(points)); legend('True','95% CI','GP Estimated','Last Training Data');
        xlabel('Time (days)'); ylabel(['k',num2str(s)]);
    elseif s==2 || s==3
        figure(f); subplot(2,3,s); hold on; plot(time,param_tr(s,(1:length(param_tr))),'r');
        patch([t_str';flipud(t_str')],[b(:,1);flipud(b(:,2))],'c','FaceAlpha',0.1);
        plot(t_str,a,'b'); xline(time(points)); xlabel('Time (days)'); ylabel(['k',num2str(s)]);
    else
        figure(f); subplot(2,3,s); hold on; plot(time,param_tr(s,(1:length(param_tr))),'r');
        patch([t_str';flipud(t_str')],[b(:,1);flipud(b(:,2))],'c','FaceAlpha',0.1);
        plot(t_str,a,'b'); xline(time(points)); xlabel('Time (days)'); ylabel(['k',num2str(s+1)]);
    end
end

