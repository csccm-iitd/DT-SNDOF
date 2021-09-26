 %%%%%%%%%% Data Points Generated For 7DOF System %%%%%%%%%%
%%%%%%%%%% Subjected To Deterministic And Stochastic Forces  %%%%%%%%%%
%%%%%%%%%% Shailesh Garg  %%%%%%%%%%

clc
clear
close all

%%%%%%%%%% Time Vector For Individual Data Point %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dt
dt = 0.001; t = 0:dt:5;

%%%%%%%%%% Stiffness And Damping Simulation %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m1 = 20; m2 = 20; m3 = 10; m4 = 10; m5 = 10; m6 = 10; m7 = 5; alpha = 100;
c1 = 20; c2 = 20; c3 = 20; c4 = 20; c5 = 20; c6 = 20; c7 = 20;
k1 = 2000; k2 = 2000; k3 = 1000; k4 = 1000; k5 = 1000; k6 = 1000; k7 = 500;

t_upper = 10000; time = 0:50:t_upper; degradation = exp(-1*0.00005*time);
K = [k1*degradation; k2*degradation; k3*degradation; k4*ones(1,length(time))
    k5*degradation; k6*degradation; k7*degradation];
C = [c1*degradation; c2*degradation; c3*degradation; c4*degradation;
    c5*degradation; c6*degradation; c7*degradation];

%%%%%%%%%% Data Simulation Using Taylor 1.5 Strong %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K1_est = zeros(length(time),1); K2_est = zeros(length(time),1); K3_est = zeros(length(time),1);
K4_est = zeros(length(time),1); K5_est = zeros(length(time),1);
K6_est = zeros(length(time),1); K7_est = zeros(length(time),1);

nsim = 1; %%%%%%%%%% No. Of Simulations For A Particular Data Point %%%%%%%%%%

%%%%%%%%%% Iterations %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dps_pointer = 1; data =  zeros(length(t),7*length(time)); amp = 10; omega = 10;
F = [0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t;...
    amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t)]; % 14x1

for dp = 1:length(time)
    acceleration = zeros(7,length(t));
    
    for sim = 1:nsim
        y = zeros(14,length(t)); acc_sim = zeros(7,length(t)); fprintf('Iteration %d Simulation %d\n',dp,sim);
        k1 = K(1,dp); k2 = K(2,dp); k3 = K(3,dp); k4 = K(4,dp); k5 = K(5,dp); k6 = K(6,dp); k7 = K(7,dp);
        c1 = C(1,dp); c2 = C(2,dp); c3 = C(3,dp); c4 = C(4,dp); c5 = C(5,dp); c6 = C(6,dp); c7 = C(7,dp);
        
        for i = 1:length(t)-1
            y(:,i+1) = disp_vel_IT(y(:,i), F(:,i), K(:,dp), C(:,dp));
        end
        
        for i = 1:length(t)
            y1 = y(1,i); y2 = y(2,i); y3 = y(3,i); y4 = y(4,i); y5 = y(5,i); y6 = y(6,i); y7 = y(7,i);...
                y8 = y(8,i); y9 = y(9,i); y10 = y(10,i); y11 = y(11,i); y12 = y(12,i); y13 = y(13,i); y14 = y(14,i);
            
            acc_sim(:,i) = [-(y1*(k1+k2)-c2*y4-k2*y3+y2*(c1+c2))/m1;
                (c2*y2-y3*(k2+k3)+c3*y6+k2*y1+k3*y5-y4*(c2+c3))/m2;
                -(k4*y7-c4*y8-k3*y3-c3*y4+y5*(k3-k4)+alpha*(y5-y7)^3+y6*(c3+c4))/m3;
                (c4*y6+c5*y10-k4*y5+k5*y9+y7*(k4-k5)+alpha*(y5-y7)^3-y8*(c4+c5))/m4;
                (c5*y8-y9*(k5+k6)+c6*y12+k5*y7+k6*y11-y10*(c5+c6))/m5;
                (c6*y10-y11*(k6+k7)+c7*y14+k6*y9+k7*y13-y12*(c6+c7))/m6;
                (c7*y12-c7*y14+k7*y11-k7*y13)/m7];
        end
        acceleration = acceleration+acc_sim;
    end
    
    acc_dps{data_dps_pointer} = acceleration/nsim;
    y_dps{data_dps_pointer} = y;
    kgt_dps{data_dps_pointer} = [k1; k2; k3; k4; k5; k6; k7];
    data_dps_pointer = data_dps_pointer+1;

    data(:,dp) = acceleration(1,:)/nsim; data(:,(length(time)+dp)) = acceleration(2,:)/nsim;
    data(:,(2*length(time)+dp)) = acceleration(3,:)/nsim; data(:,(3*length(time)+dp)) = acceleration(4,:)/nsim;
    data(:,(4*length(time)+dp)) = acceleration(5,:)/nsim; data(:,(5*length(time)+dp)) = acceleration(6,:)/nsim;
    data(:,(6*length(time)+dp)) = acceleration(7,:)/nsim;
end

%%%%%%%%%% Measurement Simulation %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr = 50; measurement =  zeros(length(t),2*length(time));

for dp = 1:length(time)
    measurement(:,dp) = data(:,dp)+sqrt((1/snr)*(std(data(:,dp))^2))*randn(length(data(:,dp)),1);
    
    measurement(:,(length(time)+dp)) = data(:,(length(time)+dp))...
        +sqrt((1/snr)*(std(data(:,(length(time)+dp)))^2))*randn(length(data(:,(length(time)+dp))),1);
    
    measurement(:,(2*length(time)+dp)) = data(:,(2*length(time)+dp))...
        +sqrt((1/snr)*(std(data(:,(2*length(time)+dp)))^2))*randn(length(data(:,(2*length(time)+dp))),1);
    
    measurement(:,(3*length(time)+dp)) = data(:,(3*length(time)+dp))...
        +sqrt((1/snr)*(std(data(:,(3*length(time)+dp)))^2))*randn(length(data(:,(3*length(time)+dp))),1);
    
    measurement(:,(4*length(time)+dp)) = data(:,(4*length(time)+dp))...
        +sqrt((1/snr)*(std(data(:,(4*length(time)+dp)))^2))*randn(length(data(:,(4*length(time)+dp))),1);
    
    measurement(:,(5*length(time)+dp)) = data(:,(5*length(time)+dp))...
        +sqrt((1/snr)*(std(data(:,(5*length(time)+dp)))^2))*randn(length(data(:,(5*length(time)+dp))),1);
    
    measurement(:,(6*length(time)+dp)) = data(:,(6*length(time)+dp))...
        +sqrt((1/snr)*(std(data(:,(6*length(time)+dp)))^2))*randn(length(data(:,(6*length(time)+dp))),1);
end

%%%%%%%%%% UKF Application Using Euler Maruyama %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dps_pointer = 1;

%%%%%%%%%% Weights Calculations %%%%%%%%%%

m0 = zeros(21,1); n = length(m0); al = 0.001; beta = 2; kappa = 0;
la = al^2*(n+kappa)-n; wm = zeros(2*n+1,1); wc = zeros(2*n+1,1);

for j=1:2*n+1
    if j==1
        wm(j) = la/(n+la);
        wc(j) = la/(n+la)+(1-al^2+beta);
    else
        wm(j) =1/(2*(n+la));
        wc(j) = wm(j);
    end
end

m1 = 20; m2 = 20; m3 = 10; m4 = 10; m5 = 10; m6 = 10; m7 = 5; alpha = 90;
c1 = 15; c2 = 15; c3 = 15; c4 = 15; c5 = 15; c6 = 15; c7 = 15;
k1 = 2100; k2 = 1900; k3 = 1050; k4 = 1050; k5 = 950; k6 = 950; k7 = 475;
sig1 = 0.15; sig2 = 0.125; sig3 = 0.125; sig4 = 0.125; sig5 = 0.125; sig6 = 0.1; sig7 = 0.075;

for dp = 1:length(time)
    
    %%%%%%%%%% Noisy Force Simulation %%%%%%%%%%
    
    amp = 10; omega = 10;
    F = [0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t;...
    amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t); 0*t; amp*sin(omega*t)]; % 14x1
    
    snr = 20;
    for i = 1:14
        F(i,:) = F(i,:) + sqrt((1/snr)*(std(F(i,:))^2))*randn(1,length(F(i,:)));
    end
    
    %%%%%%%%%% UKF Algorithm %%%%%%%%%%
    
    z(1,:) = measurement(:,dp); z(2,:) = measurement(:,(length(time)+dp)); z(3,:) = measurement(:,(2*length(time)+dp));
    z(4,:) = measurement(:,(3*length(time)+dp)); z(5,:) = measurement(:,(4*length(time)+dp));
    z(6,:) = measurement(:,(5*length(time)+dp)); z(7,:) = measurement(:,(6*length(time)+dp));
    
    initial_dps{data_dps_pointer} = [k1; k2; k3; k4; k5; k6; k7];
    
    m0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; k1; k2; k3; k4; k5; k6; k7];
    p0 = diag(m0)*diag(m0)'; mk(:,1) = m0; pk = p0.*0.0001+0.000001*eye(size(p0));
    
    R = diag(cov(z')).*eye(7)./25;
    
    for i = 2:length(z)
        fprintf('Iteration %d, step %d\n',dp,i);
        yhikmi = [mk(:,i-1) mk(:,i-1)+sqrt(n+la)*chol(pk,'lower') mk(:,i-1)-sqrt(n+la)*chol(pk,'lower')];
        
        yhik = zeros(size(yhikmi));
        for j=1:length(yhikmi)
            yhik(:,j) = disp_vel_EM(yhikmi(:,j), F(:,i-1));
        end
        
        mkm = 0;
        for j=1:length(yhik)
            mkm = mkm + wm(j)*yhik(:,j);
        end
        
        pkm = 0;
        for j=1:length(yhik)            
            pkm = pkm + wc(j)*(yhik(:,j)-mkm)*(yhik(:,j)-mkm)';
        end
        qc = diag([0, (sig1*dt^(1/2))/m1, 0, (sig2*dt^(1/2))/m2, 0, (sig3*dt^(1/2))/m3...
            , 0, (sig4*mk(7,i-1)*dt^(1/2))/m4, 0, (sig5*dt^(1/2))/m5, 0, (sig6*dt^(1/2))/m6, 0, (sig7*dt^(1/2))/m7,...
            0, 0, 0, 0, 0, 0, 0]);
        Q = qc*qc';
        pkm = pkm+Q;
        
        yhikm = [mkm mkm+sqrt(n+la)*chol(pkm,'lower') mkm-sqrt(n+la)*chol(pkm,'lower')]; zk = zeros(size(z,1),length(yhikm));
        for j=1:length(yhikm)
            y1 = yhikm(1,j); y2 = yhikm(2,j); y3 = yhikm(3,j); y4 = yhikm(4,j); y5 = yhikm(5,j); y6 = yhikm(6,j); y7 = yhikm(7,j);
            y8 = yhikm(8,j); y9 = yhikm(9,j); y10 = yhikm(10,j); y11 = yhikm(11,j); y12 = yhikm(12,j); y13 = yhikm(13,j); y14 = yhikm(14,j);
            k1 = yhikm(15,j); k2 = yhikm(16,j); k3 = yhikm(17,j); k4 = yhikm(18,j);
            k5 = yhikm(19,j); k6 = yhikm(20,j); k7 = yhikm(21,j);
            zk(:,j) = [-(y1*(k1+k2)-c2*y4-k2*y3+y2*(c1+c2))/m1;
                (c2*y2-y3*(k2+k3)+c3*y6+k2*y1+k3*y5-y4*(c2+c3))/m2;
                -(k4*y7-c4*y8-k3*y3-c3*y4+y5*(k3-k4)+alpha*(y5-y7)^3+y6*(c3+c4))/m3;
                (c4*y6+c5*y10-k4*y5+k5*y9+y7*(k4-k5)+alpha*(y5-y7)^3-y8*(c4+c5))/m4;
                (c5*y8-y9*(k5+k6)+c6*y12+k5*y7+k6*y11-y10*(c5+c6))/m5;
                (c6*y10-y11*(k6+k7)+c7*y14+k6*y9+k7*y13-y12*(c6+c7))/m6;
                (c7*y12-c7*y14+k7*y11-k7*y13)/m7];
        end
        
        muk = 0;
        for j=1:length(zk)
            muk = muk + wm(j)*zk(:,j);
        end
        
        sk = 0;
        for j=1:length(zk)
            sk = sk + wc(j)*(zk(:,j)-muk)*(zk(:,j)-muk)';
        end
        sk = sk + R;
        
        ck = 0;
        for j=1:length(zk)
            ck = ck + wc(j)*(yhikm(:,j)-mkm)*(zk(:,j)-muk)';
        end
        
        kk = ck*inv(sk); %#ok<*MINV>
        mk(:,i) = mkm + kk*(z(:,i)-muk);
        pk = pkm - kk*sk*kk';
    end
    
    data_dps{data_dps_pointer} = mk;
    force_dps{data_dps_pointer} = F;
    measurement_dps{data_dps_pointer} = z;
    data_dps_pointer = data_dps_pointer+1;
        
    K1_est(dp) = mk(15,length(z)); K2_est(dp) = mk(16,length(z)); K3_est(dp) = mk(17,length(z)); K4_est(dp) = mk(18,length(z));...
        K5_est(dp) = mk(19,length(z)); K6_est(dp) = mk(20,length(z)); K7_est(dp) = mk(21,length(z));
    k1 = mk(15,length(z)); k2 = mk(16,length(z)); k3 = mk(17,length(z)); k4 = mk(18,length(z)); k5 = mk(19,length(z)); k6 = mk(20,length(z)); k7 = mk(21,length(z));
end

%%%%%%%%%% RMSE %%%%%%%%%%

results(1) = sqrt(sum(((K(1,:)-K1_est').^2))/length(time));
results(2) = sqrt(sum(((K(2,:)-K2_est').^2))/length(time));
results(3) = sqrt(sum(((K(3,:)-K3_est').^2))/length(time));
results(4) = sqrt(sum(((K(5,:)-K5_est').^2))/length(time));
results(5) = sqrt(sum(((K(6,:)-K6_est').^2))/length(time));
results(6) = sqrt(sum(((K(7,:)-K7_est').^2))/length(time));
disp('-------------------------'); 
fprintf('\nStiffness	RMSE\n');
fprintf('   k1    	%3.3f\n',results(1));
fprintf('   k2    	%3.3f\n',results(2));
fprintf('   k3    	%3.3f\n',results(3));
fprintf('   k5    	%3.3f\n',results(4));
fprintf('   k6    	%3.3f\n',results(5));
fprintf('   k7    	%3.3f\n',results(6));

%%%%%%%%%% Plots %%%%%%%%%%

figure('Name','UKF'); set(gcf,'Position',[1000,400,900,500]);

subplot(3,2,1); plot(time,K(1,:),'.r'); hold on; plot(time,K1_est,'.b');
legend('True','Estimated'); xlabel('Time(Days)'); ylabel('K1');
axis([0 t_upper 0.9*K(1,length(time)) 1.1*K(1,1)]);

subplot(3,2,2); plot(time,K(2,:),'.r'); hold on; plot(time,K2_est,'.b');
xlabel('Time(Days)'); ylabel('K2');
axis([0 t_upper 0.9*K(2,length(time)) 1.1*K(2,1)]);

subplot(3,2,3); plot(time,K(3,:),'.r'); hold on; plot(time,K3_est,'.b');
xlabel('Time(Days)'); ylabel('K3');
axis([0 t_upper 0.9*K(3,length(time)) 1.1*K(3,1)]);

subplot(3,2,4); plot(time,K(5,:),'.r'); hold on; plot(time,K5_est,'.b');
xlabel('Time(Days)'); ylabel('K5');
axis([0 t_upper 0.9*K(5,length(time)) 1.1*K(5,1)]);

subplot(3,2,5); plot(time,K(6,:),'.r'); hold on; plot(time,K6_est,'.b');
xlabel('Time(Days)'); ylabel('K6');
axis([0 t_upper 0.9*K(6,length(time)) 1.1*K(6,1)]);

subplot(3,2,6); plot(time,K(7,:),'.r'); hold on; plot(time,K7_est,'.b');
xlabel('Time(Days)'); ylabel('K7');
axis([0 t_upper 0.9*K(7,length(time)) 1.1*K(7,1)]);

param_tr = [K(1,:); K(2,:); K(3,:); K(5,:); K(6,:); K(7,:)];
param_est = [K1_est'; K2_est'; K3_est'; K5_est'; K6_est'; K7_est']; 

for i = 10
    figure('Name','UKF'); set(gcf,'Position',[1000,200,900,500]);
    subplot(7,2,1); plot(t,force_dps{i}(2,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f1');
    legend('Measurement','Actual');
    subplot(7,2,2); plot(t,force_dps{i}(4,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f2');
    legend('Measurement','Actual');
    subplot(7,2,3); plot(t,force_dps{i}(6,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f3');
    legend('Measurement','Actual');
    subplot(7,2,4); plot(t,force_dps{i}(8,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f4');
    legend('Measurement','Actual');
    subplot(7,2,5); plot(t,force_dps{i}(10,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f5');
    legend('Measurement','Actual');
    subplot(7,2,6); plot(t,force_dps{i}(12,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f6');
    legend('Measurement','Actual');
    subplot(7,2,7); plot(t,force_dps{i}(14,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f7');
    legend('Measurement','Actual');
    
    subplot(7,2,8); hold on; plot(t,measurement_dps{i}(1,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(1,:),'r'); xlabel('Time (t)'); ylabel('A1');
    legend('Measurement','Actual');
    subplot(7,2,9); hold on; plot(t,measurement_dps{i}(2,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(2,:),'r'); xlabel('Time (t)'); ylabel('A2');
    legend('Measurement','Actual');
    subplot(7,2,10); hold on; plot(t,measurement_dps{i}(3,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(3,:),'r'); xlabel('Time (t)'); ylabel('A2');
    legend('Measurement','Actual');
    subplot(7,2,11); hold on; plot(t,measurement_dps{i}(4,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(4,:),'r'); xlabel('Time (t)'); ylabel('A2');
    legend('Measurement','Actual');
    subplot(7,2,12); hold on; plot(t,measurement_dps{i}(5,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(5,:),'r'); xlabel('Time (t)'); ylabel('A2');
    legend('Measurement','Actual');
    subplot(7,2,13); hold on; plot(t,measurement_dps{i}(6,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(6,:),'r'); xlabel('Time (t)'); ylabel('A2');
    legend('Measurement','Actual');
    subplot(7,2,14); hold on; plot(t,measurement_dps{i}(7,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(7,:),'r'); xlabel('Time (t)'); ylabel('A2');
    legend('Measurement','Actual');
    

    figure('Name','UKF'); set(gcf,'Position',[1000,200,900,500]);
    subplot(7,2,1); plot(t,y_dps{i}(1,:),'r'); hold on; plot(t,data_dps{i}(1,:),'b');
    xlabel('Time (t)'); ylabel('y1'); legend('True','Filtered Result');
    subplot(7,2,2); plot(t,y_dps{i}(3,:),'r'); hold on; plot(t,data_dps{i}(3,:),'b');
    xlabel('Time (t)'); ylabel('y3'); legend('True','Filtered Result');
    subplot(7,2,3); plot(t,y_dps{i}(5,:),'r'); hold on; plot(t,data_dps{i}(5,:),'b');
    xlabel('Time (t)'); ylabel('y5'); legend('True','Filtered Result');
    subplot(7,2,4); plot(t,y_dps{i}(7,:),'r'); hold on; plot(t,data_dps{i}(7,:),'b');
    xlabel('Time (t)'); ylabel('y7'); legend('True','Filtered Result');
    subplot(7,2,5); plot(t,y_dps{i}(9,:),'r'); hold on; plot(t,data_dps{i}(9,:),'b');
    xlabel('Time (t)'); ylabel('y9'); legend('True','Filtered Result');
    subplot(7,2,6); plot(t,y_dps{i}(11,:),'r'); hold on; plot(t,data_dps{i}(11,:),'b');
    xlabel('Time (t)'); ylabel('y11'); legend('True','Filtered Result');
    subplot(7,2,7); plot(t,y_dps{i}(13,:),'r'); hold on; plot(t,data_dps{i}(13,:),'b');
    xlabel('Time (t)'); ylabel('y13'); legend('True','Filtered Result');
    
    subplot(7,2,8); plot(t,y_dps{i}(2,:),'r'); hold on; plot(t,data_dps{i}(2,:),'b');
    xlabel('Time (t)'); ylabel('y2'); legend('True','Filtered Result');
    subplot(7,2,9); plot(t,y_dps{i}(4,:),'r'); hold on; plot(t,data_dps{i}(4,:),'b');
    xlabel('Time (t)'); ylabel('y4'); legend('True','Filtered Result');
    subplot(7,2,10); plot(t,y_dps{i}(6,:),'r'); hold on; plot(t,data_dps{i}(6,:),'b');
    xlabel('Time (t)'); ylabel('y6'); legend('True','Filtered Result');
    subplot(7,2,11); plot(t,y_dps{i}(8,:),'r'); hold on; plot(t,data_dps{i}(8,:),'b');
    xlabel('Time (t)'); ylabel('y8'); legend('True','Filtered Result');
    subplot(7,2,12); plot(t,y_dps{i}(10,:),'r'); hold on; plot(t,data_dps{i}(10,:),'b');
    xlabel('Time (t)'); ylabel('y10'); legend('True','Filtered Result');
    subplot(7,2,13); plot(t,y_dps{i}(12,:),'r'); hold on; plot(t,data_dps{i}(12,:),'b');
    xlabel('Time (t)'); ylabel('y12'); legend('True','Filtered Result');
    subplot(7,2,14); plot(t,y_dps{i}(14,:),'r'); hold on; plot(t,data_dps{i}(14,:),'b');
    xlabel('Time (t)'); ylabel('y14'); legend('True','Filtered Result');

    
    figure('Name','UKF'); set(gcf,'Position',[1000,200,900,500]);
    subplot(7,1,1); plot(t,kgt_dps{i}(1)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(1)*ones(1,length(t)),'c'); plot(t,data_dps{i}(15,:),'b');
    legend('True','Ground Truth','Filtered Result');
    xlabel('Time (t)');ylabel('k1');
    
    subplot(7,1,2); plot(t,kgt_dps{i}(2)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(2)*ones(1,length(t)),'c'); plot(t,data_dps{i}(16,:),'b');
    legend('True','Ground Truth','Filtered Result');
    xlabel('Time (t)');ylabel('k2');
    
    subplot(7,1,3); plot(t,kgt_dps{i}(3)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(3)*ones(1,length(t)),'c'); plot(t,data_dps{i}(17,:),'b');
    legend('True','Ground Truth','Filtered Result');
    xlabel('Time (t)');ylabel('k2');
    
    subplot(7,1,4); plot(t,kgt_dps{i}(5)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(5)*ones(1,length(t)),'c'); plot(t,data_dps{i}(19,:),'b');
    legend('True','Ground Truth','Filtered Result');
    xlabel('Time (t)');ylabel('k2');
    
    subplot(7,1,5); plot(t,kgt_dps{i}(6)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(6)*ones(1,length(t)),'c'); plot(t,data_dps{i}(20,:),'b');
    legend('True','Ground Truth','Filtered Result');
    xlabel('Time (t)');ylabel('k2');
    
    subplot(7,1,6); plot(t,kgt_dps{i}(7)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(7)*ones(1,length(t)),'c'); plot(t,data_dps{i}(21,:),'b');
    legend('True','Ground Truth','Filtered Result');
    xlabel('Time (t)');ylabel('k2');
end

save('workspace_saved_dps_7dof')

%%%%%%%%%% Functions %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = disp_vel_IT(y, F, K, C)
global dt

m1 = 20; m2 = 20; m3 = 10; m4 = 10; m5 = 10; m6 = 10; m7 = 5; alpha = 100;
sig1 = 0.1; sig2 = 0.1; sig3 = 0.1; sig4 = 0.1; sig5 = 0.1; sig6 = 0.1; sig7 = 0.1;

k1 = K(1); k2 = K(2); k3 = K(3); k4 = K(4); k5 = K(5); k6 = K(6); k7 = K(7);
c1 = C(1); c2 = C(2); c3 = C(3); c4 = C(4); c5 = C(5); c6 = C(6); c7 = C(7);

f1 = F(2); f2 = F(4); f3 = F(6); f4 = F(8); f5 = F(10); f6 = F(12); f7 = F(14);

y1 = y(1); y2 = y(2); y3 = y(3); y4 = y(4); y5 = y(5); y6 = y(6); y7 = y(7);
y8 = y(8); y9 = y(9); y10 = y(10); y11 = y(11); y12 = y(12); y13 = y(13); y14 = y(14);

deltamat = [sqrt(dt)            0;
    dt^1.5/2    dt^1.5/(2*sqrt(3))];

dZ = (deltamat(2,:)*randn(2,7))'; % 7 x 1
dW = (deltamat(1,:)*randn(2,7))';  % 7 x 1

a1 = y2;
a2 = (1/m1*(-(c1+c2)*y2-(k1+k2)*y1+c2*y4+k2*y3))+(f1/m1);
a3 = y4;
a4 = (1/m2*(-(c2+c3)*y4-(k2+k3)*y3+c2*y2+k2*y1+c3*y6+k3*y5))+(f2/m2);
a5 = y6;
a6 = (1/m3*(-(c3+c4)*y6-(k3-k4)*y5-alpha*(y5-y7)^3+c3*y4+k3*y3+c4*y8-k4*y7))+(f3/m3);
a7 = y8;
a8 = (1/m4*(-(c4+c5)*y8-(k5-k4)*y7-alpha*(y7-y5)^3+c4*y6-k4*y5+c5*y10+k5*y9))+(f4/m4);
a9 = y10;
a10 = (1/m5*(-(c5+c6)*y10-(k5+k6)*y9+c5*y8+k5*y7+c6*y12+k6*y11))+(f5/m5);
a11 = y12;
a12 = (1/m6*(-(c6+c7)*y12-(k6+k7)*y11+c6*y10+k6*y9+c7*y14+k7*y13))+(f6/m6);
a13 = y14;
a14 = (1/m7*(-c7*y14-k7*y13+c7*y12+k7*y11))+(f7/m7);
a = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14]; % 14x1
    
b = [0, 0, 0, 0, 0, 0, 0;
    sig1/m1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0;
    0, sig2/m2, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0;
    0, 0, sig3/m3, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, sig4/m4*y7, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, sig5/m5, 0, 0;
    0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, sig6/m6, 0;
    0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, sig7/m7];

L0a = [ (f1)/m1 - (y1*(k1 + k2) - c2*y4 - k2*y3 + y2*(c1 + c2))/m1;
    (k2*y4)/m1 + ((c1 + c2)*((y1*(k1 + k2) - c2*y4 - k2*y3 + y2*(c1 + c2))/m1 - (f1)/m1))/m1 + (c2*((c2*y2 - y3*(k2 + k3) + c3*y6 + k2*y1 + k3*y5 - y4*(c2 + c3))/m2 + (f2)/m2))/m1 - (y2*(k1 + k2))/m1;
    (c2*y2 - y3*(k2 + k3) + c3*y6 + k2*y1 + k3*y5 - y4*(c2 + c3))/m2 + (f2)/m2;
    (k2*y2)/m2 - (c3*((k4*y7 - c4*y8 - k3*y3 - c3*y4 + y5*(k3 - k4) + alpha*(y5 - y7)^3 + y6*(c3 + c4))/m3 - (f3)/m3))/m2 + (k3*y6)/m2 - (c2*((y1*(k1 + k2) - c2*y4 - k2*y3 + y2*(c1 + c2))/m1 - (f1)/m1))/m2 - ((c2 + c3)*((c2*y2 - y3*(k2 + k3) + c3*y6 + k2*y1 + k3*y5 - y4*(c2 + c3))/m2 + (f2)/m2))/m2 - (y4*(k2 + k3))/m2;
    (f3)/m3 - (k4*y7 - c4*y8 - k3*y3 - c3*y4 + y5*(k3 - k4) + alpha*(y5 - y7)^3 + y6*(c3 + c4))/m3;
    (c4*((c4*y6 + c5*y10 - k4*y5 + k5*y9 + y7*(k4 - k5) + alpha*(y5 - y7)^3 - y8*(c4 + c5))/m4 + (f4)/m4))/m3 - (y6*(k3 - k4 + 3*alpha*(y5 - y7)^2))/m3 + (k3*y4)/m3 - (y8*(k4 - 3*alpha*(y5 - y7)^2))/m3 + ((c3 + c4)*((k4*y7 - c4*y8 - k3*y3 - c3*y4 + y5*(k3 - k4) + alpha*(y5 - y7)^3 + y6*(c3 + c4))/m3 - (f3)/m3))/m3 + (c3*((c2*y2 - y3*(k2 + k3) + c3*y6 + k2*y1 + k3*y5 - y4*(c2 + c3))/m2 + (f2)/m2))/m3;
    (c4*y6 + c5*y10 - k4*y5 + k5*y9 + y7*(k4 - k5) + alpha*(y5 - y7)^3 - y8*(c4 + c5))/m4 + (f4)/m4;
    (k5*y10)/m4 - (y8*(k5 - k4 + 3*alpha*(y5 - y7)^2))/m4 - (c4*((k4*y7 - c4*y8 - k3*y3 - c3*y4 + y5*(k3 - k4) + alpha*(y5 - y7)^3 + y6*(c3 + c4))/m3 - (f3)/m3))/m4 - (y6*(k4 - 3*alpha*(y5 - y7)^2))/m4 - ((c4 + c5)*((c4*y6 + c5*y10 - k4*y5 + k5*y9 + y7*(k4 - k5) + alpha*(y5 - y7)^3 - y8*(c4 + c5))/m4 + (f4)/m4))/m4 + (c5*((c5*y8 - y9*(k5 + k6) + c6*y12 + k5*y7 + k6*y11 - y10*(c5 + c6))/m5 + (f5)/m5))/m4;
    (c5*y8 - y9*(k5 + k6) + c6*y12 + k5*y7 + k6*y11 - y10*(c5 + c6))/m5 + (f5)/m5;
    (c5*((c4*y6 + c5*y10 - k4*y5 + k5*y9 + y7*(k4 - k5) + alpha*(y5 - y7)^3 - y8*(c4 + c5))/m4 + (f4)/m4))/m5 + (k5*y8)/m5 + (k6*y12)/m5 - ((c5 + c6)*((c5*y8 - y9*(k5 + k6) + c6*y12 + k5*y7 + k6*y11 - y10*(c5 + c6))/m5 + (f5)/m5))/m5 + (c6*((c6*y10 - y11*(k6 + k7) + c7*y14 + k6*y9 + k7*y13 - y12*(c6 + c7))/m6 + (f6)/m6))/m5 - (y10*(k5 + k6))/m5;
    (c6*y10 - y11*(k6 + k7) + c7*y14 + k6*y9 + k7*y13 - y12*(c6 + c7))/m6 + (f6)/m6;
    (k6*y10)/m6 + (k7*y14)/m6 - ((c6 + c7)*((c6*y10 - y11*(k6 + k7) + c7*y14 + k6*y9 + k7*y13 - y12*(c6 + c7))/m6 + (f6)/m6))/m6 + (c7*((c7*y12 - c7*y14 + k7*y11 - k7*y13)/m7 + (f7)/m7))/m6 + (c6*((c5*y8 - y9*(k5 + k6) + c6*y12 + k5*y7 + k6*y11 - y10*(c5 + c6))/m5 + (f5)/m5))/m6 - (y12*(k6 + k7))/m6;
    (c7*y12 - c7*y14 + k7*y11 - k7*y13)/m7 + (f7)/m7;
    (k7*y12)/m7 - (k7*y14)/m7 - (c7*((c7*y12 - c7*y14 + k7*y11 - k7*y13)/m7 + (f7)/m7))/m7 + (c7*((c6*y10 - y11*(k6 + k7) + c7*y14 + k6*y9 + k7*y13 - y12*(c6 + c7))/m6 + (f6)/m6))/m7];

L0b = [0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,(sig4*y8)/m4,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0];

Lja = [sig1/m1, 0, 0, 0, 0, 0, 0;
    -(sig1*(c1+c2))/m1^2, (c2*sig2)/(m1*m2), 0, 0, 0, 0, 0;
    0, sig2/m2, 0, 0, 0, 0, 0;
    (c2*sig1)/(m1*m2), -(sig2*(c2+c3))/m2^2, (c3*sig3)/(m2*m3), 0, 0, 0, 0;
    0, 0, sig3/m3, 0, 0, 0, 0;
    0, (c3*sig2)/(m2*m3), -(sig3*(c3+c4))/m3^2, (c4*sig4*y7)/(m3*m4), 0, 0, 0;
    0, 0, 0, (sig4*y7)/m4, 0, 0, 0;
    0, 0, (c4*sig3)/(m3*m4), -(sig4*y7*(c4+c5))/m4^2, (c5*sig5)/(m4*m5), 0, 0;
    0, 0, 0, 0, sig5/m5, 0, 0;
    0, 0, 0, (c5*sig4*y7)/(m4*m5), -(sig5*(c5+c6))/m5^2, (c6*sig6)/(m5*m6), 0;
    0, 0, 0, 0, 0, sig6/m6, 0;
    0, 0, 0, 0, (c6*sig5)/(m5*m6), -(sig6*(c6+c7))/m6^2, (c7*sig7)/(m6*m7);
    0, 0, 0, 0, 0, 0, sig7/m7;
    0, 0, 0, 0, 0, (c7*sig6)/(m6*m7), -(c7*sig7)/m7^2];

y = y + a*dt + b*dW + L0b*(dW*dt-dZ) + Lja*dZ + 0.5*L0a*dt^2;
end

function y_1 = disp_vel_EM(mk, F)
global dt

f1 = F(2); f2 = F(4); f3 = F(6); f4 = F(8); f5 = F(10); f6 = F(12); f7 = F(14);
m1 = 20; m2 = 20; m3 = 10; m4 = 10; m5 = 10; m6 = 10; m7 = 5;
c1 = 15; c2 = 15; c3 = 15; c4 = 15; c5 = 15; c6 = 15; c7 = 15;
alpha = 90;

y1 = mk(1); y2 = mk(2); y3 = mk(3); y4 = mk(4); y5 = mk(5); y6 = mk(6); y7 = mk(7); y8 = mk(8);
y9 = mk(9); y10 = mk(10); y11 = mk(11); y12 = mk(12); y13 = mk(13); y14 = mk(14);
k1 = mk(15); k2 = mk(16); k3 = mk(17); k4 = mk(18); k5 = mk(19); k6 = mk(20); k7 = mk(21);

y = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12; y13; y14; k1; k2; k3; k4; k5; k6; k7];

a1 = y2;
a2 = (1/m1*(-(c1+c2)*y2-(k1+k2)*y1+c2*y4+k2*y3))+(f1/m1);
a3 = y4;
a4 = (1/m2*(-(c2+c3)*y4-(k2+k3)*y3+c2*y2+k2*y1+c3*y6+k3*y5))+(f2/m2);
a5 = y6;
a6 = (1/m3*(-(c3+c4)*y6-(k3-k4)*y5-alpha*(y5-y7)^3+c3*y4+k3*y3+c4*y8-k4*y7))+(f3/m3);
a7 = y8;
a8 = (1/m4*(-(c4+c5)*y8-(k5-k4)*y7-alpha*(y7-y5)^3+c4*y6-k4*y5+c5*y10+k5*y9))+(f4/m4);
a9 = y10;
a10 = (1/m5*(-(c5+c6)*y10-(k5+k6)*y9+c5*y8+k5*y7+c6*y12+k6*y11))+(f5/m5);
a11 = y12;
a12 = (1/m6*(-(c6+c7)*y12-(k6+k7)*y11+c6*y10+k6*y9+c7*y14+k7*y13))+(f6/m6);
a13 = y14;
a14 = (1/m7*(-c7*y14-k7*y13+c7*y12+k7*y11))+(f7/m7);
a = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12; a13; a14; 0; 0; 0; 0; 0; 0; 0];

y_1 = y + a*dt;
end

%#ok<*MINV>