%%%%%%%%%% Data Points Generated For 2DOF System %%%%%%%%%%
%%%%%%%%%% Subjected To Deterministic And Stochastic Forces  %%%%%%%%%%
%%%%%%%%%% Shailesh Garg  %%%%%%%%%%

clc
clear
close all
rng(30)

%%%%%%%%%% Time Vector For Individual Data Point %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dt
dt = 0.001; t = 0:dt:5;

%%%%%%%%%% Stiffness And Damping Simulation %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = 0:50:10000; degradation = -1*0.00005;
K1 = 1000*exp(degradation*time); K2 = 500*exp(degradation*time);
C1 = 10*exp(degradation*time); C2 = 5*exp(degradation*time);

%%%%%%%%%% Iterations %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K1_est = zeros(length(time),1); K2_est = zeros(length(time),1);
C1_est = zeros(length(time),1); C2_est = zeros(length(time),1);

nsim = 1; %%%%%%%%%% No. Of Simulations For A Particular Data Point %%%%%%%%%%

lambda1 = 10; lambda2 = 10; omega1 = 10; omega2 = 10;
f1 = lambda1*sin(omega1*t); f2 = lambda2*sin(omega2*t);

%%%%%%%%%% Data Simulation Using Taylor 1.5 Strong %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dps_pointer = 1; data =  zeros(length(t),2*length(time));

for itr = 1:length(time)
    k1 = K1(itr); k2 = K2(itr); c1 = C1(itr); c2 = C2(itr);
    m1 = 20; m2 = 10; alpha = 100; z = zeros(2,length(t)); acc = zeros(2,length(t));
    
    for sim = 1:nsim
        y = zeros(4,length(t)); acc_sim = zeros(2,length(t));
        fprintf('Iteration %d Simulation %d\n',itr,sim);
        for i = 1:length(t)-1
            y_it = disp_vel_IT(y(1,i), y(2,i), y(3,i), y(4,i), f1(i), f2(i), k1, k2, c1, c2);
            y(:,i+1) = y_it;
        end
        m = [m1 0; 0 m2]; c = [c1+c2 -c2; -c2 c2];
        for i = 1:length(y)
            k = [k1+alpha*y(1,i)^2+k2 -k2; -k2 k2];
            acc_sim(:,i) = -inv(m)*(k*[y(1,i); y(2,i)]+c*[y(3,i); y(4,i)]);
        end
        acc = acc + acc_sim;
    end

    acc_dps{data_dps_pointer} = acc/nsim;
    y_dps{data_dps_pointer} = y;
    kgt_dps{data_dps_pointer} = [k1; k2];
    data_dps_pointer = data_dps_pointer+1;

    data(:,itr) = acc(1,:)/nsim; data(:,(length(time)+itr)) = acc(2,:)/nsim;
end

%%%%%%%%%% Measurement Simulation %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr = 50; measurement =  zeros(length(t),2*length(time));

for itr = 1:length(time)
    measurement(:,itr) = data(:,itr)+sqrt((1/snr)*(std(data(:,itr))^2))*randn(length(data(:,itr)),1);
    measurement(:,(length(time)+itr)) = data(:,(length(time)+itr))...
        +sqrt((1/snr)*(std(data(:,(length(time)+itr)))^2))*randn(length(data(:,(length(time)+itr))),1);
end

%%%%%%%%%% UKF Application Using Euler Maruyama %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m1 = 20; m2 = 10; alpha = 90; s1 = 0.05; s2 = 0.05;
k1_ini = 1050; k2_ini = 475; c1 = 7.5; c2 = 7.5;
m0 = [0; 0; 0; 0; k1_ini; k2_ini];

%%%%%%%%%% UKF Weights Calculations %%%%%%%%%%

n = size(m0,1); al = 0.001; beta = 2; kappa = 0;
la = al^2*(n+kappa)-n;
wm = zeros(2*n+1,1); wc = zeros(2*n+1,1);
for j=1:2*n+1
    if j==1
        wm(j) = la/(n+la);
        wc(j) = la/(n+la)+(1-al^2+beta);
    else
        wm(j) =1/(2*(n+la));
        wc(j) = wm(j);
    end
end

%%%%%%%%%% UKF Algorithm %%%%%%%%%%

data_dps_pointer = 1;

for itr = 1:length(time)
    lambda1 = 10; lambda2 = 10; omega1 = 10; omega2 = 10;
    f1 = lambda1*sin(omega1*t); f2 = lambda2*sin(omega2*t);
    snr = 20;
    f1 = f1 + sqrt((1/snr)*(std(f1)^2))*randn(1,length(f1));
    f2 = f2 + sqrt((1/snr)*(std(f2)^2))*randn(1,length(f2));
    
    z(1,:) = measurement(:,itr); z(2,:) = measurement(:,(length(time)+itr));
    m0 = [0; 0; 0; 0; k1_ini; k2_ini]; p0 = diag(m0)*diag(m0)'.*0.0001+0.000001*eye(length(m0));
    R = diag(cov(z')).*eye(2);
    
    mk(:,1) = m0; pk = p0;
    for i = 2:length(z)
        fprintf('Iteration %d, step %d\n',itr,i);
        yhikmi = [mk(:,i-1) mk(:,i-1)+sqrt(n+la)*chol(pk,'lower') mk(:,i-1)-sqrt(n+la)*chol(pk,'lower')];
        yhik = zeros(size(yhikmi));
        for j=1:length(yhikmi)
            yhik(:,j) = disp_vel_EM(yhikmi(1,j), yhikmi(2,j), yhikmi(3,j), yhikmi(4,j), yhikmi(5,j),...
                yhikmi(6,j), f1(i-1), f2(i-1));
        end
        
        mkm = 0;
        for j=1:length(yhik)
            mkm = mkm + wm(j)*yhik(:,j);
        end
        pkm = 0;
        for j=1:length(yhik)
            pkm = pkm + wc(j)*((yhik(:,j)-mkm)*(yhik(:,j)-mkm)');
        end
        qc = diag([0, 0, (s1*dt^(1/2))/m1, (s2*dt^(1/2))/m2 0, 0]);
        Q = qc*qc';
        pkm = pkm+Q;
        
        yhikm = [mkm mkm+sqrt(n+la)*chol(pkm,'lower') mkm-sqrt(n+la)*chol(pkm,'lower')];
        zk = zeros(size(z,1),length(yhikm));
        for j=1:length(yhikm)
            y1 = yhikm(1,j); y2 = yhikm(2,j); y3 = yhikm(3,j); y4 = yhikm(4,j);
            k1 = yhikm(5,j); k2 = yhikm(6,j);
            zk(:,j) = [(c2*y4+k2*y2-y1*(alpha*y1^2+k1+k2)-y3*(c1+c2))/m1;
                (c2*y3-c2*y4+k2*y1-k2*y2)/m2];
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
    force_dps{data_dps_pointer} = [f1; f2];
    measurement_dps{data_dps_pointer} = z;
    initial_dps{data_dps_pointer} = [k1_ini; k2_ini];
    data_dps_pointer = data_dps_pointer+1;
    
    K1_est(itr) = mk(5,length(z)); K2_est(itr) = mk(6,length(z));
    k1_ini = mk(5,length(z)); k2_ini = mk(6,length(z));
end

results(1) = sqrt(sum(((K1-K1_est').^2))/length(time));
results(2) = sqrt(sum(((K2-K2_est').^2))/length(time));
disp('-------------------------');
disp('RMSE [K1; K2] = '); disp(results);

%%%%%%%%%% Plots %%%%%%%%%%

figure('Name','UKF'); set(gcf,'Position',[1000,400,900,350]);
plot(time,K1,'.r'); hold on; plot(time,K1_est,'.b');
legend('Ground Truth','Estimated'); xlabel('Time(Days)'); ylabel('K1');
plot(time(1),K1_est(1),'ob'); plot(time(120),K1_est(120),'ob');
figure('Name','UKF'); set(gcf,'Position',[1000,400,900,350]);
plot(time,K2,'.r'); hold on; plot(time,K2_est,'.b');
legend('Ground Truth','Estimated'); xlabel('Time(Days)'); ylabel('K2');
plot(time(1),K2_est(1),'ob'); plot(time(120),K2_est(120),'ob');

for i = 1
    figure('Name','UKF'); set(gcf,'Position',[1000,200,900,500]);
    subplot(2,2,1); plot(t,force_dps{i}(1,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f1');
    legend('Measurement','Ground Truth');
    subplot(2,2,2); plot(t,force_dps{i}(2,:),'Color',[0.5 0.5 0.5]); hold on;
    plot(t,10*sin(10*t),'r'); xlabel('Time (t)'); ylabel('f2');
    legend('Measurement','Ground Truth');
    subplot(2,2,3); hold on; plot(t,measurement_dps{i}(1,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(1,:),'r'); xlabel('Time (t)'); ylabel('A1');
    legend('Measurement','Ground Truth');
    subplot(2,2,4); hold on; plot(t,measurement_dps{i}(2,:),'Color',[0.5 0.5 0.5]);
    plot(t,acc_dps{i}(2,:),'r'); xlabel('Time (t)'); ylabel('A2');
    legend('Measurement','Ground Truth');

    figure('Name','UKF'); set(gcf,'Position',[1000,200,900,500]);
    subplot(2,2,1); plot(t,y_dps{i}(1,:),'r'); hold on; plot(t,data_dps{i}(1,:),'b');
    xlabel('Time (t)'); ylabel('y1'); legend('Ground Truth','Filtered Result');
    subplot(2,2,2); plot(t,y_dps{i}(2,:),'r'); hold on; plot(t,data_dps{i}(2,:),'b');
    xlabel('Time (t)'); ylabel('y2'); legend('Ground Truth','Filtered Result');
    subplot(2,2,3); plot(t,y_dps{i}(3,:),'r'); hold on; plot(t,data_dps{i}(3,:),'b');
    xlabel('Time (t)'); ylabel('y3'); legend('Ground Truth','Filtered Result');
    subplot(2,2,4); plot(t,y_dps{i}(4,:),'r'); hold on; plot(t,data_dps{i}(4,:),'b');
    xlabel('Time (t)'); ylabel('y4'); legend('Ground Truth','Filtered Result');

    figure('Name','UKF'); set(gcf,'Position',[1000,200,900,500]);
    subplot(2,1,1); plot(t,kgt_dps{i}(1)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(1)*ones(1,length(t)),'c'); plot(t,data_dps{i}(5,:),'b');
    legend('Ground Truth','Initial Guess','Filtered Result');
    xlabel('Time (t)');ylabel('k1');
    ylim([ 0.99*min([initial_dps{i}(1),kgt_dps{i}(1)]) 1.01*max([initial_dps{i}(1),kgt_dps{i}(1)]) ]);
    subplot(2,1,2); plot(t,kgt_dps{i}(2)*ones(1,length(t)),'r'); hold on;
    plot(t,initial_dps{i}(2)*ones(1,length(t)),'c'); plot(t,data_dps{i}(6,:),'b');
    legend('Ground Truth','Initial Guess','Filtered Result');
    xlabel('Time (t)');ylabel('k2');
    ylim([ 0.99*min([initial_dps{i}(2),kgt_dps{i}(2)]) 1.01*max([initial_dps{i}(2),kgt_dps{i}(2)]) ]);
end

save('workspace_saved_dps_2dof')

%%%%%%%%%% Functions %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = disp_vel_IT(y1, y2, y3, y4, f1, f2, k1, k2, c1, c2)
global dt

m1 = 20; m2 = 10; alpha = 100; s1 = 0.1; s2 = 0.1;
r1 = randn; r2 = randn; r3 = randn; r4 = randn;

y = [y1; y2; y3; y4];
a = [ y3;
    y4;
    ((-(c1*y3+c2*y3-c2*y4+k1*y1+k2*y1-k2*y2+alpha*y1^3)/m1)+(f1/m1));
    (((c2*y3-c2*y4+k2*y1-k2*y2)/m2)+(f2/m2))];

L0a = [-(c1*y3 + c2*y3 - c2*y4 + k1*y1 + k2*y1 - k2*y2 - f1 + alpha*y1^3)/m1;
    (c2*y3 - c2*y4 + k2*y1 - k2*y2 + f2)/m2;
    (c1^2*m2*y3 + c2^2*m1*y3 - c2^2*m1*y4 + c2^2*m2*y3 - c2^2*m2*y4 + 2*c1*c2*m2*y3 - c1*c2*m2*y4 + c1*k1*m2*y1 + c1*k2*m2*y1 + c2*k1*m2*y1 + c2*k2*m1*y1 - c1*k2*m2*y2 - c2*k2*m1*y2 + c2*k2*m2*y1 - c2*k2*m2*y2 - k1*m1*m2*y3 - k2*m1*m2*y3 + k2*m1*m2*y4 - c1*m2*f1 - c2*m2*f1 + c2*m1*f2 + alpha*c1*m2*y1^3 + alpha*c2*m2*y1^3 - 3*alpha*m1*m2*y1^2*y3)/(m1^2*m2);
    -(c2^2*m1*y3 - c2^2*m1*y4 + c2^2*m2*y3 - c2^2*m2*y4 + c1*c2*m2*y3 + c2*k1*m2*y1 + c2*k2*m1*y1 - c2*k2*m1*y2 + c2*k2*m2*y1 - c2*k2*m2*y2 - k2*m1*m2*y3 + k2*m1*m2*y4 - c2*m2*f1 + c2*m1*f2 + alpha*c2*m2*y1^3)/(m1*m2^2)];

b = [0;0;s1/m1;s2/m2];
L1a = [s1/m1;0;-(s1*(c1+c2))/m1^2;(c2*s1)/(m1*m2)];
L2a = [0;s2/m2;(c2*s2)/(m1*m2);-(c2*s2)/m2^2];
dz1 = (dt^(3/2)*r1)/2 + (3^(1/2)*dt^(3/2)*r2)/6;
dz2 = (dt^(3/2)*r3)/2 + (3^(1/2)*dt^(3/2)*r4)/6;
dw = [0; 0; dt^(1/2)*r1; dt^(1/2)*r3];

for i = 1:4
    y(i) = y(i)+a(i)*dt+b(i)*dw(i)+L1a(i)*dz1+L2a(i)*dz2+0.5*L0a(i)*dt^2;
end
end

function y_1 = disp_vel_EM(y1, y2, y3, y4, k1, k2, f1, f2)
global dt

m1 = 20; m2 = 10; alpha = 90; c1 = 7.5; c2 = 7.5;

y = [y1; y2; y3; y4; k1; k2]; y_1 = zeros(6, 1);
a = [y3;
    y4;
    (-(c1*y3+c2*y3-c2*y4+k1*y1+k2*y1-k2*y2+alpha*y1^3)/m1)+(f1/m1);
    ((c2*y3-c2*y4+k2*y1-k2*y2)/m2)+(f2/m2);
    0;
    0];

for i = 1:length(a)
    y_1(i) = y(i)+a(i)*dt;
end
end

%#ok<*MINV>