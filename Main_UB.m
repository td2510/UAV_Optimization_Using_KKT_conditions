%% 
clc; clear all; close all;
tic
%addpath('D:\Matlab OPtimization\cvx');
%addpath('C:\Users\AS\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\YALMIP-master\@sdpvar');
%addpath('D:\Matlab OPtimization \YALMIP-master');
cvx_setup;
cvx_solver SDPT3;

tic

global P_s V_max sigma_sq H delta_t omega_0 P_c alpha miu q_I1 q_F1 q_I2 q_F2 w_s ....
    w_d epsilon sigma Euler eta_max S Theta Theta_0 B P_u P_h

rng('default')

% rng(14,'twister')

%% Initialization
V_max = 20; % m/s
% sigma_sq = 10^(-6); %-90 dB = -60 dBm AWGN
sigma_sq = 10^(-6);

omega_0 = 10^(-3); % -30dB

P_c = 10^(-3); % 10^-6 Watt = 10^-3 mW power consumption of Backscatter device
P_u = 5; % 10 mW 0 dBm
alpha = 2.2; % Channel coefficient
miu = 0.84; % energy harvesting efficiency
epsilon = 1e-4; % tolerance value
sigma = 0.5; % caching gain coefficience: part of file that UB cache
delta_t = 0.5; % second
Euler = 0.5772156649; % Euler constant 8.367.1 in Table of Integral
q_I1 = [0;10;10];q_F1 = [20;10;10];%???
q_I2 = [0;10;5];q_F2 = [20;10;5];%???
w_s = [5;0;0]; w_d = [15;0;0];%???

% P_s = 1000; % transmission power of source to charge the UAV 1W = 10^3 mW
T_array = [20]; % second
P_s = 10^1.6; % 16 dBm
P_h = 1e7;% 10^4 W = 10^7 mW this is P_WPT in the paper

H = 10; % m
eta_max = 0.5; % maximum value of reflection coefficient
B = 20; % Mbits
S = 70; % demanded data of destination in Mbits
Theta = exp(-Euler)*omega_0/sigma_sq;
Theta_0 = exp(-Euler)*omega_0/sigma_sq;

for i=1:length(T_array)    
    N = T_array(:,i)/delta_t;% Number of time slot
    %% Initialize trajectory q for UAV 1
% q_e1 = [0.5*(w_s(1,:)+w_d(1,:));0;5]; %??? 
 q_e1 = [0.5*(w_s(1,:)+w_d(1,:));0;10]; %??? 
 q_start = q_I1;
 N_segment = (N-1)/2;
 m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:)); % m = -1
 center = q_start(2,:) - q_start(1,:) * m; % = 10
 q_j_1= [];
 for x=q_start(1,:):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:) %0:10/19.5:10
    y = m *x + center; % y=10
%     z = m *x/2 + center; %???
%     q_j = [q_j,[x;y;z]]; %???
    q_j_1 = [q_j_1,[x;y;10]];%???
 end
q_j_1 = [q_j_1,q_e1];
%% 
%  q_start = [0.5*(w_s(1,:)+w_d(1,:));0;5]; %???
 q_start = [0.5*(w_s(1,:)+w_d(1,:));0;10]; %???
 q_e1 = q_F1; % = [20;10;10];
 N_segment = (N-1)/2;
 m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:)); % m = 1
 center = q_start(2,:) - q_start(1,:) * m; % = -10
 for x=q_start(1,:)+(q_e1(1,:)-q_start(1,:))/(N_segment):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:)% 10.5128:10/19.5:20
    y = m *x + center;
%     z = m *x/2 ; %???
%     q_j = [q_j,[x;y;z]]; %???
    q_j_1 = [q_j_1,[x;y;10]];%???
 end 
q_j_1 = [q_j_1,q_e1];
    %% Initialize trajectory q for UAV 2
% q_e1 = [0.5*(w_s(1,:)+w_d(1,:));0;5]; %??? 
 q_e1 = [0.5*(w_s(1,:)+w_d(1,:));0;5]; %??? 
 q_start = q_I2;
 N_segment = (N-1)/2;
 m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:)); % m = -1
 center = q_start(2,:) - q_start(1,:) * m; % = 10
 q_j_2= [];
 for x=q_start(1,:):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:) %0:10/19.5:10
    y = m *x + center; % y=10
%     z = m *x/2 + center; %???
%     q_j = [q_j,[x;y;z]]; %???
    q_j_2 = [q_j_2,[x;y;5]];%???
 end
q_j_2 = [q_j_2,q_e1];
%% 
 q_start = [0.5*(w_s(1,:)+w_d(1,:));0;5]; %???
 q_e1 = q_F2; % = [20;10;10];
 N_segment = (N-1)/2;
 m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:)); % m = 1
 center = q_start(2,:) - q_start(1,:) * m; % = -10
 for x=q_start(1,:)+(q_e1(1,:)-q_start(1,:))/(N_segment):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:)% 10.5128:10/19.5:20
    y = m *x + center;
%     z = m *x/2 ; %???
%     q_j = [q_j,[x;y;z]]; %???
    q_j_2= [q_j_2,[x;y;5]];%???
 end 
q_j_2 = [q_j_2,q_e1];
    %% Test 
%       plot3(q_j_1(1,:),q_j_1(2,:),q_j_1(3,:)) %???
%       hold on
%       plot3(q_j_2(1,:),q_j_2(2,:),q_j_2(3,:)) %???
    %% Initial backscatter coefficient and DTS
    tau_j = 0.5*ones(1,N);
    eta_j = eta_max*ones(1,N);
    %% Linear EH  model: Alternating optimization algorithm
    [q_L_1,tau_L_1,eta_L_1,g_x_L_1] = Linear_model_1(tau_j,eta_j,q_j_1,N);
    [q_L_2,tau_L_2,eta_L_2,g_x_L_2] = Linear_model_2(tau_j,eta_j,q_j_2,N);
    if P_h==5e5
        q_L_5e5 = q_L;
    elseif P_h==1e5
        q_L1e5 = q_L;
    elseif P_h==3e5
        q_L3e5 = q_L;
    elseif P_h==1e6
        q_L1e6 = q_L;
    elseif P_h==1e7
        q_L_11e7 = q_L_1;
        q_L_21e7 = q_L_2;
    end
end



% %% Test plot Trajectory
% figure(1)
% plot(q_L(1,:),q_L(2,:),'b-','LineWidth',1.5)
% xlabel('x (m)','FontSize',14) 
% ylabel('y (m)','FontSize',14) 
% legend('Linear EH')
% plot(0,10,'ks','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% plot(20,10,'ks','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% plot(5,0,'ko','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% plot(15,0,'k>','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% grid on
% box on
% set(gca,'FontSize',14)


% % figure(2)
plot3(q_L_1(1,:),q_L_1(2,:),q_L_1(3,:))
hold on
plot3(q_L_2(1,:),q_L_2(2,:),q_L_2(3,:))
% plot3(q_L(1,:),q_L(2,:),q_L(3,:))
% hold on
% %plot(q_NL1(1,:),q_NL1(2,:),'r--','LineWidth',1.5)
% xlabel('x (m)','FontSize',14) 
% ylabel('y (m)','FontSize',14) 
% legend('Linear EH','Non-Linear EH')
% plot(0,10,'ks','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% plot(20,10,'ks','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% plot(5,0,'ko','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% plot(15,0,'k>','MarkerSize', 8,'MarkerFaceColor','c')
% hold on
% grid on
% box on
% set(gca,'FontSize',14)

% save('/mnt/irisgpfs/users/htran/Fig518/HD_20GUsB10S_k30N802to3065to75rng10_C_size100to600_ratethreshold=2.5rng14Pk1.5P_U2.3.mat') 
% save('C:\Users\hieu.tran-dinh\OneDrive\PhD\Hieu\Papers\UAV_Backscatter\Code\V4\Fig3\Ph=1e7.mat')
toc