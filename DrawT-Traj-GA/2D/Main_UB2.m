%% 
clc; clear all; close all;
tic
addpath('D:\Program Files\Matlab OPtimization\cvx');
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
% sigma = 0.9; % caching gain coefficience: part of file that UB cache %???
sigma = 0.45;
delta_t = 0.5; % second
Euler = 0.5772156649; % Euler constant 8.367.1 in Table of Integral
q_I1 = [0;10];q_F1 = [20;10];
q_I2 = [0;5];q_F2 = [20;5]; %???
w_s = [5;0]; w_d = [15;0];

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
    %% Initialize trajectory q
    q_e1 = [0.5*(w_s(1,:)+w_d(1,:));0];
     q_j_1 = [];
     q_start = q_I1;
     N_segment = (N-1)/2;
     m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:));
     center = q_start(2,:) - q_start(1,:) * m;
     q_j_1= [];
     for x=q_start(1,:):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:)
        y = m *x + center;
        q_j_1 = [q_j_1,[x;y]];
     end
    q_j_1 = [q_j_1,q_e1];
    %% 
     q_start = [0.5*(w_s(1,:)+w_d(1,:));0];

     q_e1 = q_F1;
     N_segment = (N-1)/2;
     m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:));
     center = q_start(2,:) - q_start(1,:) * m;
     for x=q_start(1,:)+(q_e1(1,:)-q_start(1,:))/(N_segment):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:)
        y = m *x + center;
        q_j_1 = [q_j_1,[x;y]];
     end 
    q_j_1 = [q_j_1,q_e1];
    %% 
    q_e1 = [0.5*(w_s(1,:)+w_d(1,:));0];
     q_j_2 = [];
     q_start = q_I2;
     N_segment = (N-1)/2;
     m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:));
     center = q_start(2,:) - q_start(1,:) * m;
     q_j_2= [];
     for x=q_start(1,:):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:)
        y = m *x + center;
        q_j_2 = [q_j_2,[x;y]];
     end
    q_j_2 = [q_j_2,q_e1];
    %% 
     q_start = [0.5*(w_s(1,:)+w_d(1,:));0];

     q_e1 = q_F2;
     N_segment = (N-1)/2;
     m = (q_e1(2,:)-q_start(2,:)) / (q_e1(1,:)-q_start(1,:));
     center = q_start(2,:) - q_start(1,:) * m;
     for x=q_start(1,:)+(q_e1(1,:)-q_start(1,:))/(N_segment):(q_e1(1,:)-q_start(1,:))/(N_segment):q_e1(1,:)
        y = m *x + center;
        q_j_2 = [q_j_2,[x;y]];
     end 
    q_j_2 = [q_j_2,q_e1];
    %% Test 
%     plot(q_j_1(1,:),q_j_1(2,:))
%     hold on
%     plot(q_j_2(1,:),q_j_2(2,:)) %???
    %% Initial backscatter coefficient and DTS
    tau_j = 0.5*ones(1,N);
    eta_j = eta_max*ones(1,N);
    %% Linear EH  model: Alternating optimization algorithm
    [q_L_1,tau_L_1,eta_L_1,g_x_L_1,P_un_1,P_sn_1] = Linear_model2(tau_j,eta_j,q_j_1,N);
%     [q_L_2,tau_L_2,eta_L_2,g_x_L_2,P_un_2,P_sn_2] = Linear_model3(tau_j,eta_j,q_j_2,N);
end


fprintf('Optimal value of KKT_base for UAV_1 = %f\n', g_x_L_1);
% fprintf('Optimal value of KKT_base for UAV_2 = %f\n', g_x_L_2);
%% Test plot Trajectory
figure(1)
plot(q_L_1(1,:),q_L_1(2,:),'b-','LineWidth',1.5)
% hold on
% plot(q_L_2(1,:),q_L_2(2,:),'b-','LineWidth',1.5)
toc