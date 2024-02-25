%% Test plot Trajectory
clear; clc; close all;
load('Ph=1e7_2.mat')
figure(1)
plot3(q_L_11e7(1,:),q_L_11e7(2,:),q_L_11e7(3,:),'b-','LineWidth',1)
hold on
plot3(q_L_21e7(1,:),q_L_21e7(2,:),q_L_21e7(3,:),'r-','LineWidth',1)
hold on
scatter3(15, 0, 0, 'filled', 'Marker', '^', 'DisplayName', 'My Point')
hold on
scatter3(5, 0, 0, 'filled', 'Marker', 'd', 'DisplayName', 'My Point')
hold on
scatter3(0, 10, 10,'filled', 'Marker', 'o', 'DisplayName', 'My Point')
hold on
scatter3(20, 10, 10,'filled', 'Marker', 'o', 'DisplayName', 'My Point')
hold on
scatter3(0, 10, 5,'filled', 'Marker', 's', 'DisplayName', 'My Point')
hold on
scatter3(20, 10, 5,'filled', 'Marker', 's', 'DisplayName', 'My Point')

xlabel('x (m)','FontSize',14);
ylabel('y (m)','FontSize',14);
zlabel('z (m)','FontSize',14);
legend('UAV_1','UAV_2', 'Destination', 'Source', 'Initial point of UAV_1', 'Final point of UAV_1', 'Initial point of UAV_2', 'Final point of UAV_2');
grid on
box on
set(gca,'FontSize',14)