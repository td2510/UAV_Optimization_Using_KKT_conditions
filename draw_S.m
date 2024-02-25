x = [20 30 40 50 60 70];
y1 = [65.350058 81.747292 94.149883 116.185832 132.484 146.327];
y2 = [34.917570 46.139436 53.314697 72.248670 75.322956 90.915282];
y3 = [49.6616 64.5907 78.4407 89.498841 97.871520 105.359];
% y1 = [63.027675 59.598859 72.585509 73.453368 90.123698 98.607281];
% y2 = [34.917570 46.139436 53.314697 72.248670 75.322956 90.915282];
% y3 = [115.144253 115.144253 115.144253 115.144253 115.144253 115.144253];
figure(1)
hold on
plot(x, y1,'bo-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y2,'ms-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y3,'gd-','MarkerSize',8,'LineWidth',1.5);

xlabel('S (Mbits)');
ylabel('Total Throughput (Mbits)');
legend('Com','3D+OP','3D+2UAV')
grid on
box on
set(gca,'FontSize',14)
