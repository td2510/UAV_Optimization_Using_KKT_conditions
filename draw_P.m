x = [70 72 74 76 78 80];
y1 = [98.607281 104.988008 110.000903 113.611489 115.351 116.844939];
y2 = [90.915282 93.188231 95.101188 96.664959 97.395736 98.115282];
y3 = [71.523535 73.083561 74.640504 74.826430 75.166803 75.195808];
figure(1)
hold on
plot(x, y1,'bo-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y2,'ms-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y3,'gd-','MarkerSize',8,'LineWidth',1.5);

xlabel('P_{WPT} (dB)');
ylabel('Total Throughput (Mbits)');
legend('Com','3D+OP','3D+2UAV')
grid on
box on
set(gca,'FontSize',14)