x = [20 40 60 80 100];
y1 = [98.607281 132.906555 143.588061 154.114156 180.522864];
y2 = [90.915282 91.252732 94.508997 94.721053 96.601535];
y3 = [71.523535 80.026131 97.006876 103.141840 119.541905];
figure(1)
hold on
plot(x, y1,'bo-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y2,'ms-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y3,'gd-','MarkerSize',8,'LineWidth',1.5);

xlabel('T (seconds)');
ylabel('Total Throughput (Mbits)');
legend('Com','3D+OP','3D+2UAV');
grid on
box on
set(gca,'FontSize',14)