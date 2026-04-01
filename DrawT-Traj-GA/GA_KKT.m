x = [20 30 40 50 60];
y1 = [265.519721 358.098904 466.523414 573.54726 682.255544];
y2 = [295.87813313219124 399.22116407373153 508.3566089810294 618.3961754308457 727.8545336334471];
figure(1)
hold on
plot(x, y1,'bo-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y2,'ms-','MarkerSize',8,'LineWidth',1.5);

xlabel('T (seconds)');
% ylabel('Total Throughput (Mbits)');
ylabel('Tổng thông lượng (Mbits)');
legend('BCD','GA');
grid on
box on
set(gca,'FontSize',14)