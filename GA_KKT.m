x = [20 30 40 50 60];
y1 = [282.9644 416.3363 550.4389 677.9835 806.8833];
% y2 = [197.214562 250.714046 265.81311 286.31753 327.176122];
y2 = [197.214562 250.714046 298.81311 348.31753 397.176122];
figure(1)
hold on
plot(x, y1,'bo-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y2,'ms-','MarkerSize',8,'LineWidth',1.5);

xlabel('T (seconds)');
ylabel('Total Throughput (Mbits)');
legend('GA','KKT');
grid on
box on
set(gca,'FontSize',14)