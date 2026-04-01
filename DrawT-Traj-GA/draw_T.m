x = [20 30 40 50 60 70];
y1 = [265.519721 358.098904 466.523414 583.54726 711.69559 823.875742];
y2 = [162.363598 243.92355 325.27938 406.546632 487.796664 569.051344];
y3 = [150.429209 224.917392 297.857950 384.985574 463.791331 543.710646];
% y4 = [232.210036 350.091738 454.185028 572.022912 700.179668 816.884022];


figure(2)
hold on
plot(x, y1,'bo-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y2,'ms-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y3,'gd-','MarkerSize',8,'LineWidth',1.5);
% hold on
% plot(x, y4,'gd-','MarkerSize',8,'LineWidth',1.5);

xlabel('T (seconds)');
% ylabel('Total Throughput (Mbits)');
ylabel('Tổng thông lượng (Mbits)');
legend('Com','3D+OP','2D+2UAV');
% legend('Com','3D+OP','2D+2UAV','3D+2UAV');
grid on
box on
set(gca,'FontSize',14)