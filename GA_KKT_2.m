x = [20 30 40 50 60];
y1 = [1230 1560 1800 2040 2340];
% y2 = [359.343585 761.6494728 617.942676 859.266496];
y2 = [646.818453 810.2494728 1112.2968168 1546.6796928 1714.9827885];
figure(1)
hold on
plot(x, y1,'bo-','MarkerSize',8,'LineWidth',1.5);
hold on
plot(x, y2,'ms-','MarkerSize',8,'LineWidth',1.5);

xlabel('T (seconds)');
ylabel('Running time (seconds)');
legend('GA','KKT');
grid on
box on
set(gca,'FontSize',14)