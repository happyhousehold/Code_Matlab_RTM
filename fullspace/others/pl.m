h=10;
rho=1;
s=0:0.001:2*pi;
plot(rho*cos(s),rho*sin(s)-h,'r','linewidth',1.5);
axis([-8 8 -15 1]);
hold on
x=-12:0.1:12;
plot(x,zeros(size(x)),'b','linewidth',1.5);
hold off