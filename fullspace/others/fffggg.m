%% The plot of obstacle and sampling domain
figure,%circle
x=-3:0.1:3;h=0;
y=-3:0.1:3;
plot(x,(-3+h)*ones(size(x)),'b');
hold on
plot(x,(3+h)*ones(size(x)),'b');
hold on
plot(-3*ones(size(y)),y+h,'b');
hold on
plot(3*ones(size(y)),y+h,'b');
hold on
s=0:0.01:2*pi;
rho=1;
plot(rho*cos(s),rho*sin(s)+h,'r','linewidth',1.5);
axis([-18 18 -18 18]);
hold on
s=0:0.01:2*pi;
rho=10;
plot(rho*cos(s),rho*sin(s)+h,'r','linewidth',1.5);
axis([-18 18 -18 18]);hold on
s=0:0.01:2*pi;
rho=15;
plot(rho*cos(s),rho*sin(s)+h,'r','linewidth',1.5);
axis([-18 18 -18 18]);

figure,%kite
z=1.5:0.1:6.5;
plot(x,1.5*ones(size(x)),'b');
hold on
plot(x,6.5*ones(size(x)),'b');
hold on
plot(-3*ones(size(z)),z,'b');
hold on
plot(3*ones(size(z)),z,'b');
hold on
s=0:0.01:2*pi;
xs=cos(s)+0.65*cos(2*s)-0.65;  
ys=1.5*sin(s) ;
scale=1;
xs=scale*xs;
ys=scale*ys;
x_shift=0;y_shift=4;
xs=xs+x_shift;ys=ys+y_shift;
plot(xs,ys,'r','linewidth',1.5);
axis([-4 4 0 8]);


