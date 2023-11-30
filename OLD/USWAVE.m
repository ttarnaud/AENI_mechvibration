L=0.07;
x=-L/2:0.0001:L/2;
c=1550;
f1 = 1e6;
f2 = 1.1e6;
t=0;
%%
lambda1=c/f1;
lambda2 = c/f1;
k1=2*pi/lambda1;
k2=2*pi/lambda1;
w1=2*pi*f1;
w2=2*pi*f1;
I = 10*10^3;
Z=1.6*10^6;
v0 = sqrt(2*I/Z);
p0 = v0*Z;
p1=@(x,t) p0*cos(k1*x-w1*t);
p2 =@(x,t) p0*cos(k1*x+w1*t);
v1 =@(x,t) v0*cos(k1*x-w1*t);
v2 =@(x,t) -v0*cos(k1*x+w1*t);
xpointsv = -9*lambda1/4:2*lambda1:9*lambda1/4;
xpointsp = -4*lambda1/2:2*lambda1:4*lambda1/2;
xpoints = sort([xpointsv,xpointsp]);
p =@(x,t) p1(x,t)+p2(x,t);
v=@(x,t) v1(x,t)+v2(x,t);
t = linspace(0,2/f1,1000);
figure
for i=1:length(xpoints)
    plot3(xpoints(i).*ones(size(t)),t,p(xpoints(i),t)/p0,'b')
    hold on
    plot3(xpoints(i).*ones(size(t)),t,v(xpoints(i),t)/v0,'r')
end
figure(10)
subplot(4,1,1)
plot3(x,zeros(size(x)),p(x,0)/p0/2)
hold on
plot3(x,1/(4*f1)*ones(size(x)),v(x,1/(4*f1))/v0/2)
hold off
legend('Pressure','Particle velocity or displacement')
ylabel('time')
xlabel('position')
zlabel('Normalized amplitudes')
set(gca,'view',[0,0])
set(findobj('type','axes'),'fontsize',14)
set(gca,'XAxisLocation','origin')
%%
lambda1=c/f1;
lambda2 = c/f2;
k1=2*pi/lambda1;
k2=2*pi/lambda2;
w1=2*pi*f1;
w2=2*pi*f2;
I = 10*10^3;
Z=1.6*10^6;
v0 = sqrt(2*I/Z);
p0 = v0*Z;
p1=@(x,t) p0*cos(k1*x-w1*t);
p2 =@(x,t) p0*cos(k2*x+w2*t);
v1 =@(x,t) v0*cos(k1*x-w1*t);
v2 =@(x,t) -v0*cos(k2*x+w2*t);
xpointsv = -9*lambda2/4:2*lambda2:9*lambda2/4;
xpointsp = -4*lambda2/2:2*lambda2:4*lambda2/2;
xpoints = sort([xpointsv,xpointsp]);
p =@(x,t) p1(x,t)+p2(x,t);
v=@(x,t) v1(x,t)+v2(x,t);
t = linspace(0,2/f2,1000);
figure
for i=1:length(xpoints)
    plot3(xpoints(i).*ones(size(t)),t,p(xpoints(i),t)/p0,'b')
    hold on
    plot3(xpoints(i).*ones(size(t)),t,v(xpoints(i),t)/v0,'r')
end
figure(10)
subplot(4,1,2)
plot3(x,zeros(size(x)),p(x,0)/p0/2)
hold on
plot3(x,1/(4*f1)*ones(size(x)),v(x,1/(4*f1))/v0/2)
hold off
legend('Pressure','Particle velocity or displacement')
ylabel('time')
xlabel('position')
zlabel('Normalized amplitudes')
set(gca,'view',[0,0])
set(findobj('type','axes'),'fontsize',14)
set(gca,'XAxisLocation','origin')
for it = 1:length(t)
    if it==1
        pmax = p(x,t(it))/p0/2;
        vmax = v(x,t(it))/v0/2;
    else
        pmax = max([pmax;p(x,t(it))/p0/2],[],1);
        vmax = max([vmax;v(x,t(it))/v0/2],[],1);
    end
end
figure(10)
subplot(4,1,3)
plot3(x,zeros(size(x)),pmax)
hold on
plot3(x,1/(4*f1)*ones(size(x)),vmax)
hold off
legend('Pressure','Particle velocity or displacement')
ylabel('time')
xlabel('position')
zlabel('Normalized amplitudes')
set(gca,'view',[0,0])
set(findobj('type','axes'),'fontsize',14)
set(gca,'XAxisLocation','origin')
%%
lambda1=c/f1;
lambda2 = c/f2;
k1=2*pi/lambda1;
k2=2*pi/lambda2;
w1=2*pi*f1;
w2=2*pi*f2;
I = 10*10^3;
Z=1.6*10^6;
v0 = sqrt(2*I/Z);
p0 = v0*Z;
p1=@(x,t) p0*cos(k1*x-w1*t);
p2 =@(x,t) p0*cos(k1*x+w1*t);
v1 =@(x,t) v0*cos(k1*x-w1*t);
v2 =@(x,t) -v0*cos(k1*x+w1*t);
p3=@(x,t) p0*cos(k2*x-w2*t);
p4 =@(x,t) p0*cos(k2*x+w2*t);
v3 =@(x,t) v0*cos(k2*x-w2*t);
v4 =@(x,t) -v0*cos(k2*x+w2*t);
xpointsv = -9*lambda2/4:2*lambda2:9*lambda2/4;
xpointsp = -4*lambda2/2:2*lambda2:4*lambda2/2;
xpoints = sort([xpointsv,xpointsp]);
p =@(x,t) p1(x,t)+p2(x,t)+p3(x,t)+p4(x,t);
v=@(x,t) v1(x,t)+v2(x,t)+v3(x,t)+v4(x,t);
t = linspace(0,11/f2,1e5);
figure
for i=1:length(xpoints)
    plot3(xpoints(i).*ones(size(t)),t,p(xpoints(i),t)/p0,'b')
    hold on
    plot3(xpoints(i).*ones(size(t)),t,v(xpoints(i),t)/v0,'r')
end
figure(10)
subplot(4,1,4)
plot3(x,zeros(size(x)),p(x,0)/p0/2)
hold on
plot3(x,1/(4*f1)*ones(size(x)),v(x,1/(4*f1))/v0/2)
hold off
legend('Pressure','Particle velocity or displacement')
ylabel('time')
xlabel('position')
zlabel('Normalized amplitudes')
set(gca,'view',[0,0])
set(findobj('type','axes'),'fontsize',14)
set(gca,'XAxisLocation','origin')