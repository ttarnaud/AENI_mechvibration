dt = 0.1; %fs = 10Hz

%define time axis
t = 0:dt:4;

% rect 
A = 3;
w = 3;
x = A*rect(t,0,w-dt);
%x= A*sin(2*pi*t*1);

% exponential
T0 = 1;
h = exp(-t/T0);
%h = rect(t,0,T0-dt);

y = plotconv(x,h,dt);