function y = plotconv(x,h,dt)

% plotconv: graphical visualization of the convolution
% y convolved signal
% x, h signals to convolve (1xn)
% dt sampling time

% the length of the convolved signal depends on the length of the originals
% signals
lx = length(x);
lh = length(h);
clen = lx+lh-1;

% new x axis
t=(-lh+1:clen-1)/dt;

y = conv(x,h); m=min(y); M=max(y);
% build up x on the new x axis
sig1 =  [zeros(1,lh-1),x,zeros(1,lh-1)];
sig2 =  [zeros(1,lh-1),h,zeros(1,lh-1)];
y = zeros(1,clen);

figure

subplot(5,1,1);
plot(t,sig2,'.-b');title('h(t)');
axis([t(1),t(end),min(sig2),max(sig2)]);

subplot(5,1,2);
plot(t,sig1,'.-b');title('x(t)');
axis([t(1),t(end),min(sig1),max(sig1)]);

for k=0:clen-1,
    sig3 = [zeros(1,k),fliplr(h),zeros(1,lx+lh-2-k)];

    subplot(5,1,3);
    plot(t,sig3,'.-b'); title('h(t-\tau)')
    axis([t(1),t(end),min(sig3),max(sig3)]);

    subplot(5,1,4);
    plot(t,sig1.*sig3,'.-b'); title('overlap'); 
    axis([t(1),t(end),min(sig1.*sig3)-.1,max(sig1.*sig3)+.1]);

    y(k+1) = sum(sig1.*sig3);
    subplot(5,1,5);
    plot(t,[zeros(1,lh-1),y],'.-r'); title('conv(x,h)');xlabel('ms');
    axis([t(1),t(end),m,M]);
    pause(0.5);
end

return
