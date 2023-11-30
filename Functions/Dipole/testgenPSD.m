%test GenPSD
close all
fus = 1;
Fs = 15
t = 0:1/Fs:1-1/Fs;
t(end)
N = length(t)
signal = sin(2*pi*fus*t);
figure
subplot(2,2,1)
plot(t,signal)

[PSD,f] = genPSD(signal,t,Fs,0);
subplot(2,2,2)
plot(f,PSD)


N2 = 17;
t = linspace(0,1,N2);
t(end)
Fs2 = 1/(t(2)-t(1))
signal = sin(2*pi*fus*t);
subplot(2,2,3)
plot(t,signal)

[PSD,f] = genPSD(signal,t,Fs2,0);
subplot(2,2,4)
plot(f,PSD)