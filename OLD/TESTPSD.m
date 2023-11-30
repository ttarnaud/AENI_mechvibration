clear all
close all
clc
% init signal
Tend=0.029;
Fs = 2^21/Tend;
t = 0:1/Fs:Tend-1/Fs;
s0 = Itime(t,500e-6);
figure
plot(t*1e3,s0)
xlabel('time [ms]')
ylabel('normalized dipole moment')
%% calculate PSD
N = length(s0);
s0dft = fft(s0);
s0dftr = s0dft(1:N/2+1);
psds0 = (1/(N)) * abs(s0dftr).^2;
psds0(2:end-1) = 2*psds0(2:end-1);
freq = 0:Fs/N:Fs/2;
% calculate new signals
Apsd = sqrt(N*psds0/2);
Phifk = 2*pi*rand(1,length(Apsd));
Zfk = Apsd.*exp(complex(0,Phifk));
Yr = Zfk;
Yl=Yr(end-1:-1:2);
Yl=conj(Yl);
Y2=[Yr Yl];
s1=ifft(Y2,'symmetric');
Phifk = 2*pi*rand(1,N);
Zfk = s0dft.*exp(complex(0,Phifk));
s2=ifft(Zfk,'symmetric');

subplot(2,1,1)
plot(freq/10^6,10*log10(1/Fs*psds0),'b')
hold on
[ps1,fs1]=periodogram(s1,rectwin(length(s1)),length(s1),Fs);
plot(freq/10^6,10*log10(ps1),'r')
[ps2,fs2]=periodogram(s2,rectwin(length(s2)),length(s2),Fs);
plot(freq/10^6,10*log10(ps2),'g')
hold off
grid on
title('Periodogram Using FFT')
xlabel('Frequency (MHz)')
ylabel('Power (dB)')
subplot(2,1,2)
plot(t,s0)
hold on
plot(t,s1./max(abs(s1)))
plot(t,s2./max(abs(s2)));
hold off
ylabel('normalized dipole moment')
xlabel('time [s]')
%% not same result when calculated form Rxx ? :( periodogram does give the same
% Rs0s0 = xcorr(s0,'biased');
% lag = (-(N-1):1:(N-1))*1/Fs;
% 
% PSDfromRxx = 1/Fs*abs(fft(Rs0s0,N)/N);
% PSDfromRxx = PSDfromRxx(1:N/2+1);
% PSDfromRxx(2:end-1) = PSDfromRxx(2:end-1)*2;
% plot(freq/10^6,10*log10(PSDfromRxx(1:length(freq))));
% 
% periodogram(s0,rectwin(length(s0)),length(s0),Fs)
hold off