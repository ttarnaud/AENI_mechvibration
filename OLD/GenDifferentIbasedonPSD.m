% Generate 1e6 different signals
clear all
close all
clc
addpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions\noisefunctions')
% init signal
Tend=0.029;
Fs = 2^18/Tend;
t = 0:1/Fs:Tend-1/Fs;
s0 = Itime(t,500e-6);
% calculate PSD
N = length(s0);
s0dft = fft(s0);
s0dftr = s0dft(1:N/2+1);
psds0 = (1/(N)) * abs(s0dftr).^2;
psds0(2:end-1) = 2*psds0(2:end-1);
freq = 0:Fs/N:Fs/2;
filterdt = 100*1e-6;
windowSize = filterdt*Fs;
b = (1/windowSize)*ones(1,int16(windowSize));
if input('calculate new signals? \n')
% calculate new signals
tic
parfor i=1:100
    stemp = [];
    for j=1:2
    Phifk = 2*pi*rand(1,N);
    Zfk = s0dft.*exp(complex(0,Phifk));
    stemp=horzcat(stemp,ifft(Zfk,'symmetric'));
    end
    %signals(i,:) = stemp;
    signals(i,:) = filter(b,1,stemp);
    disp(num2str(i))
end
toc
t2 = [t,t+t(end)];
disp('saving')
%save('100randomIlength58ms.mat','signals','t2','-v7.3')
end
clear('signals')
disp('start second')
tic

for i=1:4
    stemp = [];
    Phifk = 2*pi*rand(1,N);
    Zfk = s0dft.*exp(complex(0,Phifk));
    %stemp=horzcat(stemp,ifft(Zfk,'symmetric'));
    signals(i,:) = ifft(Zfk,'symmetric');
    %signals(i,:) = filter(b,1,signals(i,:));
    disp(num2str(i))
end
toc
if input('save?\n')
disp('saving')
save('200randomIlength29ms.mat','signals','t','-v7.3')
end
% generate 'alpha' noise
an = alphanoise(1, N,1);
an = an./max(an);
% alpha funciton
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0;
Tau = 0.005;
AlphaDelay = 0;
Ifun = Alphafun(t,Tau,AlphaDelay);

%% plot 
figure
subplot(1,2,1)
plot(freq,10*log10(psds0),'b')
hold on
s2=signals(1,:);
[ps2,fs2]=periodogram(s2,rectwin(length(s2)),length(s2),Fs,'power');
plot(fs2,10*log10(ps2),'g')
[Pxx, fan] = periodogram(an, rectwin(length(an)),length(an), Fs,'power');
plot(fan,10*log10(Pxx),'r')
[Pxxaf, faf] = periodogram(Ifun, rectwin(length(Ifun)),length(Ifun), Fs,'power');
plot(faf,10*log10(Pxxaf),'y')
hold off
grid on
set(gca,'xscale','log')
title('Periodogram')
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')
subplot(1,2,2)
plot(t,an,'r','DisplayName','pink noise (1/f)')
hold on
t2 = [t,t+t(end)];
plot(t,s2,'g','DisplayName','modified ECD CA1');
plot(t,s0,'b','DisplayName','source ECD CA1')
plot(t,Ifun,'y','DisplayName','alpha function')
hold off
legend('show')
% pwelch periodogram
figure
winlen = floor(N/10);
window = hanning(winlen, 'periodic');
noverlap = winlen/2;
nfft = winlen;
[Pxx, fan] = pwelch(an, window, noverlap, nfft, Fs, 'onesided');
subplot(2,1,1)
[ps0,fs0]=pwelch(s0, window, noverlap, nfft, Fs, 'onesided');
plot(fs0,10*log10(ps0),'b')
hold on
s2=signals(1,:);
[ps2,fs2]=pwelch(s2, window, noverlap, nfft, Fs, 'onesided');
plot(fs2,10*log10(ps2),'g')
plot(fan,10*log10(Pxx),'r')
hold off
grid on
set(gca,'xscale','log')
title('pwelch')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
subplot(2,1,2)
plot(t,an,'r')

hold on
t2 = [t,t+t(end)];
plot(t,s2,'g');
plot(t,s0,'b')
hold off