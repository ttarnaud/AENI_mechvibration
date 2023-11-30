%%
clear all;
close all;
%%
fs=10000; %Hz
t=0:1/fs:0.5-1/fs; %0.5 seconds

%frequency of different notes in Hz
Do_freq=262;
Re_freq=294;
Mi_freq=330;
Fa_freq=349;
Sol_freq=392;
La_freq=440;
Si_freq=494;
Do_freq_high=523;

%Sinusoid of 0.5s
Do=sin(2*pi*t*Do_freq); 
Re=sin(2*pi*t*Re_freq); 
Mi=sin(2*pi*t*Mi_freq); 
Fa=sin(2*pi*t*Fa_freq); 
Sol=sin(2*pi*t*Sol_freq); 
La=sin(2*pi*t*La_freq); 
Si=sin(2*pi*t*Si_freq); 
Do_high=sin(2*pi*t*Do_freq_high); 

%%
%plotting first 1000 points
figure;
subplot(8,1,1);plot(t(1:1000),Do(1:1000));
subplot(8,1,2);plot(t(1:1000),Re(1:1000));
subplot(8,1,3);plot(t(1:1000),Mi(1:1000));
subplot(8,1,4);plot(t(1:1000),Fa(1:1000));
subplot(8,1,5);plot(t(1:1000),Sol(1:1000));
subplot(8,1,6);plot(t(1:1000),La(1:1000));
subplot(8,1,7);plot(t(1:1000),Si(1:1000));
subplot(8,1,8);plot(t(1:1000),Do_high(1:1000));
%%
%play music scale
music_scale=[Do Re Mi Fa Sol La Si Do_high];
sound(music_scale,fs);
%%
sound([music_scale fliplr(music_scale)],fs)
%%
%brother jacob
a=La;
b=Si;
c=Do;
d=Re;
e=Mi;
f=Fa;
g=Sol;
g_low=sin(2*pi*t*Sol_freq/2);
brother_jacob=[c d e c zeros(1,500) ...
    c d e c zeros(1,500) ...
    e f g zeros(1,5000) ...
    e f g zeros(1,5000) ...
    g(1:end/2) a(1:end/2) g(1:end/2) f(1:end/2) e c(1:end/2) zeros(1,2500)...
    g(1:end/2) a(1:end/2) g(1:end/2) f(1:end/2)  e c(1:end/2) zeros(1,2500)...
    c g_low c zeros(1,5000)...
    c g_low c];
time=0:1/fs:length(brother_jacob)/fs-1/fs;
figure,plot(time,brother_jacob);
%%
%play song
sound(brother_jacob,fs)


%%
%white noise generation
y=brother_jacob;
noise=randn(size(y));
ynoise=y+noise;

time_points=1000;
figure;
subplot(3,1,1); plot(time(1:time_points),y(1:time_points));
subplot(3,1,2); plot(time(1:time_points),noise(1:time_points));
subplot(3,1,3); plot(time(1:time_points),ynoise(1:time_points));
%%
%play with noise
sound(ynoise,fs)

%%
%filter noise using convolution with kernel
k=ones(1,5)/5;
yfilt=conv(k,ynoise);
figure;
subplot(3,1,1); plot(time(1:time_points),y(1:time_points));ylim([-2 2]);
subplot(3,1,2); plot(time(1:time_points),ynoise(1:time_points));ylim([-2 2]);
subplot(3,1,3); plot(time(1:time_points),yfilt(1:time_points));ylim([-2 2]);

%%
%play filtered sound 
sound(yfilt,fs)
%%
%Fourier Transformation
fs=200;
t=0:1/fs:0.5-1/fs;
f=40.45;
x=sin(2*pi*f*t);

Fx=fft(x);
N=length(x);
freq= 0 : fs/N : fs/2-fs/N;
abs_Fx=abs(Fx);

figure,
plot(freq,abs_Fx(1:length(freq)));
title('Fourier Spectrum')
xlabel('freq (Hz)');
ylabel('Fx');

%%
%fft with higher frequency resolution
NFFT=2^12;
Fx_longer=fft(x,NFFT);
f=0:fs/NFFT:fs/2-fs/NFFT;

figure;
plot(f,abs(Fx_longer(1:length(f))));
title('Fourier Spectrum')
xlabel('freq (Hz)');
ylabel('Fx');
%%
% full spectrum
f=0:fs/NFFT:fs-fs/NFFT;
figure;
plot(f,abs(Fx_longer(1:length(f))));
title('Fourier Spectrum')
xlabel('freq (Hz)');
ylabel('Fx');
%%
% negative side before positive side
f_pos=0:fs/NFFT:fs/2-fs/NFFT;
f_neg=-fs/NFFT:-fs/NFFT:-fs/2;
f=[fliplr(f_neg) f_pos];
Fx_shift=fftshift(Fx_longer);
figure;
plot(f,abs(Fx_shift));
title('Fourier Spectrum')
xlabel('freq (Hz)');
ylabel('Fx');
%%
% fourier transform of brother jacob
fs=10000; % Reset to original sampling rate
x=brother_jacob;
NFFT=2^(nextpow2(size(x,2)));
Fx_longer=fft(x,NFFT);
f=0:fs/NFFT:fs/2-fs/NFFT;
 
figure;
plot(f,abs(Fx_longer(1:length(f))));
title('Fourier Spectrum')
xlabel('freq (Hz)');
ylabel('Fx');
%%
% fourier transform of brother jacob with noise
x=ynoise;
NFFT=2^(nextpow2(size(x,2)));
Fx_longer=fft(x,NFFT);
f=0:fs/NFFT:fs/2-fs/NFFT;
 
figure;
plot(f,abs(Fx_longer(1:length(f))));
title('Fourier Spectrum')
xlabel('freq (Hz)');
ylabel('Fx');

%%
% fourier transform of filtered brother jacob with noise
x=yfilt;
NFFT=2^(nextpow2(size(x,2)));
Fx_longer=fft(x,NFFT);
f=0:fs/NFFT:fs/2-fs/NFFT;
 
figure;
plot(f,abs(Fx_longer(1:length(f))));
title('Fourier Spectrum')
xlabel('freq (Hz)');
ylabel('Fx');

