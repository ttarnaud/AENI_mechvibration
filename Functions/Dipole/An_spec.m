clear all
close all
Tau = 0.005;
AlphaDelay = 0;
fus = 1e3;
Fs1 = 40e3;
Fs2 = 40e6;
fun = @(t) sin(2*pi*fus*t);
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0;

t = 0:1/Fs1:1;
t2 = 0:1/Fs2:1;

sininit = fun(t);
sin_best = fun(t2);
sin_interp = interp1(t,sininit,t2);
%%
%sininit = sininit(1:end-1);
% [sin_p1,sin_f1] = periodogram(sininit,[],[],40e3,'power');
% [sin_p2,sin_f2] = periodogram(sin_best,rectwin(length(sin_best)),length(sin_best),40e6,'power');
% [sin_p3,sin_f3] = periodogram(sin_interp,rectwin(length(sin_interp)),length(sin_interp),40e6,'power');

[sin_p1,sin_f1] = genPSD(sininit,t,Fs1,'db/freq');
[sin_p2,sin_f2] = genPSD(sin_best,t2,Fs2,'db/freq');
[sin_p3,sin_f3] = genPSD(sin_interp,t2,Fs2,'db/freq');

Alphainit = Alphafun(t,Tau,AlphaDelay);
Alpha_best = Alphafun(t2,Tau,AlphaDelay);
Alpha_interp = interp1(t,Alphainit,t2);

% Alphainit = Alphainit(1:end-1);

% [A_p1,A_f1] = periodogram(Alphainit,rectwin(length(Alphainit)),length(Alphainit),40e3,'power');
% [A_p2,A_f2] = periodogram(Alpha_best,rectwin(length(Alpha_best)),length(Alpha_best),40e6,'power');
% [A_p3,A_f3] = periodogram(Alpha_interp,rectwin(length(Alpha_interp)),length(Alpha_interp),40e6,'power');
[A_p1,A_f1] = genPSD(Alphainit,t,Fs1,'db/freq');
[A_p2,A_f2] = genPSD(Alpha_best,t2,Fs2,'db/freq');
[A_p3,A_f3] = genPSD(Alpha_interp,t2,Fs2,'db/freq');



%%
figure
plot(sin_f1,sin_p1,'displayName','Fs = 40 kHz')
hold on
plot(sin_f3,sin_p3,'displayName','40 kHz interpolated to 40 MHz')
plot(sin_f2,sin_p2,'displayName','Fs = 40 MHz')
set(gca,'xscale','log')
legend('show')
title('sinus')
hold off
ylabel('power/frequency  (DB/FS)')

figure
plot(A_f1,A_p1,'displayName','Fs = 40 kHz')
hold on

plot(A_f3,A_p3,'displayName','40 kHz interpolated to 40 MHz')
plot(A_f2,A_p2,'displayName','Fs = 40 MHz')
set(gca,'xscale','log')
legend('show')
hold off

% normalized to 1
figure
plot(sin_f1,sin_p1-max(sin_p1),'displayName','Fs = 40 kHz')
hold on
plot(sin_f3,sin_p3-max(sin_p1),'displayName','40 kHz interpolated to 40 MHz')
plot(sin_f2,sin_p2-max(sin_p1),'displayName','Fs = 40 MHz')
set(gca,'xscale','log')
legend('show')
title('sinus')
hold off
ylabel('power/frequency  (DB/FS)')

figure
plot(A_f1,A_p1-max(A_p1),'displayName','Fs = 40 kHz')
hold on

plot(A_f3,A_p3-max(A_p1),'displayName','40 kHz interpolated to 40 MHz')
plot(A_f2,A_p2-max(A_p1),'displayName','Fs = 40 MHz')
set(gca,'xscale','log')
legend('show')
hold off
%%
figure 
plot(A_f2,A_p2-A_p3,'displayName','Fs = 40 MHz')
set(gca,'xscale','log','yscale','log')

% 
%% noise function
noiseinit = Itime(t,0.5e3,'hipPyr25us1');
noise2 = Itime(t2,0.5e3,'hipPyr25us1');
noise_interp = interp1(t,noiseinit,t2);
[N_p1,N_f1] = genPSD(noiseinit,t,Fs1,'db');
[N_p2,N_f2] = genPSD(noise2,t2,Fs2,'db');
[N_p3,N_f3] = genPSD(noise_interp,t2,Fs2,'db');

%%
figure
plot(N_f1,N_p1,'displayName','Fs = 40 kHz')
hold on
plot(N_f3,N_p3,'displayName','40 kHz interpolated to 40 MHz')
plot(N_f2,N_p2,'displayName','Fs = 40 MHz')
set(gca,'xscale','log')
legend('show')
title('NoiseDM')
hold off
ylabel('power/frequency  (DB/FS)')



% normalized to 1
figure
plot(N_f2,N_p2-max(N_p2),'displayName','Fs = 40 MHz')
hold on
plot(N_f3,N_p3-max(N_p3),'displayName','40 kHz interpolated to 40 MHz')

plot(N_f1,N_p1-max(N_p1),'displayName','Fs = 40 kHz')
set(gca,'xscale','log')
legend('show')
title('NoiseDM')
hold off
ylabel('power/frequency  (DB/FS)')
%we see high powerlaw trend but still highpower for 10^6. largely due to
%interpolation> Ofcourse frequency content at this level is not known thus
%speculation of real value. According to literature at low frequency domain
%a small powerlaw of alpha = 1-2 is followed with possible bumb at 10^1. 
%with our generated signals we excpect a higher powerlaw when f=10^3 is
%exceeded. 

%% what is effect of filter in getCurrentArray

Tend = 0.1;
Tsim = t;
% init signal
nTend = ceil(Tsim(end)/Tend);
Tend = nTend*Tend;
Fs = 2^(17+nextpow2(nTend))/Tend; %this way Fs is at least 4MHz
tt = 0:1/Fs:Tend;
s0 = Itime(tt,0.5e-3,'hipPyr25us1');  %upsampled signal with linear interp. create high freaquency components are questionqble 1/f relationship
[N_p1,N_f1] = genPSD(s0,tt,Fs,'db');
% calculate PSD
L = length(s0);
PSD = fft(s0);
filterdt = 100*1e-6;
windowSize = filterdt*Fs;
b = (1/windowSize)*ones(1,int16(windowSize));
Phifk = 2*pi*rand(1,L);   %random phases
Zfk = PSD.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
Iarray_temp = ifft(Zfk,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
Iarray_temp = Iarray_temp./max(abs(Iarray_temp));
[N_p4,N_f4] = genPSD(Iarray_temp,tt,Fs,'db');
Iarray = interp1(tt,Iarray_temp,Tsim);
[N_p2,N_f2] = genPSD(Iarray,Tsim,Fs1,'db');
Iarray2 = filter(b,1,Iarray);
Iarray2 = Iarray2./max(abs(Iarray2));
[N_p3,N_f3] = genPSD(Iarray2,Tsim,Fs1,'db');

figure
subplot(2,4,[1,4])
plot(tt,s0,'displayName','raw signal')
hold on
plot(tt,Iarray_temp,'displayName','Iarray_temp')
plot(Tsim,Iarray,'displayName','Iarray')
plot(Tsim,Iarray2,'displayName','Iarray2')
hold off
legend('show')
subplot(2,4,5)
plot(N_f2,N_p2-max(N_p2),'displayName','Iarray')
set(gca,'xscale','log')
legend('show')
subplot(2,4,6)
plot(N_f3,N_p3-max(N_p3),'displayName','Iarray2')
set(gca,'xscale','log')
legend('show')
subplot(2,4,7)
plot(N_f1,N_p1-max(N_p1),'displayName','s0')
set(gca,'xscale','log')
legend('show')
subplot(2,4,8)
plot(N_f4,N_p4-max(N_p4),'displayName','Iarray_temp')
set(gca,'xscale','log')
legend('show')
title('NoiseDM')
hold off
ylabel('power/frequency  (DB/FS)')

%% higher Fs t2

Tend = 0.1;
Tsim = t2;
% init signal
nTend = ceil(Tsim(end)/Tend);
Tend = nTend*Tend;
Fs = 2^(17+nextpow2(nTend))/Tend; %this way Fs is at least 4MHz
tt = 0:1/Fs:Tend;
s0 = Itime(tt,0.5e-3,'hipPyr25us1');  %upsampled signal with linear interp. create high freaquency components are questionqble 1/f relationship
[N_p1,N_f1] = genPSD(s0,tt,Fs,'db');
% calculate PSD
L = length(s0);
PSD = fft(s0);
filterdt = 100*1e-6;
windowSize = filterdt*Fs;
b = (1/windowSize)*ones(1,int16(windowSize));
Phifk = 2*pi*rand(1,L);   %random phases
Zfk = PSD.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
Iarray_temp = ifft(Zfk,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
Iarray_temp = Iarray_temp./max(abs(Iarray_temp));
[N_p4,N_f4] = genPSD(Iarray_temp,tt,Fs,'db');
Iarray = interp1(tt,Iarray_temp,Tsim);
[N_p2,N_f2] = genPSD(Iarray,Tsim,Fs2,'db');
Iarray2 = filter(b,1,Iarray);
Iarray2 = Iarray2./max(abs(Iarray2));
[N_p3,N_f3] = genPSD(Iarray2,Tsim,Fs2,'db');

figure
subplot(2,4,[1,4])
plot(tt,s0,'displayName','raw signal')
hold on
plot(tt,Iarray_temp,'displayName','Iarray_temp')
plot(Tsim,Iarray,'displayName','Iarray')
plot(Tsim,Iarray2,'displayName','Iarray2')
hold off
legend('show')
subplot(2,4,5)
plot(N_f1,N_p1-max(N_p1),'displayName','fft:Mothersignal')

set(gca,'xscale','log')
legend('show')
subplot(2,4,6)
plot(N_f4,N_p4-max(N_p4),'displayName','Mothersignal+random phases')


set(gca,'xscale','log')
legend('show')
title('NoiseDM')
subplot(2,4,7)
plot(N_f2,N_p2-max(N_p2),'displayName','Current correct Tsim')
set(gca,'xscale','log')
legend('show')
title('NoiseDM')
subplot(2,4,8)
plot(N_f3,N_p3-max(N_p3),'displayName','filtered')
set(gca,'xscale','log')
legend('show')
title('NoiseDM')
hold off
ylabel('power/frequency  (DB/FS)')
%% self created PSD (see ArtificialPSDtest script) 
%here the artificial PSD was finetuned and validated + effect of fiulter
%was checked. but can be omitted now
Tsim=t2;
N = length(Tsim);
Fs = 1/(Tsim(2)-Tsim(1));
dF = Fs/N;
posf = 0:dF:Fs/2; 



wn = 14; 
theta = 0.7*rand()+0.3; 
w2 = 10^3;
alpha1 = rand()+2
alpha2 = 2*rand()+2
noise2 =@(x) 10.^(randn(1,length(x)));
noise3 = @(x) lognrnd(0,0.5,1,length(x))
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
fun = @(x) y(x).*noise3(x);
L = N;

if mod(N,2)==0
    abs_freqs=[posf,fliplr(posf(2:end-1))]; 
    PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
        fun(posf(end)),fun(fliplr(posf(2:end-1)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end-1) = 2*PSDpos(2:end-1); 
    
else
    abs_freqs=[posf,fliplr(posf(2:end))]; 
    PSD = [fun(posf(1)),fun(posf(2:end))/2,...
        fun(fliplr(posf(2:end)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end) = 2*PSDpos(2:end);
end
s0dft = sqrt(N*PSD);

filterdt = 100*1e-6;
windowSize = filterdt*Fs;
b = (1/windowSize)*ones(1,int16(windowSize));
Phifk = 2*pi*rand(1,L);   %random phases
Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
Iarray_temp = ifft(Zfk,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
Iarray_temp = Iarray_temp./max(abs(Iarray_temp));
[N_p4,N_f4] = genPSD(Iarray_temp,Tsim,Fs,'db');
Iarray = interp1(Tsim,Iarray_temp,Tsim);
[N_p2,N_f2] = genPSD(Iarray,Tsim,Fs,'db');
Iarray2 = filter(b,1,Iarray);
Iarray2 = Iarray2./max(abs(Iarray2));
[N_p3,N_f3] = genPSD(Iarray2,Tsim,Fs,'db');

figure
subplot(2,4,[1,4])
plot(Tsim,Iarray_temp,'displayName','Iarray_temp')
hold on
plot(Tsim,Iarray,'displayName','Iarray')
plot(Tsim,Iarray2,'displayName','Iarray2')
hold off
legend('show')
subplot(2,4,5)
plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')

set(gca,'xscale','log')
legend('show')
subplot(2,4,6)
plot(N_f4,N_p4,'displayName','Mothersignal+random phases')


set(gca,'xscale','log')
legend('show')
title('NoiseDM')
subplot(2,4,7)
plot(N_f2,N_p2,'displayName','Current correct Tsim')
set(gca,'xscale','log')
legend('show')
title('NoiseDM')
subplot(2,4,8)
plot(N_f3,N_p3,'displayName','filtered')
set(gca,'xscale','log')
legend('show')
title('NoiseDM')
hold off
ylabel('power/frequency  (DB/FS)')


%% final test on how it is in getIarray
clear all
close all

Fs1 = 40e3;
Fs2 = 40e6;


t = 0:1/Fs1:1;
t2 = 0:1/Fs2:0.1;



Tsim=t2;
N = 2 %number of currenttraces required
L = length(Tsim);
Fs = 1/(Tsim(2)-Tsim(1));
dF = Fs/L;
posf = 0:dF:Fs/2; 



wn = 14; 
theta = 0.7*rand()+0.3; 
w2 = 10^3;
alpha1 = rand()+2
alpha2 = 2*rand()+2
% noise2 =@(x) 10.^(randn(1,length(x)));
noise3 = @(x) lognrnd(0,0.5,1,length(x));
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
fun = @(x) y(x).*noise3(x);


if mod(L,2)==0
    abs_freqs=[posf,fliplr(posf(2:end-1))]; 
    PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
        fun(posf(end)),fun(fliplr(posf(2:end-1)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end-1) = 2*PSDpos(2:end-1); 
    
else
    abs_freqs=[posf,fliplr(posf(2:end))]; 
    PSD = [fun(posf(1)),fun(posf(2:end))/2,...
        fun(fliplr(posf(2:end)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end) = 2*PSDpos(2:end);
end
s0dft = sqrt(L*PSD);


Phifk = 2*pi*rand(N,L);   %random phases
Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained

Iarray_temp = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)

Iarray_temp = Iarray_temp./max(abs(Iarray_temp),[],2);


for iN=1:N
    
[N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray_temp(iN,:),Tsim,Fs,'db');
end


iNs = randperm(N,min(10,N));
figure
subplot(2,2,[1,2])
for iN=iNs
plot(Tsim,Iarray_temp(iN,:),'displayName','Iarray_temp')
hold on
end
legend('show')
ylim([-1.5,1.5])
subplot(2,2,3)
plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')

set(gca,'xscale','log')
legend('show')
subplot(2,2,4)
for iN=iNs
plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
hold on
end
hold off


set(gca,'xscale','log')
legend('show')
title('NoiseDM')


%% final test on how it is in getIarray Best case
% change of artificial PSD to be zero for high f
clear all
% close all

Fs1 = 40e3;
Fs2 = 40e6;


t = 0:1/Fs1:1;
t2 = 0:1/Fs2:0.1;



Tsim=t2;
N = 2 %number of currenttraces required
L = length(Tsim);
Fs = 1/(Tsim(2)-Tsim(1));
dF = Fs/L;
posf = 0:dF:Fs/2; 



wn = 14; 
theta = 0.7*rand()+0.3; 
w2 = 10^3;
alpha1 = rand()+2
alpha2 = 2*rand()+2
alpha3 = 10;
w3 = 10^4;
% noise2 =@(x) 10.^(randn(1,length(x)));
noise3 = @(x) lognrnd(0,0.5,1,length(x));
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2).*(1+(x/w3).^alpha3));
fun = @(x) y(x).*noise3(x).*double(x<10^4);


if mod(L,2)==0
    abs_freqs=[posf,fliplr(posf(2:end-1))]; 
    PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
        fun(posf(end)),fun(fliplr(posf(2:end-1)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end-1) = 2*PSDpos(2:end-1); 
    
else
    abs_freqs=[posf,fliplr(posf(2:end))]; 
    PSD = [fun(posf(1)),fun(posf(2:end))/2,...
        fun(fliplr(posf(2:end)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end) = 2*PSDpos(2:end);
end
s0dft = sqrt(L*PSD);


Phifk = 2*pi*rand(N,L);   %random phases
Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained

Iarray_temp = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)

Iarray_temp = Iarray_temp./max(abs(Iarray_temp),[],2);


for iN=1:N
    
[N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray_temp(iN,:),Tsim,Fs,'db');
end


iNs = randperm(N,min(10,N));
figure
subplot(2,2,[1,2])
for iN=iNs
plot(Tsim,Iarray_temp(iN,:),'displayName','Iarray_temp')
hold on
end
legend('show')
ylim([-1.5,1.5])
subplot(2,2,3)
plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')

set(gca,'xscale','log')
legend('show')
subplot(2,2,4)
for iN=iNs
plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
hold on
end
hold off


set(gca,'xscale','log')
legend('show')
title('NoiseDM')

%% final test on how it is in getIarray worst case
% change of artificial PSD to be only first powerlaw followed for high f
clear all
close all

Fs1 = 40e3;
Fs2 = 40e6;


t = 0:1/Fs1:1;
t2 = 0:1/Fs2:0.1;



Tsim=t2;
N = 2 %number of currenttraces required
L = length(Tsim);
Fs = 1/(Tsim(2)-Tsim(1));
dF = Fs/L;
posf = 0:dF:Fs/2;



wn = 14; 
theta = 0.7*rand()+0.3; 

alpha1 = rand()+2


% noise2 =@(x) 10.^(randn(1,length(x)));
noise3 = @(x) lognrnd(0,0.5,1,length(x));
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1));
fun = @(x) y(x).*noise3(x);


if mod(L,2)==0
    abs_freqs=[posf,fliplr(posf(2:end-1))]; 
    PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
        fun(posf(end)),fun(fliplr(posf(2:end-1)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end-1) = 2*PSDpos(2:end-1); 
    
else
    abs_freqs=[posf,fliplr(posf(2:end))]; 
    PSD = [fun(posf(1)),fun(posf(2:end))/2,...
        fun(fliplr(posf(2:end)))/2];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end) = 2*PSDpos(2:end);
end
s0dft = sqrt(L*PSD);


Phifk = 2*pi*rand(N,L);   %random phases
Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained

Iarray_temp = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)

Iarray_temp = Iarray_temp./max(abs(Iarray_temp),[],2);


for iN=1:N
    
[N_p4(iN,:),N_f4(iN,:)] = genPSD(Iarray_temp(iN,:),Tsim,Fs,'db');
end


iNs = randperm(N,min(10,N));
figure
subplot(2,2,[1,2])
for iN=iNs
plot(Tsim,Iarray_temp(iN,:),'displayName','Iarray_temp')
hold on
end
legend('show')
ylim([-1.5,1.5])
subplot(2,2,3)
plot(posf,10*log10(PSDpos),'displayName','fft:Mothersignal')

set(gca,'xscale','log')
legend('show')
subplot(2,2,4)
for iN=iNs
plot(N_f4(iN,:),N_p4(iN,:),'displayName','Mothersignal+random phases')
hold on
end
hold off


set(gca,'xscale','log')
legend('show')
title('NoiseDM')

