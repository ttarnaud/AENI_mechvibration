% code copied from getIarray in Functions/subfunctions
close all
clear all
% the currents are derived from a aritificial PSD  based on
% observed trends in literature for low frequency(<100Hz) and model
% ouputs (Neuron) High frequency. A plausible sprectrum is created.
% This way minimizes interpolation errors. see An_spec script for
% more details. the powerlaws are fixed
rng(1)
fus = 1e6;
dt = 1/(20*fus);
Tend = 0.025;
Tsim = 0:dt:Tend;
N=5;
gv=0.3;
mcolors = lines(7);
I = 10;
sigmaI = 10;


L = length(Tsim);
Fs = 1/(Tsim(2)-Tsim(1));
dF = Fs/L;
posf = 0:dF:Fs/2;

wn = 14;
theta = 0.6;         %max is reached between wn and theta*wn
w2 = 10^3;                      %breakpoint higher powerlaw
alpha1 = 2;
alpha2 = 2;
%noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
%were tested see script An_spec
noise3 = @(x) lognrnd(0,0.5,1,length(x));
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
fun = @(x) y(x).*noise3(x); %this function determines artificial PSD right hand sided

%store parameters
extra_opt.alpha1 = alpha1; extra_opt.alpha2 = alpha2; extra_opt.theta = theta;
if mod(L,2)==0
    
    PSD = [fun(posf(1)),fun(posf(2:end-1))/2,...
        fun(posf(end)),fun(fliplr(posf(2:end-1)))/2]; %same way as components would be stored by using fft
    
    abs_freqs=[posf,fliplr(posf(2:end-1))];
    PSDpos = PSD(1:length(posf));
    PSDpos(2:end-1) = 2*PSDpos(2:end-1);

    
else
    
    PSD = [fun(posf(1)),fun(posf(2:end))/2,...
    fun(fliplr(posf(2:end)))/2];

    PSDpos = PSD(1:length(posf));
    PSDpos(2:end) = 2*PSDpos(2:end);
    abs_freqs=[posf,fliplr(posf(2:end))];
end

s0dft = sqrt(L*PSD); %PSD to magnitude of rotor

Phifk = 2*pi*rand(N,L);   %random phases
Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
Iarray = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
Iarray = Iarray./max(abs(Iarray),[],2);
Imax = normrnd(I,sigmaI,size(Iarray,1),1);

Iarray=Imax.*Iarray;

%% Creating figure
figure

subplot(2,1,1)
plot(posf,10*log10(PSDpos),'color',[gv,gv,gv],'linewidth',1)
hold on
x = logspace(0,7,1000);
plot(x,10*log10(y(x)),'linewidth',2,'color',mcolors(5,:))
xlabel('frequency [Hz]')
ylabel('PSD [dB]')
set(gca,{'xscale','box'},{'log','off'})


subplot(2,1,2)
for iN=1:N
    plot(Tsim*1000,Iarray(iN,:),'HandleVisibility','off','color',[gv,gv,gv],'linewidth',1)
    hold on
end

Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI
Ifun2 = alphaTrainFun(Tsim,'aTInput_220308.mat');

plot(Tsim*1000,I*Ifun(Tsim),'linewidth',2,'color',mcolors(2,:),'DisplayName','alpha-fun')
plot(Tsim*1000,I*Ifun2,'linewidth',2,'color',min(mcolors(2,:)+0.2,1),'DisplayName','alpha-train')
plot(nan,nan,'color',[gv,gv,gv],'linewidth',2,'DisplayName','background dipoles')
%legend('show','box','off')
%ylim([-1.1,1.1])
set(gca,'box','off')
xlabel('time [ms]')
ylabel('I(t) [\muA]') 
xlim([0,Tend*1000])




set(gcf,{'units','color','position','paperunits','papersize'},...
    {'centimeters',[1,1,1],[1,3,7.8,8.5],'centimeters',[6.8,5.5]})
set(findall(gcf,'type','axes'),'fontsize',9)

%% UltrasonicField
figure(10)
Afun =@(x,k) exp(-abs(x)/k);
x = -70:0.1:70;
for k=[5,20]
    plot(x,Afun(x,k),'color',min(mcolors(5,:)+0.2*double(k==20),1),'linewidth',2,'alpha',0.5)
    hold on
end
hold off
xlabel('|{\bf r} - {\bf r}_{dp,DOI} |')
ylabel('displacement profile')
set(gca,'box','off')
ax = gca;
ax.YAxis.Visible = 'off'; 

set(gcf,{'units','color','position','paperunits','papersize'},...
    {'centimeters',[1,1,1],[1,3,1+2*6.8,3],'centimeters',[1+2*6.8,3]})
set(findall(gcf,'type','axes'),'fontsize',9)


    
