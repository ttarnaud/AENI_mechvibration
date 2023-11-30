close all
clc
dc = 0.04;
prf = 2000;
fus = 1e6;
Fs = fus*10;
dt = 1/(Fs);
Tsim = 0.05;
t = 0:dt:Tsim-dt;
L = length(t);
f = Fs*(0:(L/2))/L;

% target functions
afun = @(t,tau) t/tau.*exp(1-t/tau);
squarefun = t<=dc/prf;
squarewave = mod(t,1/prf)<=dc/prf;
uswave = sin(2*pi*fus*t);


% construct bio noise
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

if mod(L,2)==0

    PSD = [fun(f(1)),fun(f(2:end-1))/2,...
        fun(f(end)),fun(fliplr(f(2:end-1)))/2]; %same way as components would be stored by using fft
%     if debug_flag
%         abs_freqs=[f,fliplr(f(2:end-1))];
%         PSDpos = PSD(1:length(f));
%         PSDpos(2:end-1) = 2*PSDpos(2:end-1);
%     end

else

    PSD = [fun(f(1)),fun(f(2:end))/2,...
        fun(fliplr(f(2:end)))/2];
%     if debug_flag
%         PSDpos = PSD(1:length(f));
%         PSDpos(2:end) = 2*PSDpos(2:end);
%         abs_freqs=[f,fliplr(f(2:end))];
%     end
end

s0dft = sqrt(L*PSD); %PSD to magnitude of rotor

Phifk = 2*pi*rand(1,L);   %random phases
Zfk = s0dft.*exp(complex(0,Phifk)); %add random phases while same frequency content maintained
Ibionoise = ifft(Zfk,[],2,'symmetric'); %symmetric to get real output however also some artefacts but not important here (see calcgetfourier signal)
Ibionoise = Ibionoise./max(abs(Ibionoise),[],2);




%% plot time domain
subplot(2,2,1)
plot(t,uswave)
hold on
plot(t,squarewave);
hold off
subplot(2,2,3)
plot(t,squarewave.*uswave)
xlabel('time')

subplot(2,2,2)
plot(t,afun(t,0.005))
yyaxis right
plot(t,afun(t,0.005).*uswave.*squarewave);
hold off
subplot(2,2,4)
plot(t,Ibionoise)
xlabel('time')

%% fft
SW = fft(squarewave);
UW = fft(uswave);
SF = fft(squarefun);
AF = fft(afun(t,0.005));
SWUW = fft(squarewave.*uswave);
SFUW = fft(squarefun.*uswave);
AFSWUW = fft(afun(t,0.005).*squarewave.*uswave);
AFSWUWbionoise = fft(afun(t,0.005).*squarewave.*uswave+1e3*Ibionoise);


SWP1 = singlefft(SW,L);
UWP1 = singlefft(UW,L);
SFP1 = singlefft(SF,L);
SWUWP1 = SWUWP1fun(t,dc,prf,fus,L);
SFUWP1 = SFUWP1fun(t,dc,prf,fus,L);
AFSWUWP1 = AFSWUWP1fun(t,dc,prf,fus,L,afun(t,0.005));
AFSWUWPbn1 = AWSWUWPbn1fun(t,dc,prf,fus,L,1e3*Ibionoise,afun(t,0.005));


idx = abs(f-fus)<1e-6;
fprintf('signal at fus:\n')
fprintf('  squarewave: %5.2e\n',SWP1(idx))
fprintf('  ultrsoundwave: %5.2e\n',UWP1(idx))
fprintf('  squarefun: %5.2e\n',SFP1(idx))
fprintf('  squarewave * ultrasoundwave: %5.2e\n',SWUWP1(idx))
fprintf('  squarefun * ultrasoundwave: %5.2e\n',SFUWP1(idx))
fprintf('  alphafun * squarewave * ultrasoundwave: %5.2e\n',AFSWUWP1(idx))
fprintf('  alphafun* squarewave * ultrasoundwave + bionoise: %5.2e\n',AFSWUWPbn1(idx))



figure()
subplot(3,2,1)
plot(f,UWP1)
ylabel('uswave')

subplot(3,2,3)
plot(f,SWP1)
ylabel('square wave')
yyaxis right
plot(f,SFP1)
ylabel('square fun')

subplot(3,2,5)
plot(f,UWP1)
ylabel('uswave')
hold on
plot(f,SWUWP1)
ylabel('sq wave . uswave')
yyaxis right
plot(f,SFUWP1)
ylabel('sq fun . uswave')
xlabel('frequency')

subplot(3,2,2)
plot(f,AFSWUWP1)
ylabel('afun * sqwave * uswave')
subplot(3,2,4)
plot(f,AFSWUWPbn1)
ylabel('plus noise')
xlabel('frequency')

%% Normalized

nSWP1 = SWP1/max(SWP1);
nUWP1 = UWP1/max(UWP1);
nSFP1 = SFP1/max(SFP1);
nSWUWP1 = SWUWP1/max(SWUWP1);
nSFUWP1 = SFUWP1/max(SFUWP1);
nAFSWUWP1 = AFSWUWP1/max(AFSWUWP1);
nAFSWUWPbn1 = AFSWUWPbn1/max(AFSWUWPbn1);

fprintf('normalized signal at fus:\n')
fprintf('  squarewave: %5.2e\n',nSWP1(idx))
fprintf('  ultrsoundwave: %5.2e\n',nUWP1(idx))
fprintf('  squarefun: %5.2e\n',nSFP1(idx))
fprintf('  squarewave * ultrasoundwave: %5.2e\n',nSWUWP1(idx))
fprintf('  squarefun * ultrasoundwave: %5.2e\n',nSFUWP1(idx))
fprintf('  alphafun * squarewave * ultrasoundwave: %5.2e\n',nAFSWUWP1(idx))
fprintf('  alphafun* squarewave * ultrasoundwave + bionoise: %5.2e\n',nAFSWUWPbn1(idx))


figure()
subplot(3,2,1)
plot(f,nUWP1)
ylabel('uswave')

subplot(3,2,3)
plot(f,nSWP1)
ylabel('square wave')
yyaxis right
plot(f,nSFP1)
ylabel('square fun')

subplot(3,2,5)
plot(f,nUWP1)
ylabel('uswave')
hold on
plot(f,nSWUWP1)
ylabel('sq wave . uswave')
yyaxis right
plot(f,nSFUWP1)
ylabel('sq fun . uswave')
xlabel('frequency')

subplot(3,2,2)
plot(f,nAFSWUWP1)
ylabel('afun * sqwave * uswave')
subplot(3,2,4)
plot(f,nAFSWUWPbn1)
ylabel('plus noise')
xlabel('frequency')

% centered at fus +fband
fband = 10e4;
xlims = [fus-fband,fus+fband];
figure()
subplot(3,2,1)
plot(f,nUWP1)
ylabel('uswave')
xlim(xlims)

subplot(3,2,3)
plot(f,nSWP1)
ylabel('square wave')
yyaxis right
plot(f,nSFP1)
ylabel('square fun')
xlim(xlims)

subplot(3,2,5)
plot(f,nUWP1)
ylabel('uswave')
hold on
plot(f,nSWUWP1)
ylabel('sq wave . uswave')
yyaxis right
plot(f,nSFUWP1)
ylabel('sq fun . uswave')
xlabel('frequency')
xlim(xlims)

subplot(3,2,2)
plot(f,nAFSWUWP1)
ylabel('afun * sqwave * uswave')
xlim(xlims)
subplot(3,2,4)
plot(f,nAFSWUWPbn1)
ylabel('plus noise')
xlabel('frequency')
xlim(xlims)


%% effect dc

dcs = [0.04,0.1,1];
% centered at fus +fband
fband = 10e4;
xlims = [fus-fband,fus+fband];
figure()
for idc=1:length(dcs)

    subplot(3,3,0+idc)
    plot(f,SWUWP1fun(t,dcs(idc),prf,fus,L))
    ylabel('sq wave . uswave')
    yyaxis right
    plot(f,SFUWP1fun(t,dcs(idc),prf,fus,L))
    ylabel('sq fun . uswave')
    xlabel('frequency')
    xlim(xlims)

    subplot(3,3,3+idc)
    plot(f,AFSWUWP1fun(t,dcs(idc),prf,fus,L,afun(t,0.005)))
    ylabel('afun * sqwave * uswave')
    xlim(xlims)
    
    subplot(3,3,6+idc)
    plot(f,AWSWUWPbn1fun(t,dcs(idc),prf,fus,L,1e3*Ibionoise,afun(t,0.005)))
    ylabel('1e3 plus noise')
    xlabel('frequency')
    xlim(xlims)
end


% centered at fus +fband
fband = 2e3;
xlims = [fus-fband,fus+fband];
figure()
for idc=1:length(dcs)

    subplot(3,3,0+idc)
    plot(f,SWUWP1fun(t,dcs(idc),prf,fus,L))
    ylabel('sq wave . uswave')
    yyaxis right
    plot(f,SFUWP1fun(t,dcs(idc),prf,fus,L))
    ylabel('sq fun . uswave')
    xlabel('frequency')
    xlim(xlims)

    subplot(3,3,3+idc)
    plot(f,AFSWUWP1fun(t,dcs(idc),prf,fus,L,afun(t,0.005)))
    ylabel('afun * sqwave * uswave')
    xlim(xlims)
    
    subplot(3,3,6+idc)
    plot(f,AWSWUWPbn1fun(t,dcs(idc),prf,fus,L,0.5e4*Ibionoise,afun(t,0.005)))
    ylabel(' 0.5e4 plus noise')
    xlabel('frequency')
    xlim(xlims)
end

% centered at prf +fband
fband = 2e3;
xlims = [prf-fband,prf+fband];
figure()
for idc=1:length(dcs)

    subplot(3,3,0+idc)
    plot(f,SWUWP1fun(t,dcs(idc),prf,fus,L))
    ylabel('sq wave . uswave')
    yyaxis right
    plot(f,SFUWP1fun(t,dcs(idc),prf,fus,L))
    ylabel('sq fun . uswave')
    xlabel('frequency')
    xlim(xlims)

    subplot(3,3,3+idc)
    plot(f,AFSWUWP1fun(t,dcs(idc),prf,fus,L,afun(t,0.005)))
    ylabel('afun * sqwave * uswave')
    xlim(xlims)
    
    subplot(3,3,6+idc)
    plot(f,AWSWUWPbn1fun(t,dcs(idc),prf,fus,L,0.5e4*Ibionoise,afun(t,0.005)))
    ylabel(' 1e4 plus noise')
    xlabel('frequency')
    xlim(xlims)
end

%% effect prf

prfs = [500, 2e3, 1e4];
dc = 0.04;
% centered at fus +fband
fband = 10e4;
xlims = [fus-fband,fus+fband];
figure()
for iprf=1:length(prfs)

    subplot(3,3,0+iprf)
    plot(f,SWUWP1fun(t,dc,prfs(iprf),fus,L))
    ylabel('sq wave . uswave')
    yyaxis right
    plot(f,SFUWP1fun(t,dc,prfs(iprf),fus,L))
    ylabel('sq fun . uswave')
    xlabel('frequency')
    xlim(xlims)

    subplot(3,3,3+iprf)
    plot(f,AFSWUWP1fun(t,dc,prfs(iprf),fus,L,afun(t,0.005)))
    ylabel('afun * sqwave * uswave')
    xlim(xlims)
    
    subplot(3,3,6+iprf)
    plot(f,AWSWUWPbn1fun(t,dc,prfs(iprf),fus,L,1e3*Ibionoise,afun(t,0.005)))
    ylabel('1e3 plus noise')
    xlabel('frequency')
    xlim(xlims)
end


% centered at fus +fband
fband = 2e3;
xlims = [fus-fband,fus+fband];
figure()
for iprf=1:length(prfs)


    subplot(3,3,0+iprf)
    plot(f,SWUWP1fun(t,dc,prfs(iprf),fus,L))
    ylabel('sq wave . uswave')
    yyaxis right
    plot(f,SFUWP1fun(t,dc,prfs(iprf),fus,L))
    ylabel('sq fun . uswave')
    xlabel('frequency')
    xlim(xlims)

    subplot(3,3,3+iprf)
    plot(f,AFSWUWP1fun(t,dc,prfs(iprf),fus,L,afun(t,0.005)))
    ylabel('afun * sqwave * uswave')
    xlim(xlims)
    
    subplot(3,3,6+iprf)
    plot(f,AWSWUWPbn1fun(t,dc,prfs(iprf),fus,L,0.5e4*Ibionoise,afun(t,0.005)))
    ylabel(' 0.5e4 plus noise')
    xlabel('frequency')
    xlim(xlims)
end

% centered at prf +fband
fband = 2e3;
xlims = [prf-fband,prf+fband];
figure()
for iprf=1:length(prfs)
    xlims = [prfs(iprf)-fband,prfs(iprf)+fband];
    subplot(3,3,0+iprf)
    plot(f,SWUWP1fun(t,dc,prfs(iprf),fus,L))
    ylabel('sq wave . uswave')
    yyaxis right
    plot(f,SFUWP1fun(t,dc,prfs(iprf),fus,L))
    ylabel('sq fun . uswave')
    xlabel('frequency')
    xlim(xlims)

    subplot(3,3,3+iprf)
    plot(f,AFSWUWP1fun(t,dc,prfs(iprf),fus,L,afun(t,0.005)))
    ylabel('afun * sqwave * uswave')
    xlim(xlims)
    
    subplot(3,3,6+iprf)
    plot(f,AWSWUWPbn1fun(t,dc,prfs(iprf),fus,L,0.5e4*Ibionoise,afun(t,0.005)))
    ylabel(' 1e4 plus noise')
    xlabel('frequency')
    xlim(xlims)
end
%%
target = Pintm(abs(f-prf)<1e-6)/max(Pintm);

myfun = @(A) target*2*pi*(prf-fus)./A-sin((2*pi*(prf-fus))./A);
myfun2 = @(A) target-sinc((2*pi*(prf-fus))./A);
fsolve(myfun,prf/0.04*2)
fsolve(myfun2,prf/0.04*2)

myfun = @(A) target*(prf-fus)./A-sin(((prf-fus))./A);
myfun2 = @(A) target-sinc(((prf-fus))./A);
fsolve(myfun,prf/0.04*2)
fsolve(myfun2,prf/0.04*2)

function P1 = singlefft(Y,L)
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
end

function SWP1 = SWP1fun(t,dc,prf,L)
squarewave = mod(t,1/prf)<=dc/prf;
SWUW = fft(squarewave);
SWP1 = singlefft(SWUW,L);
end

function SFP1 = SFP1fun(t,dc,prf,L)
squarefun = t<=dc/prf;
SFUW = fft(squarefun);
SFP1 = singlefft(SFUW,L);
end

function SWUWP1 = SWUWP1fun(t,dc,prf,fus,L)
squarewave = mod(t,1/prf)<=dc/prf;
uswave = sin(2*pi*fus*t);
SWUW = fft(squarewave.*uswave);
SWUWP1 = singlefft(SWUW,L);
end

function SFUWP1 = SFUWP1fun(t,dc,prf,fus,L)
squarefun = t<=dc/prf;
uswave = sin(2*pi*fus*t);
SFUW = fft(squarefun.*uswave);
SFUWP1 = singlefft(SFUW,L);
end

function AFSWUWP1 = AFSWUWP1fun(t,dc,prf,fus,L,afun)
squarewave = mod(t,1/prf)<=dc/prf;
uswave = sin(2*pi*fus*t);
Y = fft(afun.*squarewave.*uswave);
AFSWUWP1 = singlefft(Y,L);
end
function AWSWUWPbn1 = AWSWUWPbn1fun(t,dc,prf,fus,L,Ibionoise,afun)
squarewave = mod(t,1/prf)<=dc/prf;
uswave = sin(2*pi*fus*t);
Y = fft(afun.*squarewave.*uswave+Ibionoise);
AWSWUWPbn1 = singlefft(Y,L);
end