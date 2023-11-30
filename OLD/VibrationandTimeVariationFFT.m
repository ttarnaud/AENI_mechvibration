% Run simulations for moving dipole Script outdated see function
% investBiologicalNoise
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
if exist('Fignr','var')
    clearvars -except 'Fignr'
    Fignr = Fignr+1;
else
    clear all
    Fignr = 10;
end
CM = 'jet';
GENGIF = false;
close all
plotIarray_flag = 1;
% Electrical properties
Npercolumn = 1e4; 
I = 1*1e-3*Npercolumn; %µA
sigma = 0.33;    %S/m
d = 500*10^-6; %distance between current source and current sink
dI = I*d;
stdI = 1;
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+1/10*double(t(1)~=t0).*gaussmf(t,[0.0005,t0+0.0005]);
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0;
Tau = 0.0005;


% US options
fus = 10e5; %Hz
wus = 2*pi*fus; %rad/s
Aus = 100e-9; %m
Phaseus = 0; %rad
Dirus = [0,0,1]; Dirus = Dirus/norm(Dirus);
USwave = @(t) Dirus.*Aus.*sin(wus.*t+Phaseus);

% Simulation settings
resUS = 4;
dt = (resUS*fus)^-1;
Tend =50*Tau;
Tsim = 0:dt:Tend;

% Geometric  properties
SolutionType = '3SphereS8.2R25';
[Options,RSphere,RPOI] = getSettings(SolutionType,1);

% POIs on a sphere if required currently not used

Theta = [0:pi/100:pi];
dPhi = pi/100;
[phi,theta]=meshgrid(0:dPhi:2*pi,Theta);
phisc = phi(:); thetasc = theta(:);
POIsSphere = RSphere*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];




TSTART = tic;
%% dipole in center and yz-direction
OrienDipole = [0,0,1];
POIs = [0,0,1;0,1,0;0,1/sqrt(2),1/sqrt(2)];
POIs = RSphere*POIs./vecnorm(POIs,2,2);
VRTOT = zeros(size(POIs,1),length(Tsim));
Totaldps = 1e2;
dps_run = 1e2;
irun_end = Totaldps/dps_run;
% Select Sources and sinks + select vibrator
[CSource,CSink] = GendpPos('fibonaccisingle','radial',d,RSphere-0.005,Totaldps-1);
CSource = vertcat(CSource,(0.065+d/2)*OrienDipole);
CSink = vertcat(CSink,(0.065-d/2)*OrienDipole);

for irun = 1:irun_end
% Gen I vals over time
tic
Iarray = getIarray('randomi','normal',Tsim,dps_run,I,d,stdI,[]);
toc
if irun==irun_end
IOI = I.*Alphafun(Tsim,Tau,0);
%IOI = I.*Itime(Tsim,d);
idxOscDip = dps_run;
Iarray(idxOscDip,:) = IOI;
else 
    idxOscDip = dps_run+1;
end
idxdps = (irun-1)*dps_run+1:irun*dps_run;
if plotIarray_flag
    figure
    PlotIarray(Iarray,Tsim)
end
disp('start simulation')
Settings = horzcat(Options,{'POI',POIs,'ShowSphere','wpoi','resUS',resUS});
VR=SimDipoleOsc(Tend,CSource(idxdps,:),CSink(idxdps,:),Iarray,USwave,1/fus,idxOscDip,Settings);
VRTOT = VRTOT+VR;
end
toc(TSTART)
VR = VRTOT;
%%
% function plotinfo(Tsim,Dirus,Aus,wus,Phaseus,POIs,VR,Alphafun,Tau)
%VR = VROLD+rednoise(size(VR,1),size(VR,2));
Display_flag = 1;
inputSignal = [Tsim;IOI];
if Display_flag
figure

POIs = [1,2,3];
subplot(4,3,2);
plot([Tsim(1:100)*1e6.*ones(3,1)]',[Dirus'.*Aus.*sin(wus.*Tsim(1:100)+Phaseus)]'*1e9)
legend({'x component','y component','z component'})
ylabel('displacement [nm]')
xlabel('Time [µs]')

subplot(4,3,3);
plot(Tsim*1e3,IOI)
ylabel('Amplitude dipole [µA]')
xlabel('Time [ms]')


for ifig =1:3
    ax{ifig} = subplot(4,3,3*ifig+1);
    nr = POIs(ifig);
    plot(Tsim*1e3,VR(nr,:))
    if ifig == 3
        xlabel('Time [ms]')
    end
    if ifig == 1
        title('Measured potential at POI [µV]')
    end
    ylabel({['\bf POI = ',num2str(nr)],''})
    
    subplot(4,3,3*ifig+2)
    Fs = 1/(Tsim(2)-Tsim(1));
    L = length(Tsim);
    %N = 2^nextpow2(L); % to get higher resolution willl giv sinc function
    N = L;
    Y = fft(VR(nr,:),N);
    freqs = 0:Fs/N:Fs/2-Fs/N;
    plot(freqs,abs(Y(1:length(freqs))))
    %plot(freqs,P1);
    xlim([0,1e3])
    if ifig == 1
        title('Single-Sided Amplitude Spectrum')
    end
    if ifig == 3
        xlabel('f [Hz]')
    end
    ylabel('|Fx|')
    %set(gca,'xscale','log','yscale','log')
    subplot(4,3,3*ifig+3)
    
    plot(freqs,abs(Y(1:length(freqs))))
    xlim([fus-1000,fus+1000])
    ylabel('|Fx|')
    set(gca,'XTickLabel',arrayfun(@num2str,get(gca,'xTick')*1e-6,'UniformOutput',false))
    if ifig == 1
        title('Amplitude Spectrum near 1 MHz')
    end
    if ifig == 3
        xlabel('f [MHz]')
    end
end
end
%%
DOIspos = (CSource(DOIindices,:)+CSink(DOIindices,:))/2;
dDOIPOI0 = vecnorm(DOIspos-POIs,2,2);
dDOIPOI1 = vecnorm(DOIspos+USwave(1/10*1/fus)-POIs,2,2);
mirror = double(dDOIPOI1>dDOIPOI0);

calcFourierGetsignal(VR,Tsim,[1,2,3],fus,2*1e3,inputSignal,mirror,Display_flag)
toc(TSTART)
%% functions
function PlotIarray(Iarray,Tsim)
subplot(2,1,1)
for i=1:10
    idx = randi(size(Iarray,1));
    plot(Tsim,Iarray(idx,:))
    hold on
end
hold off
ylabel('Current in dipole [µA]')
title('Selection of Currents in dipoles')
subplot(2,1,2)
plot(Tsim*1e3,Iarray(1,:))
hold on
plot(Tsim*1e3,Iarray(end,:))
hold off
xlabel('Time [ms]')
ylabel('Current in dipole [µA]')
end