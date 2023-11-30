function [RMS,SNR,TSTOP] = investBiologicalNoise(varargin)
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
TSTART = tic;
%% Initialise default variables
PLOT = 0;
Display = 0;
plotIarray_flag = 0;
showSphere='';
ParallelCompute_flag = 1;
ThermalNoiseAmp = 0;
scale_flag = 1;
% Electrical properties
Npercolumn = 1e4; 
meanI = 1*1e-3*Npercolumn; %µA
sigma = 0.33;    %S/m
d = 500*10^-6; %distance between current source and current sink
dI = meanI*d;
stdI = 3;
%Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+1/10*double(t(1)~=t0).*gaussmf(t,[0.0005,t0+0.0005]);
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0;
Tau = 0.005;
AlphaDelay = 0;
Ifun = @(t) Alphafun(t,Tau,AlphaDelay);

% US options
fus = 10e5; %Hz
Aus = 10e-9; %m
Phaseus = 0; %rad
Dirus = [0,0,1];
fbandwidth = 1e3;

% Simulation settings
resUS = 10;
Tend = 0.025;


% Geometric  properties
SolutionType = '3SphereS8.2R25';

% POIs
POIs = [0,0,1;0,1,0;0,1/sqrt(2),1/sqrt(2)];

% Dipoles
OrienDipole = [0,0,1];
%posDp = [0,0,0.065];
Totaldps = 1e2;
dps_run = 1e2;
dpDistribution = 'fibonaccisingle';
dpOrientation = 'radial';
dpI_time = 'randominewgen';
dpI_space = 'normal';
%% change default variables
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'Plot'))
            PLOT = varargin{find(strcmpi(varargin,'Plot'))+1};
        end
        if any(strcmpi(varargin,'meanI'))
            meanI = varargin{find(strcmpi(varargin,'meanI'))+1};
        end
        if any(strcmpi(varargin,'stdI'))
            stdI = varargin{find(strcmpi(varargin,'stdI'))+1};
        end
        if any(strcmpi(varargin,'Ifun'))
            Ifun = varargin{find(strcmpi(varargin,'Ifun'))+1};
        end
        if any(strcmpi(varargin,'fus'))
            fus = varargin{find(strcmpi(varargin,'fus'))+1};
        end
        if any(strcmpi(varargin,'Aus'))
            Aus = varargin{find(strcmpi(varargin,'Aus'))+1};
        end
        if any(strcmpi(varargin,'Dirus'))
            Dirus = varargin{find(strcmpi(varargin,'Dirus'))+1};
        end
        if any(strcmpi(varargin,'resUS'))
            resUS = varargin{find(strcmpi(varargin,'resUS'))+1};
        end
        if any(strcmpi(varargin,'Tend'))
            Tend = varargin{find(strcmpi(varargin,'Tend'))+1};
        end
        if any(strcmpi(varargin,'SolutionType'))
            SolutionType = varargin{find(strcmpi(varargin,'SolutionType'))+1};
        end
        if any(strcmpi(varargin,'POIs'))
            POIs = varargin{find(strcmpi(varargin,'POIs'))+1};
        end
        if any(strcmpi(varargin,'Totaldps'))
            Totaldps = varargin{find(strcmpi(varargin,'Totaldps'))+1};
        end
        if any(strcmpi(varargin,'dps_run'))
            dps_run = varargin{find(strcmpi(varargin,'dps_run'))+1};
        end
        if any(strcmpi(varargin,'dpDistribution'))
            dpDistribution = varargin{find(strcmpi(varargin,'dpDistribution'))+1};
        end
        if any(strcmpi(varargin,'dpOrientation'))
            dpOrientation = varargin{find(strcmpi(varargin,'dpOrientation'))+1};
        end
        if any(strcmpi(varargin,'dpI_time'))
            dpI_time = varargin{find(strcmpi(varargin,'dpI_time'))+1};
        end
        if any(strcmpi(varargin,'dpI_space'))
            dpI_space = varargin{find(strcmpi(varargin,'dpI_space'))+1};
        end
        if any(strcmpi(varargin,'posDp'))
            posDp = varargin{find(strcmpi(varargin,'posDp'))+1};
        end
        if any(strcmpi(varargin,'OrienDipole'))
            OrienDipole = varargin{find(strcmpi(varargin,'OrienDipole'))+1};
        end
        if any(strcmpi(varargin,'ParallelCompute'))
            ParallelCompute_flag = varargin{find(strcmpi(varargin,'ParallelCompute'))+1};
        end
        if any(strcmpi(varargin,'ThermalNoiseAmp'))
            ThermalNoiseAmp = varargin{find(strcmpi(varargin,'ThermalNoiseAmp'))+1};
        end
        if any(strcmpi(varargin,'scale'))
            scale_flag = varargin{find(strcmpi(varargin,'scale'))+1};            
        end
    end
end

%% Update simulation settings
% US-wave
wus = 2*pi*fus; %rad/s
Dirus = Dirus/norm(Dirus);
USwave = @(t) Dirus.*Aus.*sin(wus.*t+Phaseus);
% Simulation settings
Tend = ceil(fus*Tend)/fus; % recalculate Tend such that frequency spectrum contains fus
dt = (resUS*fus)^-1;
Tsim = 0:dt:Tend;
% Geometric properties
[Options,RSphere] = getSettings(SolutionType,1);
% POIs
if scale_flag
POIs = RSphere*POIs./vecnorm(POIs,2,2);
end
% Plot settings
if PLOT
    plotIarray_flag = 1;
    showSphere = 'off';
    Display = 1;
end
% Select Sources and sinks + select vibrator
if strcmpi(SolutionType,'Mouse4Sphere~fair0')
dRSphere = 0.0005;
else
dRSphere = 0.005;
end
posDp = 0;
[CSource,CSink] = GendpPos(dpDistribution,dpOrientation,d,RSphere-dRSphere,Totaldps-1);
CSource = vertcat(CSource,(posDp+d/2)*OrienDipole);
CSink = vertcat(CSink,(posDp-d/2)*OrienDipole);
DOIindices = Totaldps;
TSTOP.startSim = toc(TSTART);

%% Start simulations
VR = zeros(size(POIs,1),length(Tsim));
VRosc = zeros(size(POIs,1),length(Tsim));
VRnoise = zeros(size(POIs,1),length(Tsim));
irun_end = Totaldps/dps_run;
if irun_end-round(irun_end,0)~=0
    error('irun_end must be integer')
end
if logical(mod(dps_run,4))
    warning('dps_run preferably multiple of 4')
end

for irun = 1:irun_end
% Generate I(t)
Iarray = getIarray(dpI_time,dpI_space,Tsim,dps_run,meanI,d,stdI,[]);
TSTOP.Itgen(irun) = toc(TSTART);
DOIflag = DOIindices>(irun-1)*dps_run & DOIindices<=irun*dps_run;
if any(DOIflag)
    IOI = meanI.*Ifun(Tsim);
    idxOscDip = DOIindices(DOIflag)-(irun-1)*dps_run;
    Iarray(idxOscDip,:) = IOI;
else 
    idxOscDip = dps_run+1;
end
idxdps = (irun-1)*dps_run+1:irun*dps_run;
if plotIarray_flag
    figure
    PlotIarray(Iarray,Tsim)
end
if Display
disp(['start simulation ',num2str(irun),' of ',num2str(irun_end)])
end
Settings = horzcat(Options,{'POI',POIs,'ShowSphere',showSphere,'resUS',resUS,'ParallelCompute',ParallelCompute_flag,'scale',scale_flag});
[VRrun,VRosc_run,VRnoise_run] = SimDipoleOsc(Tend,CSource(idxdps,:),CSink(idxdps,:),Iarray,USwave,1/fus,idxOscDip,Settings);
TSTOP.VRgen(irun) = toc(TSTART);
VR = VR+VRrun; VRosc = VRosc+VRosc_run; VRnoise = VRnoise+VRnoise_run;
end
ThermalNoise = ThermalNoiseAmp.*rand(size(VR));
VR = VR + ThermalNoise; VRosc = VRosc; VRnoise = VRnoise + ThermalNoise;
TSTOP.endsim = toc(TSTART);

%% plotFigure
% function plotinfo(Tsim,Dirus,Aus,wus,Phaseus,POIs,VR,Alphafun,Tau)
%VR = VROLD+rednoise(size(VR,1),size(VR,2));
POIsidx = 1:min(size(POIs,1),3);
inputSignal = [Tsim;IOI];
if PLOT
figure

subplot(4,3,2);
plot([Tsim(1:100)*1e6.*ones(3,1)]',[Dirus'.*Aus.*sin(wus.*Tsim(1:100)+Phaseus)]'*1e9)
legend({'x component','y component','z component'})
ylabel('displacement [nm]')
xlabel('Time [µs]')

subplot(4,3,3);
plot(Tsim*1e3,IOI)
ylabel('Amplitude dipole [µA]')
xlabel('Time [ms]')


for ifig =1:length(POIsidx)
    ax{ifig} = subplot(4,3,3*ifig+1);
    nr = POIsidx(ifig);
    plot(Tsim*1e3,VR(nr,:))
    if ifig == length(POIsidx)
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
%% Signal reconstruction and metric determination
% define if signal needs to be mirrored
DOIspos = (CSource(DOIindices,:)+CSink(DOIindices,:))/2;
dDOIPOI0 = vecnorm(DOIspos-POIs,2,2);
dDOIPOI1 = vecnorm(DOIspos+USwave(1/10*1/fus)-POIs,2,2);
mirror = double(dDOIPOI1>dDOIPOI0);

    
[reconsSignalVR,powerSignalVR,timeReconS] = calcFourierGetsignal(VR,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,PLOT);
[reconsSignalVRosc,powerSignalVRosc] = calcFourierGetsignal(VRosc,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,0);
[reconsSignalVRnoise,powerSignalVRnoise] = calcFourierGetsignal(VRnoise,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,0);
TSTOP.signalreconstruction = toc(TSTART);
RMS = struct;
for iPOI=1:length(reconsSignalVR)
    % SNR
    Powerosc = calcPower(reconsSignalVRosc{iPOI});
    Powernoise = calcPower(reconsSignalVRnoise{iPOI});
    SNR(1,iPOI) = 10*log10(Powerosc/Powernoise);
    SNR(2,iPOI) = 10*log10(powerSignalVRosc(iPOI)/powerSignalVRnoise(iPOI));
    %RMS
    [RMS(iPOI).RMSnorm,RMS(iPOI).RMSosc] = ...
        calcRMS(reconsSignalVR{iPOI},reconsSignalVRosc{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT);   

end

%% functions
    function PlotIarray(Iarray,Tsim)
    subplot(2,1,1)
    for ipArray=1:10
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
    function Power = calcPower(Signal)
        s0dft = fft(Signal);
        psds0 = (1/(length(Signal))) * abs(s0dft).^2;
        Power = sum(psds0);
    end
    function [RMS,RMSosc,RMSabs] = calcRMS(Signal,Signalosc,InputSignal,timeReconS,mirror,disp)
        InS = interp1(InputSignal(1,:),InputSignal(2,:),timeReconS);
        % Check if signal lengths are the same
        if length(Signal)~=length(Signalosc) || length(Signal)~=length(InS)
            error('signal lengths not the same')
        end
        Signal = (-1)^mirror.*Signal;
        Signalosc = (-1)^mirror.*Signalosc;
        
        % RMS
            RMS = sqrt(mean((Signal/max(abs(Signal))-InS/max(abs(InS))).^2));
        %RMSabs
            RMSabs = sqrt(mean((abs(Signal)/max(abs(Signal))-abs(InS)/max(abs(InS))).^2));
        %RMSosc
            RMSosc = sqrt(mean((Signal/max(abs(Signal))-Signalosc/max(abs(Signalosc))).^2));
        if disp
            figure
            subplot(2,1,1)
            plot(timeReconS,(Signal/max(abs(Signal))),'DisplayName','CompleteSignal')
            hold on
            plot(timeReconS,InS/max(abs(InS)),'DisplayName','SampledInput')
            plot(timeReconS,Signalosc/max(abs(Signalosc)),'DisplayName','VibrationOnly')
            hold off
            ylabel('normalized signals')
            legend('show')
            subplot(2,1,2)
            plot(timeReconS,(abs(Signal)/max(abs(Signal))),'DisplayName','CompleteSignal')
            hold on
            plot(timeReconS,abs(InS)/max(abs(InS)),'DisplayName','Sampled Input')
            hold off
            ylabel('normalized & rectified signals')
            xlabel('time [ms]')
            legend('show')
            
        end        
    end
end