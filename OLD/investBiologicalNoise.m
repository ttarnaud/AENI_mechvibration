function [RMS,SNR,TSTOP,Out] = investBiologicalNoise(varargin)
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
TSTART = tic;
%% Initialise default variables
PLOT = 1;
Display = 0;
plotIarray_flag = 0;
plotUSwaves_flag = 1;
showSphere='';
ParallelCompute_flag = 0;
ThermalNoiseAmp = 0;
scale_flag = 1;
input_posDp_flag = 0;
sel_oscil_from_gendPos_flag = 0;
S4l_flag = 0;
inerp3_method = 'linear';
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
resUS = 20;
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
        if any(strcmpi(varargin,'showSphere'))
            showSphere = varargin{find(strcmpi(varargin,'showSphere'))+1};
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
        if any(strcmpi(varargin,'Phaseus'))
            Phaseus = varargin{find(strcmpi(varargin,'Phaseus'))+1};
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
            input_posDp_flag = 1;
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
        if any(strcmpi(varargin,'ROI_OSC'))
            ROI_OSC = varargin{find(strcmpi(varargin,'ROI_OSC'))+1};
            sel_oscil_from_gendPos_flag = 1;
        end
        if any(strcmpi(varargin,'Sim4life'))
            S4l_Input = varargin{find(strcmpi(varargin,'Sim4life'))+1};
            S4l_flag = S4l_Input.S4l_flag;
        end
    end
end

%% Update simulation settings
% US-wave
% Input is velocity generated by sim4life
if S4l_flag
    fus = S4l_Input.fus;
    wus = 2*pi*fus;
else
wus = 2*pi*fus; %rad/s
Dirus = Dirus/norm(Dirus);
USwave_varlengths = [size(Dirus,1), size(Aus,1), length(wus), size(Phaseus,1)];
end
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
    %showSphere = 'off';
    Display = 1;
end
% Select Sources and sinks + select vibrator
if strcmpi(SolutionType,'Mouse4Sphere~fair0')
dRSphere = 0.0005;
else
dRSphere = 0.005;
end
if ~input_posDp_flag
posDp = [0,0,RSphere-dRSphere];
end

if strcmpi(dpDistribution,'s4l_hotspot')
    [CSource,CSink] = GendpPos(dpDistribution,dpOrientation,d,RSphere-dRSphere,Totaldps-size(posDp,1),'S4l_Input',S4l_Input);
else    
    [CSource,CSink] = GendpPos(dpDistribution,dpOrientation,d,RSphere-dRSphere,Totaldps-size(posDp,1));
end
Totaldps = size(CSource,1) + size(posDp,1);

if sel_oscil_from_gendPos_flag
    CSource = vertcat(posDp+d/2*OrienDipole,CSource);
    CSink = vertcat(posDp-d/2*OrienDipole,CSink);
    OSCindices = find(inhull(CSource,ROI_OSC)&inhull(CSink,ROI_OSC));
elseif S4l_flag
    CSource = vertcat(posDp+d/2.*OrienDipole,CSource);
    CSink = vertcat(posDp-d/2.*OrienDipole,CSink);
    posDPall = (CSource+CSink)/2;
    % check if all dipoles are in brain simulation domain
    if any(vecnorm(CSource,2,2)>S4l_Input.RBrain) || any(vecnorm(CSink,2,2)>S4l_Input.RBrain)
        error('dipoles are not in simulation domain s4l')
    end
    switch S4l_Input.loc.unit
        case 'm'
            f_loc = 1;
        case 'mm'
            f_loc = 1e-3;
        otherwise
            error('unit not implemented')
    end
    S4l_Input.loc.X = S4l_Input.loc.X*f_loc;
    S4l_Input.loc.Y = S4l_Input.loc.Y*f_loc;
    S4l_Input.loc.Z = S4l_Input.loc.Z*f_loc;
    S4l_maxX = max(S4l_Input.loc.X(:)); S4l_minX = min(S4l_Input.loc.X(:));
    S4l_maxY = max(S4l_Input.loc.Y(:)); S4l_minY = min(S4l_Input.loc.Y(:));
    S4l_maxZ = max(S4l_Input.loc.Z(:)); S4l_minZ = min(S4l_Input.loc.Z(:));
    S4l_maxloc = [S4l_maxX, S4l_maxY, S4l_maxZ];
    S4l_minloc = [S4l_minX, S4l_minY, S4l_minZ];
    % output of s4l can be cropped for memory reasons find which are
    % contained within boundary
    OSCindices = find(~(any([S4l_maxloc-CSource]<0,2)|any([S4l_maxloc-CSink]<0,2)|...
        any([S4l_minloc-CSink]>0,2)|any([S4l_minloc-CSource]>0,2)));
    
    vampx = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.ampx,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),inerp3_method);
    vampy = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.ampy,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),inerp3_method);
    vampz = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.ampz,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),inerp3_method); 
    vphasex = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.phasex,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),inerp3_method);
    vphasey = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.phasey,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),inerp3_method);
    vphasez = interp3(S4l_Input.loc.X,S4l_Input.loc.Y,S4l_Input.loc.Z,S4l_Input.veloc.phasez,...
        posDPall(:,1),posDPall(:,2),posDPall(:,3),inerp3_method);
    
else
    posOsc = [posDp]; % position of oscilators
    CSource = vertcat(posDp+d/2.*OrienDipole,CSource);
    CSink = vertcat(posDp-d/2.*OrienDipole,CSink);
    OSCindices = 1:size(posOsc,1);    
end

DOIindice = 1;
if Display
    fprintf('\n Dipole of interest is in position %i \n',DOIindice)
end

if ~any(DOIindice==OSCindices)
    warning('DOI not an oscillator')
end
%Create USwave_all
USwave_allsame_flag = 0;
if ~S4l_flag
    if length(OSCindices) == 1
        % if dimensions of USwaves bigger display warning
        if any(USwave_varlengths-1)
            warning('USwaves more dimension than selected oscillators')
        end
        % because there is only one oscillator we take always first
        % setting set as USwave
        USwave_all = @(t,idx) Dirus(idx==OSCindices,:).*Aus(idx==OSCindices,:).*sin(wus(idx==OSCindices).*t+Phaseus(idx==OSCindices,:));
        USwave_allsame_flag = 1;
    elseif isempty(OSCindices)
        warning('no oscillator selected')
        USwave_all = @(t,idx) Dirus.*Aus.*sin(wus.*t+Phaseus);
        % eventhough no oscillator selected still declare USwave_all
        plotUSwaves_flag = 0;
    else
        % check if dimensions USwave == 1
        if ~any(USwave_varlengths-1)
            warning('all oscillators have same USwave')
            USwave_all = @(t,idx) Dirus.*Aus.*sin(wus.*t+Phaseus);
            USwave_allsame_flag = 1;
        else
            % extract max dimension and check if matches number of oscillators
            [max_USwave_varlengths, idx_mUv] = max(USwave_varlengths);
            if length(OSCindices)>max_USwave_varlengths
                error('size oscillators does not match dimensions USwaves')
            else
                % throw warning if less osccilators than dimensions
                if length(OSCindices)<max_USwave_varlengths
                    warning('USwaves more dimension than selected oscillators')
                end
                % check if all same dimensions if so do nothing
                if  any(USwave_varlengths-USwave_varlengths(1))
                    if any(USwave_varlengths(~idx_mUv)-1)
                        error('length of US variables should either be the same for all or one should be bigger while others are one')
                    else
                        switch idx_mUv
                            case 1
                                Aus = repmat(Aus,max_USwave_varlengths,1);
                                wus = repmat(wus,max_USwave_varlengths,1);
                                Phaseus = repmat(Phaseus,max_USwave_varlengths,1);
                            case 2
                                Dirus = repmat(Dirus,max_USwave_varlengths,1);
                                wus = repmat(wus,max_USwave_varlengths,1);
                                Phaseus = repmat(Phaseus,max_USwave_varlengths,1);
                            case 3
                                Dirus = repmat(Dirus,max_USwave_varlengths,1);
                                Aus = repmat(Aus,max_USwave_varlengths,1);
                                Phaseus = repmat(Phaseus,max_USwave_varlengths,1);
                            case 4
                                Dirus = repmat(Dirus,max_USwave_varlengths,1);
                                Aus = repmat(Aus,max_USwave_varlengths,1);
                                wus = repmat(wus,max_USwave_varlengths,1);
                        end
                    end
                end
                USwave_interm = @(t,idx) Dirus(idx,:).*Aus(idx,:).*sin(wus(idx).*t+Phaseus(idx,:));
                USwave_all = @(t,idx) USwave_interm(t,find(idx==OSCindices));
                
            end
        end
    end
else
    USwave_all = @(t,idx) 1/wus.*[vampx(idx).*sin(wus.*t+vphasex(idx)),vampy(idx).*sin(wus.*t+vphasey(idx)),vampz(idx).*sin(wus.*t+vphasez(idx))];
    % velocity vectors are given als input v(r,t) = Re(v(r)exp(jwt))==>
    % v(r,t) =
    % [|vx(r)|cos(wt+thetax),|vy(r)|cos(wt+thetay),|vz(r)|cos(wt+thetaz)]
    % x(r,t) =  int(v(r,t),t) = 1/w*
    % [|vx(r)|sin(wt+thetax),|vy(r)|sin(wt+thetay),|vz(r)|sin(wt+thetaz)]
    
end

TSTOP.startSim = toc(TSTART);

%% Start simulations
VR = zeros(size(POIs,1),length(Tsim));
VRDOI = zeros(size(POIs,1),length(Tsim));
VRosc = zeros(size(POIs,1),length(Tsim));
VRnoise = zeros(size(POIs,1),length(Tsim));
irun_end = ceil(Totaldps/dps_run);
if (Totaldps/dps_run)-round((Totaldps/dps_run),0)~=0
    warning('dps preferably devider of Totaldps')
end
if logical(mod(dps_run,4))
    warning('dps_run preferably multiple of 4')
end

for irun = 1:irun_end
    if irun==irun_end
        dpsrun_end = Totaldps;
    else
        dpsrun_end = irun*dps_run;
    end
    dps_thisrun = dpsrun_end-(irun-1)*dps_run;
    dpsrun_start = (irun-1)*dps_run;
% Generate I(t)
Iarray = getIarray(dpI_time,dpI_space,Tsim,dps_thisrun,meanI,d,stdI,[]);
TSTOP.Itgen(irun) = toc(TSTART);
DOIflag = DOIindice>dpsrun_start & DOIindice<=dpsrun_end;
if any(DOIflag)
    IOI = meanI.*Ifun(Tsim);
    idxDOIindice = DOIindice(DOIflag)-dpsrun_start;
    Iarray(idxDOIindice,:) = IOI;
end

if any(OSCindices>dpsrun_start & OSCindices<=dpsrun_end)
    idxOscDip = OSCindices-dpsrun_start;
    USwave = @(t,idx) USwave_all(t,idx+dpsrun_start);
else 
    idxOscDip = dps_thisrun+1;
    USwave = @(t,idx) USwave_all(t,1);
end


idxdps = dpsrun_start+1:dpsrun_end;
if plotIarray_flag
    figure
    PlotIarray(Iarray,Tsim)
end
if Display
disp(['start simulation ',num2str(irun),' of ',num2str(irun_end)])
end
Settings = horzcat(Options,{'POI',POIs,'ShowSphere',showSphere,'resUS',resUS,'ParallelCompute',ParallelCompute_flag,'scale',scale_flag,'Display',Display});
[VRrun,VRDOI_run,VRosc_run,VRnoise_run] = SimDipoleOsc(Tend,CSource(idxdps,:),CSink(idxdps,:),Iarray,USwave,1/fus,idxOscDip,Settings);
TSTOP.VRgen(irun) = toc(TSTART);
VR = VR+VRrun; VRDOI = VRDOI + VRDOI_run; VRosc = VRosc+VRosc_run; VRnoise = VRnoise+VRnoise_run;
end
ThermalNoise = ThermalNoiseAmp.*rand(size(VR));
VR = VR + ThermalNoise; VRDOI = VRDOI; VRosc = VRosc+ThermalNoise; VRnoise = VRnoise + ThermalNoise;
TSTOP.endsim = toc(TSTART);

%% plotFigure
% function plotinfo(Tsim,Dirus,Aus,wus,Phaseus,POIs,VR,Alphafun,Tau)
%VR = VROLD+rednoise(size(VR,1),size(VR,2));
POIsidx = 1:min(size(POIs,1),3);

inputSignal = [Tsim;IOI];
if PLOT
figure

subplot(4,3,2);
plot([Tsim(1:100)*1e6.*ones(3,1)]',USwave_all(Tsim(1:100)',1)*1e9)
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
%% plot USwaves
xyzplot = 0;
if plotUSwaves_flag
%     if length(OSCindices)>8
%         figdim = [3,3];
%         selOSCindices = OSCindices([1,randperm(length(OSCindices)-1,7)+1],:);
%     else
%         if length(OSCindices)>5
%         figdim = [3,3];
%         elseif length(OSCindices)>3
%         figdim = [2,3];
%         elseif length(OSCindices)>1
%         figdim = [2,2];
%         else figdim = [1,2]; 
%         end
%         selOSCindices = OSCindices;
%     end
    if length(OSCindices)>3
        figdim = [2,2];
        selOSCindices = OSCindices([1,randperm(length(OSCindices)-1,2)+1],:);
    else
        if length(OSCindices)>1
        figdim = [2,2];
        else
            figdim = [1,2];
        end
        selOSCindices = OSCindices;
    end
    figure
    fignr = get(gcf,'number');
    subplot(figdim(1),figdim(2),1)
    varargin2 = {'DOI',DOIindice,'OSC',1:length(selOSCindices),'scale',scale_flag,'ndipole_flag',1};
    PotentialSphere_Multi(CSource(selOSCindices,:),CSink(selOSCindices,:),10,sigma,RSphere,20,fignr,varargin2);
    for iplot = 2:min(figdim(1)*figdim(2),length(OSCindices)+1)
        
        subplot(figdim(1),figdim(2),iplot)
        if xyzplot
        plot([Tsim(1:100)*1e6.*ones(3,1)]',USwave_all(Tsim(1:100)',selOSCindices(iplot-1))*1e9)
        legend({'x component','y component','z component'})
        ylabel('displacement [nm]')
        xlabel('Time [µs]')
        else
        if S4l_flag
            dAmpx = vampx(selOSCindices(iplot-1))/wus;
            dAmpy = vampy(selOSCindices(iplot-1))/wus;
            dAmpz = vampz(selOSCindices(iplot-1))/wus;
            phasex = vphasex(selOSCindices(iplot-1));
            phasey = vphasey(selOSCindices(iplot-1));
            phasez = vphasez(selOSCindices(iplot-1));
        else
            if USwave_allsame_flag
                idx_sel = 1;
            else
                idx_sel = selOSCindices(iplot-1);
            end
            if size(Aus,2)==3
                dAmpx = Aus(idx_sel ,1).*Dirus(idx_sel ,1);
                dAmpy = Aus(idx_sel ,2)*Dirus(idx_sel ,2);
                dAmpz = Aus(idx_sel ,3)*Dirus(idx_sel ,3);
            elseif size(Aus,2)==1
                dAmpx = Aus(idx_sel ,1).*Dirus(idx_sel ,1);
                dAmpy = Aus(idx_sel ,1).*Dirus(idx_sel ,2); 
                dAmpz = Aus(idx_sel ,1).*Dirus(idx_sel ,3);
            else 
                error('possible?')
            end
            if size(Phaseus,2)==3
                phasex = Phaseus(idx_sel ,1);
                phasey = Phaseus(idx_sel ,2);
                phasez = Phaseus(idx_sel ,3);
            elseif size(Aus,2)==1
                phasex = Phaseus(idx_sel ,1);
                phasey = phasex; phasez = phasey;
            else 
                error('possible?')
            end
        end
        genvibrplot(CSource(selOSCindices(iplot-1),:),CSink(selOSCindices(iplot-1),:),...
            dAmpx,dAmpy,dAmpz,phasex,phasey,phasez,iplot-1)
        end
    end
    
end

%% Signal reconstruction and metric determination
% define if signal needs to be mirrored
mirror = zeros(size(POIs,1),1);
DOIspos = (CSource(DOIindice,:)+CSink(DOIindice,:))/2;
dDOIPOI0 = vecnorm(DOIspos+USwave_all(0,DOIindice )-POIs,2,2);
dDOIPOI1 = vecnorm(DOIspos+USwave_all(1/10*1/fus,DOIindice )-POIs,2,2);
dDOIPOI2 = vecnorm(DOIspos+USwave_all(2/10*1/fus,DOIindice )-POIs,2,2);
dCSPOI0 = vecnorm(CSource(DOIindice,:)+USwave_all(0,DOIindice )-POIs,2,2);
%dCSPOI0 = vecnorm(CSource(DOIindice,:)-POIs,2,2);
dCSPOI1 = vecnorm(CSource(DOIindice,:)+USwave_all(1/10*1/fus,DOIindice )-POIs,2,2);
dCSPOI2 = vecnorm(CSource(DOIindice,:)+USwave_all(2/10*1/fus,DOIindice )-POIs,2,2);
dCSinkPOI0 = vecnorm(CSink(DOIindice,:)+USwave_all(0,DOIindice )-POIs,2,2);
%dCSinkPOI0 = vecnorm(CSink(DOIindice,:)-POIs,2,2);
dCSinkPOI1 = vecnorm(CSink(DOIindice,:)+USwave_all(1/10*1/fus,DOIindice )-POIs,2,2);
dCSinkPOI2 = vecnorm(CSink(DOIindice,:)+USwave_all(2/10*1/fus,DOIindice )-POIs,2,2);
if ~any(dDOIPOI0-dDOIPOI1)
    error('mirroring detection false')
end
%mirror_idx = dDOIPOI1>dDOIPOI0 & dCSPOI0<dCSinkPOI0 | dDOIPOI1<dDOIPOI0 & dCSPOI0>dCSinkPOI0;
mirror_idx = dCSPOI0>dCSinkPOI0;
mirror(mirror_idx) = 1;
if PLOT
figure
for iplot = 1:min(size(POIs,1),3)
    tvals = linspace(0,1/fus,100)';
    dDOIPOIvals = vecnorm(DOIspos+USwave_all(tvals,DOIindice )-POIs(iplot,:),2,2);
    mean_dDPv = mean(dDOIPOIvals);
    dDPv = dDOIPOIvals-mean_dDPv;
    phase = atan2(sin(wus.*tvals(2)),dDPv(2)/dDPv(1)-cos(wus*tvals(2)));
    amplitude = dDPv(1)/sin(phase);
    
    subplot(3,1,iplot)
    plot(tvals,dDOIPOIvals)
    hold on
    %yyaxis right
    plot(tvals,amplitude.*sin(wus.*tvals+phase)+mean_dDPv,'DisplayName',['phase: ',num2str(phase*180/pi),' ampl: ',num2str(amplitude)])
    legend('show')
    hold off
end
end
% determine theta of USwave
tvals = linspace(0,1/fus,1e8)'; tvals = tvals(1:end-1);
USwave_sel = USwave_all(tvals,DOIindice);
USwave_sel = USwave_sel-mean(USwave_sel);
amp_vals = max(USwave_sel);
USwave_sel = USwave_sel(1:floor(length(tvals)/2),:);
USwave_t0 = USwave_sel(1,:);
dUSwave_t1t0= USwave_sel(2,:)-USwave_sel(1,:);

[~,cross_zero] = min(abs(USwave_sel));
theta_vals(1) = wus*(1/(2*fus)*double(USwave_t0(1)>0)-tvals(cross_zero(1)));
theta_vals(2) = wus*(1/(2*fus)*double(USwave_t0(2)>0)-tvals(cross_zero(2)));
theta_vals(3) = wus*(1/(2*fus)*double(USwave_t0(3)>0)-tvals(cross_zero(3)));
theta_vals(cross_zero==1) = 0; theta_vals(cross_zero==1 & dUSwave_t1t0<0) = pi
% USwave_t0 = USwave_sel(t0);USwave_t1 = USwave_sel(t1);
% theta_vals = atan2(sin(wus.*t1).*ones(1,3),(USwave_t1./USwave_t0-cos(wus.*t1)));
% % nanmean is for special case when theta is zero ==> sin(theta) is 0 gives
% % devision by zero = nan
% amp_vals = nanmean(vertcat(abs(USwave_t0./sin(theta_vals)),abs(USwave_t1./sin(wus.*t1+theta_vals))));
avg_vect = nansum(amp_vals.*exp(complex(0,theta_vals)))./nansum(amp_vals);
theta_global = atan2(imag(avg_vect),real(avg_vect));
fprintf('\nphase of US wave targeting DOI is: %5.2f°\n', theta_global*180/pi)
% Redefine VR
VRnoise_stat = VRnoise;
VRnoise_osc = VRosc-VRDOI;
VRnoise_all = VR-VRDOI;

%reconstruct based on complete signal    
[reconsSignalVR,powerSignalVR,timeReconS] = calcFourierGetsignal(VR,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,PLOT);
% reconstruct based on only DOI
[reconsSignalVRDOI,powerSignalVRDOI] = calcFourierGetsignal(VRDOI,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct osc noise
[reconsSignalVRnoise_osc,powerSignalVRnoise_osc] = calcFourierGetsignal(VRnoise_osc,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct static noise
[reconsSignalVRnoise_stat,powerSignalVRnoise_stat] = calcFourierGetsignal(VRnoise_stat,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);
% reconstruct static noise all
[reconsSignalVRnoise_all,powerSignalVRnoise_all] = calcFourierGetsignal(VRnoise_all,Tsim,POIsidx,fus,fbandwidth,inputSignal,mirror,theta_global,0);

TSTOP.signalreconstruction = toc(TSTART);
RMS = struct;
SNR = struct;
for iPOI=1:length(reconsSignalVR)
    % SNR
    PowerDOI = calcPower(reconsSignalVRDOI{iPOI});
    Powernoise_osc = calcPower(reconsSignalVRnoise_osc{iPOI});
    Powernoise_stat = calcPower(reconsSignalVRnoise_stat{iPOI});
    Powernoise_all = calcPower(reconsSignalVRnoise_all{iPOI});
    SNR.DOI_noiseall(1,iPOI) = 10*log10(PowerDOI/Powernoise_all);
    SNR.DOI_noiseall(2,iPOI) = 10*log10(powerSignalVRDOI(iPOI)/powerSignalVRnoise_all(iPOI));
    SNR.DOI_noiseosc(1,iPOI) = 10*log10(PowerDOI/Powernoise_osc);
    SNR.DOI_noiseosc(2,iPOI) = 10*log10(powerSignalVRDOI(iPOI)/powerSignalVRnoise_osc(iPOI));
    SNR.DOI_noisestat(1,iPOI) = 10*log10(PowerDOI/Powernoise_stat);
    SNR.DOI_noisestat(2,iPOI) = 10*log10(powerSignalVRDOI(iPOI)/powerSignalVRnoise_stat(iPOI));
    %RMS
    [RMS(iPOI).RMS,RMS(iPOI).RMSDOI,~,RMS(iPOI).Q2,RMS(iPOI).Q2DOI] = ...
        calcRMS(reconsSignalVR{iPOI},reconsSignalVRDOI{iPOI},inputSignal,timeReconS{iPOI},mirror(iPOI),PLOT);
    if PLOT
        title(['POI ',num2str(iPOI)])
    end

end
Out.VR = VR;
Out.reconsSignalVR = reconsSignalVR;
Out.timereconS = timeReconS;
Out.reconsSignalVRDOI = reconsSignalVRDOI;
%% functions
    function PlotIarray(Iarray,Tsim)
    subplot(2,1,1)
    for ipArray=1:10
        idx = randi(size(Iarray,1));
        plot(Tsim*1e3,Iarray(idx,:),'k')
        hold on
    end
    hold off
    ylabel('Current in dipole [µA]')
    title('Selection of Currents in dipoles')
    subplot(2,1,2)
    plot(Tsim*1e3,Iarray(end,:),'k')
    hold on
    plot(Tsim*1e3,Iarray(1,:),'r')
    hold off
    xlabel('Time [ms]')
    ylabel('Current in dipole [µA]')
    end
    function Power = calcPower(Signal)
        s0dft = fft(Signal);
        psds0 = (1/(length(Signal))) * abs(s0dft).^2;
        Power = sum(psds0);
    end
    function [RMS,RMSosc,RMSabs, Q2, Q2osc, Q2abs] = calcRMS(Signal,Signalosc,InputSignal,timeReconS,mirror,disp)
        InS = interp1(InputSignal(1,:),InputSignal(2,:),timeReconS);
        % Check if signal lengths are the same
        if length(Signal)~=length(Signalosc) || length(Signal)~=length(InS)
            error('signal lengths not the same')
        end
        Signal = (-1)^mirror.*Signal;
        Signalosc = (-1)^mirror.*Signalosc;
        % redefined normalization in such way that sign diffrence doesn't
        % affect RMS
        %Signal_peakval = max(abs(Signal));
        Signal_peakval = Signal(max(abs(Signal))==abs(Signal));
        %InS_peakval = max(abs(InS));
        InS_peakval = InS(max(abs(InS))==abs(InS));
        %Signalosc_peakval = max(abs(Signalosc));
        Signalosc_peakval = Signalosc(max(abs(Signalosc))==abs(Signalosc));
        
        % RMS 
            RMS = sqrt(mean((Signal/Signal_peakval-InS/InS_peakval).^2));
        %RMSabs
            RMSabs = sqrt(mean((abs(Signal/Signal_peakval)-abs(InS/InS_peakval)).^2));
        %RMSosc
            RMSosc = sqrt(mean((Signal/Signal_peakval-Signalosc/Signalosc_peakval).^2));
        % Q2 
            Q2 = sqrt(mean((Signal/Signal_peakval-InS/InS_peakval).^2))./sqrt(sum((InS/InS_peakval).^2));
        %Q2abs
            Q2abs = sqrt(mean((abs(Signal/Signal_peakval)-abs(InS/InS_peakval)).^2))./sqrt(sum((abs(InS/InS_peakval)).^2));
        %Q2osc
            Q2osc = sqrt(mean((Signal/Signal_peakval-Signalosc/Signalosc_peakval).^2))./sqrt(sum((Signalosc/Signalosc_peakval).^2));
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