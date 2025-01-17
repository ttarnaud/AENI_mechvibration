function [VRMAT,VRMATDOI,VRMATDOIstat,VRMATOSC,VRMATOSCstat,VRMATstatnoise,extra_opt_intm]=...
    SimDipoleOsc_IsP(VRMATsP,VRMATsPstat,Tend,CSource,CSink,USperiod,idxOscDip,...
    dpI_time,dpI_space,dps_thisrun,d,meanI,stdI,Ifun,all_osc_flag,varargin)
% with this function we determine the potentials measured at the POIs.
% based on already determined potentials during one period 
%Output:  VRMAT: sum of all signals
%         VRMATDOI: only signal generated by DOI
%         VRMATDOIstat: signal of DOI if would stand still
%         VRMATosc: only signal of oscillators (not DOI)
%         VRMATOSCstat: only signal of oscillators (not DOI) if would stand still 
%         VRMATstatnoise:static noise if everything would stand still
%         
%        VRMAT =  staticomponents + vibrcomponents of all signals(DOI+OSC+Static)
%        VRMATDOI = static + vibr of DOI
%        VRMATDOIstat = static of DOI
%        VRMATOSC = static+vibr of OSC
%        VRMATOSCstat = statoc of OSC
%        VRMATstat noise = stat of all signals
%
%Input:   VRMATsP [POIs x Length(tsimUS) x length (CSourse)] same data as
%           VRMAT
%         VRMATsPstat [POIs x Length(tsimUS) x length (CSourse)] same data as
%           VRMATstat noise
%         
%         Tend ==> end of simulation
%         CSource: Currentsource locations.
%         CSink: Currentsink locaitons
%         AmpI:  Currents throuhg dipoles in uA!!!
%         USwave: vibrations of each dipole
%         USperiod: period of one vibration
%         idxOscDip: idices of vibrators (not all are vibrating perse ==>
%         decrease comp time)



ONOFFstr = {'OFF','ON'};
%default settings
display_flag = 0;        %display input settings and progresss
Show_sphere = 'OFF';     % create sphere with positions, 'wpoi','nopoi','OFF'
SphereRes = 30;          % number of points in sphere
sigma = 0.33; %S/m       % conductivity of brain matter
RSphere = 0.07; %m       % size of brain
resUS = 100;             % resUS wave
POIs = [];               %declaration POIs
scale_flag = 1;            % scale POIs to outer sphere
idxDOIindice = 1;          % idx of DOI
pulsed_flag = 0;
pulsed_prp = 0;
pulsed_dc = 0;


if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'display'))
            display_flag = varargin{find(strcmpi(varargin,'display'))+1};
        end
        if any(strcmpi(varargin,'RSphere'))
            RSphere = varargin{find(strcmpi(varargin,'RSphere'))+1};
        end
        if any(strcmpi(varargin,'sigma'))
            sigma = varargin{find(strcmpi(varargin,'sigma'))+1};
        end
        if any(strcmpi(varargin,'resUS'))
            resUS = varargin{find(strcmpi(varargin,'resUS'))+1};
        end
        if any(strcmpi(varargin,'idxI'))
            idxI = varargin{find(strcmpi(varargin,'idxI'))+1};
        end
        if any(strcmpi(varargin,'POI'))
            idxPOI = find(strcmpi(varargin,'POI'))+1;
            POIs = varargin{idxPOI};
        end
        if any(strcmpi(varargin,'ShowSphere'))
            Show_sphere = varargin{find(strcmpi(varargin,'ShowSphere'))+1};
        end
        if any(strcmpi(varargin,'scale'))
            scale_flag = varargin{find(strcmpi(varargin,'scale'))+1};
        end
        if any(strcmpi(varargin,'DOI'))
            idxDOIindice = varargin{find(strcmpi(varargin,'DOI'))+1};
        end
        if any(strcmpi(varargin,'pulsed'))
            pulsedInfo = varargin{find(strcmpi(varargin,'pulsed'))+1};
            pulsed_flag = pulsedInfo.flag;
            pulsed_prp = pulsedInfo.prp;
            pulsed_dc = pulsedInfo.dc;
        end
        
    end
else
    error('incorrect input')
end
if display_flag
    fprintf('SETTINGS:\n')
    fprintf('\t Show sphere: %s\n',Show_sphere);
    fprintf('\t Sphere Resolution: %d\n',SphereRes);
    fprintf('\t Sigma: %3.2g\n',sigma);
    fprintf('\t Brain size: %3.2g\n',RSphere);
    fprintf('\t ResUS: %d\n',resUS);
    fprintf('\t Scale POIs: %s\n',ONOFFstr{scale_flag+1});   
end
try
    Aus = varargin{find(strcmpi(varargin,'Aus'))+1};
catch
    error('Aus needs to be included')
end

% simulate first single oscilation cycle: only necessary to simulate one
% period
dt = USperiod/resUS;
TsimUS = 0:dt:USperiod-dt;
if isempty(TsimUS)
    TsimUS = 1;
    disp('resolution not high enough for US modeling')
end
Tsim = 0:dt:Tend;
f = ceil(length(Tsim)/length(TsimUS));   %periods needed to contain whole Tsim

% Generate I(t)
[I,extra_opt_intm] = getIarray(dpI_time,dpI_space,Tsim,dps_thisrun,d,meanI,stdI,[]);


% Check if DOI is in this run if so assign IOI as current
DOIflag = idxDOIindice<=dps_thisrun;
if any(DOIflag)
    IOI = meanI.*Ifun(Tsim);
    I(idxDOIindice,:) = IOI;
end


if pulsed_flag
    pulseidx = double(mod(Tsim,pulsed_prp)<=(pulsed_dc*pulsed_prp));
    pulseidx(f*length(TsimUS)) = 0;
else
    pulseidx = ones(1,f*length(TsimUS));
end

% plot initial settings
%prior to simulation a plot is shown of where dipoles are located
if strcmpi(Show_sphere,'nopoi') || strcmpi(Show_sphere,'wpoi')
%     if strcmpi(SolutionType,'3Sphere')
%         RSpherePLOT = RSphere + tSkull + tScalp;
%     elseif strcmpi(SolutionType,'4Sphere')
%         RSpherePLOT = RSphere + tSkull + tScalp + tAir;
%     else
%         RSpherePLOT = RSphere;
%     end
    if isempty(get(groot,'Children'))
        Fignr = 1;
    else
        Fignr = get(gcf,'number')+1;
    end
    
    switch lower(Show_sphere)        
        case 'nopoi'
            varargin2 = horzcat(varargin,{'DOI',idxDOIindice,'OSC',idxOscDip,'scale',scale_flag});
            if exist('idxPOI','var')
               varargin2{idxPOI} = []; %exclude idsPOI => not plotted
            end            
            PotentialSphere_Multi(CSource,CSink,10,sigma,RSphere,SphereRes,Fignr,varargin2); %I(:,floor(size(I,2)/2))
        case 'wpoi'
            %POIsPLOT = POIs.*RSpherePLOT./RSphere;
            varargin2 = horzcat(varargin,{'DOI',idxDOIindice,'OSC',idxOscDip,'scale',scale_flag});
            PotentialSphere_Multi(CSource,CSink,10,sigma,RSphere,SphereRes,Fignr,varargin2);
    end
end



        

% VRMAT = zeros(size(POIs,1),length(Tsim));
% VRMATDOI = zeros(size(POIs,1),length(Tsim));
% VRMATDOIstat = zeros(size(POIs,1),length(Tsim));
% VRMATOSC = zeros(size(POIs,1),length(Tsim));
% VRMATOSCstat = zeros(size(POIs,1),length(Tsim));
% VRMATstatnoise = zeros(size(POIs,1),length(Tsim));
%warning('parfor off')

idx = 1:size(CSource,1);
idx2 = find((1:size(CSource,1))==idxOscDip);
if ~isempty(idx)
    if isempty(idx2)        
        VRstatnoise = VRMATsPstat(idx,:,:);        
        VRstatnoise = repmat(VRstatnoise,1,f,1); 
        VRMATstatnoise = permute(sum(I(idx,:).*VRstatnoise(:,1:length(I(1,:)),:),1),[3,2,1]);
        VRMAT =  VRMATstatnoise;
    elseif numel(idx2)~=numel(idx)        
        VRstatnoise = VRMATsPstat(idx,:,:);
        VRVibration = VRstatnoise;
        VRVibration(idx2,:,:) = Aus(idx2).*(VRMATsP(idx2,:,:)-VRMATsPstat(idx2,:,:))+VRMATsPstat(idx2,:,:);
        VRVibration = repmat(VRVibration,1,f,1);
        VRstatnoise = repmat(VRstatnoise,1,f,1);
        if pulsed_flag
            VRVibration = pulseidx.*(VRVibration-VRstatnoise)+VRstatnoise;
        end
        VRMAT =  permute(sum(I(idx,:).*VRVibration(:,1:length(I(1,:)),:),1),[3,2,1]);
        VRMATstatnoise = permute(sum(I(idx,:).*VRstatnoise(:,1:length(I(1,:)),:),1),[3,2,1]);
    else
        VRVibration = Aus(idx).*(VRMATsP(idx,:,:)-VRMATsPstat(idx,:,:))+VRMATsPstat(idx,:,:);
        VRstatnoise = VRMATsPstat(idx,:,:);
        VRVibration = repmat(VRVibration,1,f,1);
        VRstatnoise = repmat(VRstatnoise,1,f,1);
        if pulsed_flag
            VRVibration = pulseidx.*(VRVibration-VRstatnoise)+VRstatnoise;
        end
        VRMAT =  permute(sum(I(idx,:).*VRVibration(:,1:length(I(1,:)),:),1),[3,2,1]);
        VRMATstatnoise = permute(sum(I(idx,:).*VRstatnoise(:,1:length(I(1,:)),:),1),[3,2,1]);
    end
else
    VRMAT = zeros(size(POIs,1),length(Tsim));
    VRMATstatnoise = zeros(size(POIs,1),length(Tsim));
end



idx = find((1:size(CSource,1)==idxDOIindice));
if ~isempty(idx)
VRVibration = Aus(idx).*(VRMATsP(idx,:,:)-VRMATsPstat(idx,:,:))+VRMATsPstat(idx,:,:);
VRstatnoise = VRMATsPstat(idx,:,:);
VRVibration = repmat(VRVibration,1,f,1);
VRstatnoise = repmat(VRstatnoise,1,f,1);
if pulsed_flag
    VRVibration = pulseidx.*(VRVibration-VRstatnoise)+VRstatnoise;
end
VRMATDOI =  permute(sum(I(idx,:).*VRVibration(:,1:length(I(1,:)),:),1),[3,2,1]);
VRMATDOIstat = permute(sum(I(idx,:).*VRstatnoise(:,1:length(I(1,:)),:),1),[3,2,1]);
else
    VRMATDOI = zeros(size(POIs,1),length(Tsim));
    VRMATDOIstat = zeros(size(POIs,1),length(Tsim));
end

if  all_osc_flag
    VRMATOSC = VRMAT-VRMATDOI;
    VRMATOSCstat = VRMATstatnoise-VRMATDOIstat;
else
    idx = find((1:size(CSource,1))==idxOscDip & ~(1:size(CSource,1)==idxDOIindice));
    if ~isempty(idx)
        VRVibration = Aus(idx).*(VRMATsP(idx,:,:)-VRMATsPstat(idx,:,:))+VRMATsPstat(idx,:,:);
        VRstatnoise = VRMATsPstat(idx,:,:);
        VRVibration = repmat(VRVibration,1,f,1);
        VRstatnoise = repmat(VRstatnoise,1,f,1);
        if pulsed_flag
            VRVibration = pulseidx.*(VRVibration-VRstatnoise)+VRstatnoise;
        end
        VRMATOSC = permute(sum(I(idx,:).*VRVibration(:,1:length(I(1,:)),:),1),[3,2,1]);
        VRMATOSCstat = permute(sum(I(idx,:).*VRstatnoise(:,1:length(I(1,:)),:),1),[3,2,1]);
    else
        VRMATOSC = zeros(size(POIs,1),length(Tsim));
        VRMATOSCstat = zeros(size(POIs,1),length(Tsim));
    end
end
% 
% for idp = 1:size(CSource,1) % loop over dipoles
%     % Check if dipole lies within brain region
%     if norm(CSource(idp,:))>=RSphere || norm(CSink(idp,:))>=RSphere
%         error('dipole not within brain region')
%     end
%     VRVibration = Aus(idp)*VRMATsP(:,:,idp); %contains actually static and vibration component
%     VRstatnoise = Aus(idp)*VRMATsPstat(:,:,idp); %contains only static component
%     if display_flag
%         pdpsstep0 = 10^floor(log10(size(CSource,1)-1)); % not to print each dipole
%         if idp == 1
%             pdpsstep = pdpsstep0;
%             fprintf('\t Dipole %i \n',idp)
%         elseif idp >= pdpsstep
%             fprintf('\t Dipole %i \n',idp)
%             pdpsstep = pdpsstep+pdpsstep0;
%         end
%     end
%     VRVibration = repmat(VRVibration,1,f);
%     VRstatnoise = repmat(VRstatnoise,1,f);
%     VRMAT = VRMAT + I(idp,:).*VRVibration(:,1:length(I(idp,:)));
%     VRMATstatnoise = VRMATstatnoise + I(idp,:).*VRstatnoise(:,1:length(I(idp,:)));
%     if any(idp == idxDOIindice)
%         VRMATDOI = VRMATDOI + I(idp,:).*VRVibration(:,1:length(I(idp,:)));
%         VRMATDOIstat = VRMATDOIstat+I(idp,:).*VRstatnoise(:,1:length(I(idp,:)));
%     end
%     if any(idp==idxOscDip) && ~any(idp == idxDOIindice)
%         VRMATOSC = VRMATOSC + I(idp,:).*VRVibration(:,1:length(I(idp,:)));
%         VRMATOSCstat = VRMATOSCstat + I(idp,:).*VRstatnoise(:,1:length(I(idp,:)));
%     end
% end
% save memory
VRVibration = []; VRstatnoise = [];
VRMATsP = []; VRMATsPstat = [];
I = [];
VRMAT = VRMAT.*1e6; % outcome is given in �V
VRMATDOI = VRMATDOI.*1e6;
VRMATOSC = VRMATOSC.*1e6;
VRMATDOIstat = VRMATDOIstat.*1e6;
VRMATOSCstat = VRMATOSCstat.*1e6;
VRMATstatnoise = VRMATstatnoise*1e6;

    

end
