function [RMSall,SNRall,Outall,TSTOPall,Totaldps] = calcStrengthCurve_mosc_DBvsnormal(Input)
% Use this function for HPC measurements %function is outdated se
% calcErrorGrid
disp('start function')
if ischar(Input)
Input_temp = load(['./Inputs/',Input,'.mat']);
Input = Input_temp.Input;
end

lim_totaldip_flag = 1;
lim_Amp_flag = 1;
Decomp_struct_flag = 1;
Amp_Input = Input.Amp;
Totaldps = Input.Totaldps; dps_run = Input.dps_run; 
Settings = Input.Settings;  SettingsStr = Input.SettingsStr;
Fignr = Input.Fignr;
HPC_flag = Input.HPC_flag;
S4l_flag = Input.S4l_flag; S4l_fn = Input.Filename_S4l; S4l_loc = Input.Location_S4l;

if lim_totaldip_flag
fprintf('\n !limitation on # dipoles!')
Totaldps = logspace(1,4,7);
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
    Totaldps(i) = round(Totaldps(i),orders(i));
end
dps_run = min(100*ones(size(Totaldps)),Totaldps);
end

Model = Input.Model;
switch lower(Model)
    case 'human_scalp'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.082;
        RBrain = 0.07;
        dRBrain = 0.005;
%     case 'human_cortical'
%         SolutionType = '4SphereS8.7R25';
%         RPOI = 0.07;
%         RBrain = 0.07;
%         dRBrain = 0.005;
%     case 'mouse_fair0'
%         SolutionType = 'Mouse4Sphere~fair0';
%         RPOI = 0.0059;
%         RBrain = 0.0046;
%         dRBrain = 0.0005;
    otherwise 
        error('false input model + remember currently only database human scalp')        
end  

Thresholds = [0.1];
Amp = Amp_Input;

DEBUG = 0;
PLOT_flag = 0;
LOCALRUN = 1;
PLOT_s4l = 0;
Outall = struct();
RMSall = struct();
SNRall = struct();
TSTOPall = struct();

%txtbxstr = gentxtbxstr(Settings);

%% Debug mode
if DEBUG
    disp('DEBUG modus: limited number of dipoles !!!!')
    disp('')
    Totaldps = 300;
    dps_run = 100;
    %Amp_Input = Amp_Input([1,end]);
    HPC_flag = 0;
    PLOT_flag = 1;
    PLOT_s4l = 1;
    Amp = 1e-6;
    S4l_fn = 'ExportedData_SHHM_AW164_cropx_11slices';
    S4l_loc = 'D:\no backup\Sim4Life\3spherical\CenterFocus\Wicasim3_sHHM\';
elseif LOCALRUN
    disp('LOCAL RUN!!!!!')
    disp('___________________________')
    Amp = logspace(2,5,3).*1e-9;
    Totaldps = 2:10;%round(logspace(1,2,2));%round(logspace(0,5+double(isettings>=9 & isettings<=12),6+2*+double(isettings>=9 & isettings<=12)));
    orders = -floor(log10(Totaldps));
    for i=1:length(Totaldps)
        Totaldps(i) = round(Totaldps(i),orders(i));
    end
    dps_run = min(1e2,Totaldps);
    HPC_flag = 0;
    PLOT_flag = 0;
    PLOT_s4l = 0;
    
    S4l_fn = 'ExportedData_SHHM_AW164_cropx_11slices';
    S4l_loc = 'D:\no backup\Sim4Life\3spherical\CenterFocus\Wicasim3_sHHM\';
    if isdeployed
        S4l_loc = './PressureFields\';
        S4l_fn = 'ExportedData_MM_changeACSkull.mat';
    end
end

if lim_Amp_flag
fprintf('\n !limitation on # Amps!\n')
Amp(2:end-1) = [];
end

if HPC_flag
    try
        c = parcluster('local');
        pool = parpool(c.NumWorkers);
    end
    S4l_loc = './PressureFields\';
else
    if ~isdeployed
    addpath(genpath('./Functions'));
    addpath(genpath('./DataBases'));
    addpath('D:\users\rschoeters\Documents\Imec USEEG\Sim4life')
    end
end

disp(['Startingmodel: ',Model]);
[X,Y] = meshgrid(Totaldps,Amp);
Y = Y*1e9;
%fix rng
try
s = load('randomseed.mat');
s = s.s;
catch
    warning('new random seed')
    s = rng;
end

% load pressure field
S4l_input = struct();
[S4l_input.veloc, S4l_input.press, S4l_input.loc, S4l_input.fus, maxv,maxv_pos,maxv_dir] = ...
    calcvelocity(S4l_loc,S4l_fn,RBrain-dRBrain,PLOT_s4l,'snapshot',0,'rotate',0','totalperiods',50,'inittime',48);
S4l_input.loc.unit = 'mm';
if strcmpi(S4l_input.loc.unit,'mm')
    maxv_pos = maxv_pos*1e-3;
end    
S4l_input.S4l_flag = 1;
S4l_input.resolution=Input.S4l_res;
S4l_input.RBrain = RBrain;
t1 = fzero(@(t) vecnorm(maxv_pos+t'*maxv_dir,2,2)-RBrain,[0,1]);
t2 = fzero(@(t) vecnorm(maxv_pos+t'*maxv_dir,2,2)-RBrain,[-1,0]);
ts = [t1,t2];
[~,idxt] = min(abs(ts));
POIs = [1,0,0;zeros(size([0:pi/10:pi/2]')),sin([0:pi/10:pi/2]'),cos([0:pi/10:pi/2]');maxv_pos+ts(idxt)*maxv_dir;1,1,1];
POIs = POIs./vecnorm(POIs,2,2);
POIs = RPOI.*POIs;
Input_AdaptVamp.flag = 1;
Input_AdaptVamp.vmax = maxv;


for iAmp = 1:length(Amp)
    for idps = 1:length(Totaldps)
        rng(s)
        fprintf('\n%s\n Amp: %i/%i, idps: %i/%i\n',datestr(now),iAmp,length(Amp),idps,length(Totaldps))
        Input_AdaptVamp.Aus = Amp(iAmp);
        Input_iBN = horzcat({'posDp',maxv_pos,'OrienDipole',maxv_dir,...
            'Totaldps',Totaldps(idps),'dps_run',dps_run(idps),...
            'Aus',Amp(iAmp),'ParallelCompute',1,'SolutionType',SolutionType,...
            'POIs',POIs,'scale',0,'PLOT',PLOT_flag,'HPC_flag',HPC_flag,'Adapt_vAmp',Input_AdaptVamp,'Sim4life',S4l_input},Settings);
        
        %normal run
        [RMS,SNR,TSTOP,Out] = investBiologicalNoise(horzcat(Input_iBN,{'eval_method','normal'}));
        RMSall(iAmp,idps).RMS = RMS;
        SNRall(iAmp,idps).SNR = SNR;
        Outall(iAmp,idps).Out = Out;
        TSTOPall(iAmp,idps).TSTOP = TSTOP;
        pause(0.25)
        clear('RMS','SNR','TSTOP','Out');
        pause(0.25)
        % database
        rng(s)
        [RMS,SNR,TSTOP,Out] = investBiologicalNoise(horzcat(Input_iBN,{'eval_method','database_f'}));
        RMSall(iAmp,idps).RMS_DB = RMS;
        SNRall(iAmp,idps).SNR_DB = SNR;
        Outall(iAmp,idps).Out_DB = Out;
        TSTOPall(iAmp,idps).TSTOP_DB = TSTOP;
        pause(0.25)
        clear('RMS','SNR','TSTOP','Out');
        pause(0.25)

    end
end

%%
Itgen_mat = zeros(length(Amp),length(Totaldps),max(ceil(Totaldps./dps_run)));
VRgen_mat = zeros(length(Amp),length(Totaldps),max(ceil(Totaldps./dps_run)));
Itgen_DB_mat = zeros(length(Amp),length(Totaldps),max(ceil(Totaldps./dps_run)));
VRgen_DB_mat = zeros(length(Amp),length(Totaldps),max(ceil(Totaldps./dps_run)));
RMS_info = {};
RMS_DB_info = {};
liPOIs = size(POIs,1);
if Decomp_struct_flag
    for iAmp = 1:length(Amp)
        for idps = 1:length(Totaldps)
            for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                % Decompose RMS structure
                RMS_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMS;
                RMS_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).RMS;
                RMSDOI_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSDOI;
                RMSDOI_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).RMSDOI;
                Q2_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2;
                Q2_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).Q2;
                Q2DOI_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2DOI;
                Q2DOI_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).Q2DOI;
                RMSptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSptherm;
                RMSptherm_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).RMSptherm;
                RMSDOIptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSDOIptherm;
                RMSDOIptherm_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).RMSDOIptherm;
                Q2ptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2ptherm;
                Q2ptherm_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).Q2ptherm;
                Q2DOIptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2DOIptherm;
                Q2DOIptherm_DB_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS_DB(iPOIs).Q2DOIptherm;
                RMS_info(iAmp,idps,iPOIs) = {RMSall(iAmp,idps).RMS(iPOIs).info};
                RMS_DB_info(iAmp,idps,iPOIs) = {RMSall(iAmp,idps).RMS_DB(iPOIs).info};
                if iPOIs<=liPOIs
                %Decompose Out reconsSignal
                rSVR_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVR{iPOIs};
                rSVRDOI_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVRDOI{iPOIs};
                trS_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.timereconS{iPOIs};
                
                rSVR_DB_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out_DB.reconsSignalVR{iPOIs};
                rSVRDOI_DB_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out_DB.reconsSignalVRDOI{iPOIs};
                trS_DB_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out_DB.timereconS{iPOIs};
                end
            end
            % Decompose SNR structure
            SNR_DOI_noiseall_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOI_noiseall;
            SNR_DOI_noiseallptherm_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOI_noiseallptherm;
            SNR_DOI_noiseosc_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOI_noiseosc;
            SNR_DOI_noisestat_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOI_noisestat;
            SNR_info(iAmp,idps,:) = SNRall(iAmp,idps).SNR.info(:);
            SNR_DB_DOI_noiseall_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR_DB.DOI_noiseall;
            SNR_DB_DOI_noiseallptherm_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR_DB.DOI_noiseallptherm;
            SNR_DB_DOI_noiseosc_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR_DB.DOI_noiseosc;
            SNR_DB_DOI_noisestat_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR_DB.DOI_noisestat;
            SNR_DB_info(iAmp,idps,:) = SNRall(iAmp,idps).SNR_DB.info(:);
            
            % Decompose Out
            VR_mat(iAmp,idps,1:liPOIs,:) = Outall(iAmp,idps).Out.VR;
            VRptherm_mat(iAmp,idps,1:liPOIs,:) = Outall(iAmp,idps).Out.VRptherm;
            VRDOI_mat(iAmp,idps,1:liPOIs,:) = Outall(iAmp,idps).Out.VRDOI;
            VR_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out.mVR_POI;
            VR_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out.mVR_PsO;
            VRDOI_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out.mVRDOI_POI;
            VRDOI_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out.mVRDOI_PsO;
            VRptherm_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out.mVRptherm_POI;
            VRptherm_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out.mVRptherm_PsO;
            
            VR_DB_mat(iAmp,idps,1:liPOIs,:) = Outall(iAmp,idps).Out_DB.VR;
            VRptherm_DB_mat(iAmp,idps,1:liPOIs,:) = Outall(iAmp,idps).Out_DB.VRptherm;
            VRDOI_DB_mat(iAmp,idps,1:liPOIs,:) = Outall(iAmp,idps).Out_DB.VRDOI;
            VR_DB_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out_DB.mVR_POI;
            VR_DB_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out_DB.mVR_PsO;
            VRDOI_DB_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out_DB.mVRDOI_POI;
            VRDOI_DB_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out_DB.mVRDOI_PsO;
            VRptherm_DB_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out_DB.mVRptherm_POI;
            VRptherm_DB_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out_DB.mVRptherm_PsO;
            
            
            %Decompose Out reconsSignal
            rSVR_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out.reconsSignalmVR_POI{1};
            rSVR_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out.reconsSignalmVR_PsO{1};
            rSVRDOI_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out.reconsSignalmVRDOI_POI{1};
            rSVRDOI_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out.reconsSignalmVRDOI_PsO{1};
            trS_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out.timereconS_POI{1};
            trS_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out.timereconS_PsO{1};
            
            rSVR_DB_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out_DB.reconsSignalmVR_POI{1};
            rSVR_DB_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out_DB.reconsSignalmVR_PsO{1};
            rSVRDOI_DB_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out_DB.reconsSignalmVRDOI_POI{1};
            rSVRDOI_DB_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out_DB.reconsSignalmVRDOI_PsO{1};
            trS_DB_mat(iAmp,idps,liPOIs+1,:) = Outall(iAmp,idps).Out_DB.timereconS_POI{1};
            trS_DB_mat(iAmp,idps,liPOIs+2,:) = Outall(iAmp,idps).Out_DB.timereconS_PsO{1};
            
            %Decompose TStOP
            
            StartSim_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP.startSim;
            Itgen_mat(iAmp,idps,1:length(TSTOPall(iAmp,idps).TSTOP.Itgen)) = TSTOPall(iAmp,idps).TSTOP.Itgen;
            VRgen_mat(iAmp,idps,1:length(TSTOPall(iAmp,idps).TSTOP.VRgen)) = TSTOPall(iAmp,idps).TSTOP.VRgen;
            endsim_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP.endsim;
            signalreconstruction_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP.signalreconstruction;
            
            StartSim_DB_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP_DB.startSim;
            Itgen_DB_mat(iAmp,idps,1:length(TSTOPall(iAmp,idps).TSTOP.Itgen)) = TSTOPall(iAmp,idps).TSTOP_DB.Itgen;
            VRgen_DB_mat(iAmp,idps,1:length(TSTOPall(iAmp,idps).TSTOP.VRgen)) = TSTOPall(iAmp,idps).TSTOP_DB.VRgen;
            endsim_DB_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP_DB.endsim;
            signalreconstruction_DB_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP_DB.signalreconstruction;
            
            Outall(iAmp,idps).Out.VR = [];
            Outall(iAmp,idps).Out.VRptherm = [];
            Outall(iAmp,idps).Out.VRDOI = [];
            Outall(iAmp,idps).Out.mVR_POI = [];
            Outall(iAmp,idps).Out.mVR_PsO = [];
            Outall(iAmp,idps).Out.mVRDOI_POI = [];
            Outall(iAmp,idps).Out.mVRDOI_PsO = [];
            Outall(iAmp,idps).Out.mVRphterm_POI = [];
            Outall(iAmp,idps).Out.mVRptherm_PsO = [];
            Outall(iAmp,idps).Out.reconsSignalVR = [];
            Outall(iAmp,idps).Out.reconsSignalmVR_POI = [];
            Outall(iAmp,idps).Out.reconsSignalmVR_PsO = [];
            Outall(iAmp,idps).Out.reconsSignalVRDOI = [];
            Outall(iAmp,idps).Out.reconsSignalmVRDOI_POI = [];
            Outall(iAmp,idps).Out.reconsSignalmVRDOI_PsO = [];
            Outall(iAmp,idps).Out.timereconS = [];
            Outall(iAmp,idps).Out.timereconS_POI = [];
            Outall(iAmp,idps).Out.timereconS_PsO = [];
            
            Outall(iAmp,idps).Out_DB.VR = [];
            Outall(iAmp,idps).Out_DB.VRptherm = [];
            Outall(iAmp,idps).Out_DB.VRDOI = [];
            Outall(iAmp,idps).Out_DB.mVR_POI = [];
            Outall(iAmp,idps).Out_DB.mVR_PsO = [];
            Outall(iAmp,idps).Out_DB.mVRDOI_POI = [];
            Outall(iAmp,idps).Out_DB.mVRDOI_PsO = [];
            Outall(iAmp,idps).Out_DB.mVRphterm_POI = [];
            Outall(iAmp,idps).Out_DB.mVRptherm_PsO = [];
            Outall(iAmp,idps).Out_DB.reconsSignalVR = [];
            Outall(iAmp,idps).Out_DB.reconsSignalmVR_POI = [];
            Outall(iAmp,idps).Out_DB.reconsSignalmVR_PsO = [];
            Outall(iAmp,idps).Out_DB.reconsSignalVRDOI = [];
            Outall(iAmp,idps).Out_DB.reconsSignalmVRDOI_POI = [];
            Outall(iAmp,idps).Out_DB.reconsSignalmVRDOI_PsO = [];
            Outall(iAmp,idps).Out_DB.timereconS = [];
            Outall(iAmp,idps).Out_DB.timereconS_POI = [];
            Outall(iAmp,idps).Out_DB.timereconS_PsO = [];
            
            

        end
    end
    disp('succes');
    % calc memory difference
    To = whos('Outall','RMSall','SNRall','TSTOPall');
    dataTo = sum([To(:).bytes]);
    Ta = whos('*_mat');
    dataTa = sum([Ta(:).bytes]);
    dToTa = dataTo-dataTa; % in matlab memory differences not evan a MB but when saved 13% reduction
    
    namestosave = horzcat({Ta(:).name},{'SNR_DB_info','SNR_info','RMS_DB_info','RMS_info','Input','Outall'});
    
    disp('');
    disp(['ended, start saving']);
    if HPC_flag || DEBUG
        idxdoublepoint = find(SettingsStr==':');
        if isempty(idxdoublepoint)
            save_str_int = SettingsStr;
        else
            save_str_int = SettingsStr([1:idxdoublepoint-1,idxdoublepoint+2:end]);
        end
        save(['Results_',save_str_int,'_',datestr(now,'mm-dd-yy_HHMM'),'.mat'],namestosave{:});
    end
else
disp('');
disp(['ended, start saving']);
if HPC_flag || DEBUG
    idxdoublepoint = find(SettingsStr==':');
    if isempty(idxdoublepoint)
        save_str_int = SettingsStr;
    else
        save_str_int = SettingsStr([1:idxdoublepoint-1,idxdoublepoint+2:end]);
    end
    save(['Results_',save_str_int,'_',datestr(now,'mm-dd-yy_HHMM'),'.mat'],'RMSall','SNRall','Outall','TSTOPall','Input','-v7.3');
end
end

end
