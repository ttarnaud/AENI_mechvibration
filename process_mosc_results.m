%processing the results of full pressure fields on EM simulations
% processing results prior march 03 2020 will need to be debugged (less
% info was saved)
clear all
close all
addpath(genpath('./Functions'));
Debug_flag = 1;
PSM_flag = 1;
setnrtoplot = [];%[4]%,4,9,10,15,16,21,22]
collect_ptherm_flag = 0;
nos4l_flag = 1;
structbased = 1;
date_sel = datetime(2020,12,10,0,0,0);
% info was not saved
if datenum(date_sel)==datenum(datetime(2020,3,3,0,0,0))
%     switch lower(ROI)
%     case 'deep'
%         adj_drSphere = 0;
%         if ratioAmpflag
%             Input_dppos = [0,0,0];
%             Input_OrienDipole = [0,1,0];
%         end
%     case 'cortex'
%         adj_drSphere = 1;
%         if ratioAmpflag
%             Input_dppos = [0,RBrain-dRBrain,0];
%             Input_OrienDipole = [0,1,0];
%         end
%     otherwise
%         error('false input')
%     end
%     posDp = Input_dppos;
%     OrienDipole = Input_OrienDipole;
%     Dirus = OrienDipole;
    Dirus_way = 'radial';
    Phaseus = 0;
    Dirus = [0,1,0];
end
    
try
    s = load('randomseed.mat');
    s = s.s;
catch
    warning('new random seed')
    s = rng;
end
if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
    idxsplit = 2;
else
    idxsplit = 3;
end

%Old_dir = cd('D:\no backup\EEGUS\ResultsHPC');
%Old_dir = cd('D:\Users\rschoeters\Documents\Imec USEEG\Matlab\Results\Results_local');
%Old_dir = cd('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Results\sOSC_NCFixed\Fixed3\Crop');
%Old_dir = cd('D:\Users\rschoeters\Documents\Imec USEEG\HPC_files\1231\fixed2\Crop');
%Old_dir = cd('D:\Users\rschoeters\Documents\Imec USEEG\HPC_files\1231\Crop');
%Old_dir = cd('D:\no backup\EEGUS\HPC_files\1231\Crop');
%lall = dir('*.mat');
%cd('D:\no backup\EEGUS\HPC_files\20220228\Crop');
%lall2 = dir('*.mat');
%lall = cat(1,lall,lall2)

Old_dir = cd('D:\no backup\EEGUS\HPC_files\20220307\Crop');
lall = dir('*.mat');

% select wanted set of results based on date size and nos4l flag
% assign number of settingset to list
cd(Old_dir);
l = lall([lall.datenum]>datenum(date_sel) & [lall.bytes]>1e5); % only results after 14/02 and large enough file (to small probably error)
idx_nos4l = contains({l.name},'nos4l');
idx_sel= logical(abs((1-nos4l_flag)-idx_nos4l));   % sel indices based on nos4l flag
l = l(idx_sel);
splitnames = regexp({l(:).name},'_','split');
setnumbers = cellfun(@(x) regexp(x{idxsplit},'Settingset','split'),splitnames,'uniformOutput',false);
setnumbers = cell2mat(cellfun(@(x) str2double(x{2}),setnumbers,'uniformoutput',false));
[sns,idx_s] = sort(setnumbers);
l = l(idx_s);
l(1).nr = sns(1);
sns_cell = num2cell(sns);
[l(:).nr] = deal(sns_cell{:});

ll = length(l);
requiredAmpAll = [];
requiredAmpAllLow = [];
for i =1:ll
    ftl = fullfile(l(i).folder,l(i).name);
    fprintf('\nloading: %s\n',ftl)
    load(ftl);
    if exist('Decomp_struct_flag','var')
        if Decomp_struct_flag
            if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
                [RMSall,SNRall,Outall,TSTOPall] = recomp_struct_pr2009(RMS_mat,RMSDOI_mat,Q2_mat,Q2DOI_mat,RMSptherm_mat,...
                    RMSDOIptherm_mat,Q2ptherm_mat,Q2DOIptherm_mat,RMS_info,rSVR_mat,rSVRDOI_mat,trS_mat,...
                    SNR_DOI_noiseall_mat,SNR_DOI_noiseallptherm_mat,SNR_DOI_noiseosc_mat,SNR_DOI_noisestat_mat,...
                    SNR_info,VR_mat,VRptherm_mat,VRDOI_mat,VRnoise_stat_mat,StartSim_mat,Itgen_mat,...
                    VRgen_mat,endsim_mat,signalreconstruction_mat,POIs,Input,Outall,nos4l_flag);
                clearvars('RMS_mat','RMSDOI_mat','Q2_mat','Q2DOI_mat','RMSptherm_mat',...
                    'RMSDOIptherm_mat','Q2ptherm_mat','Q2DOIptherm_mat','RMS_info','rSVR_mat','rSVRDOI_mat','trS_mat',...
                    'SNR_DOI_noiseall_mat','SNR_DOI_noiseallptherm_mat','SNR_DOI_noiseosc_mat','SNR_DOI_noisestat_mat',...
                    'SNR_info','VR_mat','VRptherm_mat','VRDOI_mat','VRnoise_stat_mat','StartSim_mat','Itgen_mat',...
                    'VRgen_mat','endsim_mat','signalreconstruction_mat')
            else
                [RMSall,SNRall,Outall,TSTOPall] = recomp_struct(RMS_mat,RMSptherm_mat,RMSDOI_mat,RMSDOIvr_mat,RMSDOIvrptherm_mat,...
                    RMSDOIvrpOSCvr_mat,RMSDOIvrpstatnoise_mat,RMS_info,...
                    Q2_mat,Q2ptherm_mat,Q2DOI_mat,Q2DOIvr_mat,Q2DOIvrptherm_mat,Q2DOIvrpOSCvr_mat,Q2DOIvrpstatnoise_mat,...
                    rSVR_mat,rSVRDOIvr_mat,rSVRDOIptherm_mat,rSVRdvpov_mat,rSVRdvpsn_mat,rSVRDOI_mat,trS_mat,...
                    SNR_DOIvr_noiseall_mat,SNR_DOIvr_noiseallptherm_mat,SNR_DOIvr_noiseoscvr_mat,SNR_DOIvr_noisestat_mat,...
                    SNR_info,VR_mat,VRptherm_mat,VRDOI_mat,VRDOIstat_mat,VROSC_mat,VROSCstat_mat,VRstatnoise_mat,StartSim_mat,Itgen_mat,...
                    VRgen_mat,endsim_mat,signalreconstruction_mat,POIs,Input,Outall,nos4l_flag);
                clearvars('RMS_mat','RMSptherm_mat','RMSDOI_mat','RMSDOIvr_mat','RMSDOIvrptherm_mat',...
                    'RMSDOIvrpOSCvr_mat','RMSDOIvrpstatnoise_mat','RMS_info',...
                    'Q2_mat','Q2ptherm_mat','Q2DOI_mat','Q2DOIvr_mat','Q2DOIvrptherm_mat','Q2DOIvrpOSCvr_mat','Q2DOIvrpstatnoise_mat',...
                    'rSVR_mat','rSVRDOIvr_mat','rSVRDOIptherm_mat','rSVRdvpov_mat','rSVRdvpsn_mat','rSVRDOI_mat','trS_mat',...
                    'SNR_DOIvr_noiseall_mat','SNR_DOIvr_noiseallptherm_mat','SNR_DOIvr_noiseoscvr_mat','SNR_DOIvr_noisestat_mat',...
                    'VR_mat','VRptherm_mat','VRDOI_mat','VRDOIstat_mat','VROSC_mat','VROSCstat_mat','VRstatnoise_mat','StartSim_mat','Itgen_mat',...
                    'VRgen_mat','endsim_mat','signalreconstruction_mat')
                
            end
        end
    end
    
    if isfield(Input,'Dirus_way')
        Dirus_way = Input.Dirus_way;
    end
    if isfield(Input,'Dirus')
        Dirus = Input.Dirus;
    end
    if isfield(Input,'Phaseus')
        Phaseus = Input.Phaseus;
    end
    idx_sel = 1;
    for iAmp = 1:size(RMSall,1)
        yrms = RMSall(iAmp,idx_sel).RMS;
        [miniRMS,idx_miRMS] = min([yrms(:).RMS]);
        RMSMAT(i,iAmp) = miniRMS(1);
        bestPOIs{i,iAmp} = yrms(idx_miRMS).info;
        ysnr = SNRall(iAmp,end).SNR.DOIvr_NoiseAll(2,:);
        [maxiSNR,idx_miSNR] = max(ysnr);
        SNRMAT(i,iAmp) =  maxiSNR;
        bestPOIsSNR{i,iAmp} = SNRall(iAmp,idx_sel).SNR.info{idx_miSNR};
        if collect_ptherm_flag
            for inoiseamp=1:3
                if nos4l_flag
                    
                    
                    if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
                        [RMSptherm,SNRptherm] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.VR,Outall(iAmp,idx_sel).Out.VRDOI,s,inoiseamp,Input,...
                            nos4l_flag,[],Dirus_way,Outall(iAmp,idx_sel).Out.CSource,Outall(iAmp,idx_sel).Out.CSink,Input.Aus_way,...
                            Input.Amp(iAmp),Input.kAus(iAmp),Phaseus,Dirus);
                        [RMSpthermPOI,SNRpthermPOI] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.mVR_POI,Outall(iAmp,idx_sel).Out.mVRDOI_POI,s,inoiseamp,Input,...
                            nos4l_flag,[],Dirus_way,Outall(iAmp,idx_sel).Out.CSource,Outall(iAmp,idx_sel).Out.CSink,Input.Aus_way,...
                            Input.Amp(iAmp),Input.kAus(iAmp),Phaseus,Dirus);
                        [RMSpthermPsO,SNRpthermPsO] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.mVR_PsO,Outall(iAmp,idx_sel).Out.mVRDOI_PsO,s,inoiseamp,Input,...
                            nos4l_flag,[],Dirus_way,Outall(iAmp,idx_sel).Out.CSource,Outall(iAmp,idx_sel).Out.CSink,Input.Aus_way,...
                            Input.Amp(iAmp),Input.kAus(iAmp),Phaseus,Dirus);
                    else
                        [RMSptherm,SNRptherm] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.VR,Outall(iAmp,idx_sel).Out.VRDOIvr,s,inoiseamp,Input,...
                            nos4l_flag,Outall(iAmp,idx_sel).Out.Param,Dirus_way,Outall(iAmp,idx_sel).Out.CSource,Outall(iAmp,idx_sel).Out.CSink,Input.Aus_way,...
                            Input.Amp(iAmp),Input.kAus(iAmp),Phaseus,Dirus);
                    end
                    
                    
                    
                else
                    
                    if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
                        [RMSptherm,SNRptherm] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.VR,Outall(iAmp,idx_sel).Out.VRDOI,s,inoiseamp,Input,...
                            nos4l_flag,[],Outall(iAmp,idx_sel).Out.vampx,Outall(iAmp,idx_sel).Out.vampy,Outall(iAmp,idx_sel).Out.vampz,...
                            Outall(iAmp,idx_sel).Out.phasex,Outall(iAmp,idx_sel).Out.vphasey,Outall(iAmp,idx_sel).Out.vphasez);
                        [RMSpthermPOI,SNRpthermPOI] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.mVR_POI,Outall(iAmp,idx_sel).Out.mVRDOI_POI,s,inoiseamp,Input,...
                            nos4l_flag,[],Outall(iAmp,idx_sel).Out.vampx,Outall(iAmp,idx_sel).Out.vampy,Outall(iAmp,idx_sel).Out.vampz,...
                            Outall(iAmp,idx_sel).Out.phasex,Outall(iAmp,idx_sel).Out.vphasey,Outall(iAmp,idx_sel).Out.vphasez);
                        [RMSpthermPsO,SNRpthermPsO] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.mVR_PsO,Outall(iAmp,idx_sel).Out.mVRDOI_PsO,s,inoiseamp,Input,...
                            nos4l_flag,[],Outall(iAmp,idx_sel).Out.vampx,Outall(iAmp,idx_sel).Out.vampy,Outall(iAmp,idx_sel).Out.vampz,...
                            Outall(iAmp,idx_sel).Out.phasex,Outall(iAmp,idx_sel).Out.vphasey,Outall(iAmp,idx_sel).Out.vphasez);
                    else
                        [RMSptherm,SNRptherm] = ...
                            getRMSaSNRptherm(Outall(iAmp,idx_sel).Out.VR,Outall(iAmp,idx_sel).Out.VRDOIvr,s,inoiseamp,Input,...
                            nos4l_flag,Outall(iAmp,idx_sel).Out.Param,Outall(iAmp,idx_sel).Out.vampx,Outall(iAmp,idx_sel).Out.vampy,Outall(iAmp,idx_sel).Out.vampz,...
                            Outall(iAmp,idx_sel).Out.phasex,Outall(iAmp,idx_sel).Out.vphasey,Outall(iAmp,idx_sel).Out.vphasez);
                        
                    end
                end
                if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
                    RMSptherm = [RMSptherm,RMSpthermPOI,RMSpthermPsO];
                    SNRptherm = [SNRptherm,SNRpthermPOI,SNRpthermPsO];
                end
                [miniRMS,idx_miRMS] = min(RMSptherm);
                RMSpthermMAT(i,iAmp,inoiseamp) = miniRMS(1);
                bestPOIsptherm{i,iAmp,inoiseamp} = yrms(idx_miRMS).info;
            
            [maxiSNR,idx_miSNR] = max(SNRptherm);
            SNRMpthermAT(i,iAmp,inoiseamp) =  maxiSNR;
            bestPOIspthermSNR{i,iAmp,inoiseamp} = SNRall(iAmp,idx_sel).SNR.info{idx_miSNR};
            
            
            end
        end
    end
    
    
    
    if any(l(i).nr==setnrtoplot)
        if ~isfield(Input,'Dirus_way')
            Input.Dirus_way = Dirus_way;
        end
        if ~isfield(Input,'Dirus')
            Input.Dirus = Dirus;
        end
        if ~isfield(Input,'Phaseus')
            try
            Input.Phaseus = Phaseus;
            catch
                try
                Input.Phaseus = Outall(1,1).Out.Param.Phaseus;
                catch
                    Input.Phaseus = 0;
                end
            end
                
            
        end
        if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
            createfig(Input,Outall,RMSall,SNRall,POIs,Debug_flag,PSM_flag,s,nos4l_flag,[],date_sel)
        else
            createfig(Input,Outall,RMSall,SNRall,POIs,Debug_flag,PSM_flag,s,nos4l_flag,Outall(end,idx_sel).Out.Param,date_sel)
        end
    end
    requiredAmpAll = [requiredAmpAll;requiredAmp];
    requiredAmpAllLow = [requiredAmpAllLow;requiredAmplower];
    Inputall(i).Input = Input; 
    RMSALL(i).RMS = RMSall;
end
disp('finished')
%save('./Results/SCs_789_cEG_Amp_nos4l_fkApdpset_a2piCondBug_th10_v2.mat','Inputall','RMSALL','requiredAmpAll','requiredAmpAllLow')
%%
% close all
% if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
%             createfig(Input,Outall,RMSall,SNRall,POIs,Debug_flag,PSM_flag,s,nos4l_flag,[],date_sel)
%         else
%             createfig(Input,Outall,RMSall,SNRall,POIs,Debug_flag,PSM_flag,s,nos4l_flag,Outall(end,idx_sel).Out.Param,date_sel)
%         end

%% load Strenght cruves and plot of simulations Input_cEG_nos4l_v1_1-24
load('./Results/SCs_all_cEG_kAus_nos4l_a2piCondBug_v2.mat')
load('./Results/SCs_789_cEG_kAus_nos4l_a2piCondBug_th10_v2.mat')
%load('./Results/SCs_789_cEG_kAus_nos4l_a2piCondBug_th5_aIfun_v2.mat')


Totaldps = Inputall(1).Input.Totaldps;
Amp = Inputall(1).Input.Amp;

types = {'-','--',':','-.'};
Markers = {'o','x','d','s'};

iM = size(requiredAmpAll,1);
Models = {};
ROIs = {};
nLayers = [];
for i=1:iM
    Models = horzcat(Models,Inputall(i).Input.Model);
    ROIs = horzcat(ROIs,Inputall(i).Input.ROI);
    idx = find(strcmpi(Inputall(i).Input.Settings,'nLayers'));
    nLayers = horzcat(nLayers,Inputall(i).Input.Settings{idx+1});
end

uModels = unique(Models);
uROI = unique(ROIs);
unLayers = unique(nLayers);

Colors = lines(length(uModels));



figure
for i=1:iM
    idx_c = find(strcmpi(Models{i},uModels));
    idx_lt = find(strcmpi(ROIs{i},uROI));
    idx_m = find(nLayers(i)==unLayers);
    ltype = [Markers{idx_m},types{idx_lt}];
    yval = requiredAmpAll(i,:);
    yval = 1./yval*100;
    p = plot(Totaldps,yval,ltype,'color',Colors(idx_c,:),'markersize',10)
    if i==1; hold on; end
    
    set(gca,{'xscale','yscale'},{'log','log'})
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    
end

for i=1:length(uModels)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[uModels{i}(1:5),'_{',uModels{i}(7:end),'}']);
end
for i=1:length(uROI)
    plot(nan,nan,types{i},'color','k','DisplayName',uROI{i});
end
for i=1:length(unLayers)
    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',unLayers(i)));
end
hold off
legend('show','location','best')

s1_flag = 1;
s2_flag = 1;
figure
for i=1:iM
    idx_c = find(strcmpi(Models{i},uModels));
    idx_lt = find(strcmpi(ROIs{i},uROI));
    idx_m = find(nLayers(i)==unLayers);
    ltype = [Markers{idx_m},'-'];
    yval = requiredAmpAll(i,:);
    yval = 1./yval*100;
    subplot(1,2,idx_lt)
    p = plot(Totaldps,yval,ltype,'color',Colors(idx_c,:),'markersize',10)
    if s1_flag && idx_lt==1
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
    ylabel('spatial constant k [cm]')
    end
    if s2_flag && idx_lt==2
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
    end
    
    set(gca,{'xscale','yscale'},{'log','log'})
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    
end
for i=1:length(uModels)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[uModels{i}(1:5),'_{',uModels{i}(7:end),'}']);
end
for i=1:length(unLayers)
    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',unLayers(i)));
end

hold off
legend('show','location','best')
set(findall(gcf,'type','axes'),'fontsize',18)
%% load Strenght cruves and plot of simulations Input_cEG_nos4l_v1_1-24
%load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_v1.mat')
%load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_v2.mat')
%load('./Results/SCs_all_cEG_Amp_NCfixed3_sOSC_nos4l_fkApdpset_v2.mat')
%load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_a2piCondBug_v2.mat')
%load('./Results/SCs_789_cEG_Amp_nos4l_fkApdpset_a2piCondBug_th5_aIfun_v2.mat')
%load('./Results/SCs_789_cEG_Amp_nos4l_fkApdpset_a2piCondBug_th10_v2.mat')

Totaldps = Inputall(1).Input.Totaldps;
Amp = Inputall(1).Input.Amp;

types = {'-','--',':','-.'};
Markers = {'o','x','d','s'};

iM = size(requiredAmpAll,1);
Models = {};
ROIs = {};
nLayers = [];
for i=1:iM
    Models = horzcat(Models,Inputall(i).Input.Model);
    ROIs = horzcat(ROIs,Inputall(i).Input.ROI);
    idx = find(strcmpi(Inputall(i).Input.Settings,'nLayers'));
    nLayers = horzcat(nLayers,Inputall(i).Input.Settings{idx+1});
end

uModels = unique(Models);
uROI = unique(ROIs);
unLayers = unique(nLayers);

Colors = lines(length(uModels));



figure
for i=1:iM
    idx_c = find(strcmpi(Models{i},uModels));
    idx_lt = find(strcmpi(ROIs{i},uROI));
    idx_m = find(nLayers(i)==unLayers);
    ltype = [Markers{idx_m},types{idx_lt}];
    yval = requiredAmpAll(i,:);
    yval = yval*1e6;
    p = plot(Totaldps,yval,ltype,'color',Colors(idx_c,:),'markersize',10)
    if i==1; hold on; end
    
    set(gca,{'xscale','yscale'},{'log','log'})
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    
end

for i=1:length(uModels)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[uModels{i}(1:5),'_{',uModels{i}(7:end),'}']);
end
for i=1:length(uROI)
    plot(nan,nan,types{i},'color','k','DisplayName',uROI{i});
end
for i=1:length(unLayers)
    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',unLayers(i)));
end
hold off
legend('show','location','best')

s1_flag = 1;
s2_flag = 1;
figure
for i=1:iM
    idx_c = find(strcmpi(Models{i},uModels));
    idx_lt = find(strcmpi(ROIs{i},uROI));
    idx_m = find(nLayers(i)==unLayers);
    ltype = [Markers{idx_m},'-'];
    yval = requiredAmpAll(i,:);
    yval = yval*1e6;
    subplot(1,2,idx_lt)
    p = plot(Totaldps,yval,ltype,'color',Colors(idx_c,:),'markersize',10)
    if s1_flag && idx_lt==1
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
    ylabel('vibration Amp [um]')
    end
    if s2_flag && idx_lt==2
    hold on;
    title(ROIs{i})
    xlabel('nr. dipoles')
    end
    
    set(gca,{'xscale','yscale'},{'log','log'})
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    
end
for i=1:length(uModels)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[uModels{i}(1:5),'_{',uModels{i}(7:end),'}']);
end
for i=1:length(unLayers)
    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',unLayers(i)));
end

hold off
legend('show','location','best')
set(findall(gcf,'type','axes'),'fontsize',18)