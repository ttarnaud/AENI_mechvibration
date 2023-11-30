% Create results figure 2 for paper
% results from ampus study (change in spatial constant of ultrasonic field)
close all; clear all; clc

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1231\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_Amp_fkApdpset_NCfixed_nos4l_v2_12-30-20_1825';
data_rr_Amp05Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1224\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_kAus_nos4l_v2_12-28-20_1901';
data_rr_kAus05Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220304\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_kAus_nos4l_th10_v2_03-03-22_1625';
data_rr_kAus10Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220307\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_Amp_fkApdpset_NCfixed_th10_nos4l_v2_03-04-22_1746';
data_rr_Amp10Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220304_aIfun\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_kAus_nos4l_th5_v2_03-08-22_1408';
data_rr_kAus05Atrain =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220310_aIfun\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_Amp_fkApdpset_NCfixed_th5_nos4l_v2_03-09-22_2022';
data_rr_Amp05Atrain =  load(fullfile(folder_rawresult,filename_rawresult));

Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI
Atrain = @(t) alphaTrainFun(t,'./Inputs/aTInput_220308.mat');
%%
close all
data_rr = data_rr_kAus05Afun;

[M,N] = size(data_rr.Outall);
Amps = nan(M,N);
kaus = nan(M,N);
for i =1:M
    for j =1:N
        kaus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
        Amps(i,j) = data_rr.Outall(i,j).Out.Param.Aus;
    end
end
kAus = 1./kaus(:,1).*1000;
Amps = unique(Amps)*1e6;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kaus)),'\n'],kAus)
fprintf(['Amps: ', repmat('%5.2e ', 1, length(Amps)),'\n'],Amps')
Totaldps = data_rr.Input.Totaldps;
fprintf(['totaldps: ', repmat('%5.2e ', 1, length(Totaldps)),'\n'],Totaldps')

figure(10)
idp = 9;
plotmethod ='norm';
myvals = 0;
ikaus = 4;
styles = {'-','--'};

% plot alpha 
Tsim = 0:0.1:data_rr.trS_mat(ikaus,idp,end,end)*1000;
plot(Tsim,Ifun(Tsim/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on

pois = [1,2];
kaus = [4,5];
colors = flare(7);
cm = colors(:,:);
colors = cm([1,4,6,3,5,2,7],:);

for ipoi=1:length(pois)
    for ikaus=1:length(kaus)
        yvals = squeeze(data_rr.rSVR_mat(kaus(ikaus),idp,pois(ipoi),:));
        ratio = min(yvals)/max(yvals);
        if abs(ratio)>1; flip = -1;else flip = 1; end
        if strcmpi(plotmethod,'flip')
            yvals = flip*yvals;
        elseif strcmpi(plotmethod,'norm')
            [~,idx_max] = max(abs(yvals));
            yvals = yvals/yvals(idx_max);
        else
            error('incorrect value for plotflag should be either flip or norm')
        end
        lPOI = sprintf('POI_{%i}',pois(ipoi));
        if pois(ipoi)==6
            lPOI = 'mPOI';
        elseif pois(ipoi) == 7
            lPOI = 'mPsO';
        end
        plot(squeeze(data_rr.trS_mat(kaus(ikaus),idp,pois(ipoi),:))*1000,yvals,'linestyle',styles{ikaus},'color',colors(pois(ipoi),:),...
            'linewidth',1,'HandleVisibility','off')
        myvals = max(max(yvals),myvals);
    end
end

for ipoi = 1:length(pois)
    plot(nan,nan,'color',colors(pois(ipoi),:),'linewidth',1,'DisplayName',sprintf('POI_{%i}',pois(ipoi)))
end
for ikaus = 1:length(kaus)
    plot(nan,nan,'color',[0.1,0.1,0.1],'linestyle',styles{ikaus},'linewidth',1,'DisplayName',sprintf('\\kappa = %0.2f',kAus(kaus(ikaus))))
end
hold off
xlim([0,25])
xlabel('time [ms]')
ylabel('norm. Signal [-]')
set(gca,'box','off')
title(num2str(Totaldps(idp)))
l = legend('show','box','off','NumColumns',1);
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
set(l,'position',[0.5390    0.5414    0.3816    0.3476])
%%
data_rr = data_rr_kAus05Atrain;

[M,N] = size(data_rr.Outall);
Amps = nan(M,N);
kaus = nan(M,N);
for i =1:M
    for j =1:N
        kaus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
        Amps(i,j) = data_rr.Outall(i,j).Out.Param.Aus;
    end
end
kAus = 1./kaus(:,1).*1000;
Amps = unique(Amps)*1e6;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kaus)),'\n'],kAus)
fprintf(['Amps: ', repmat('%5.2e ', 1, length(Amps)),'\n'],Amps')

Totaldps = data_rr.Input.Totaldps;
fprintf(['totaldps: ', repmat('%5.2e ', 1, length(Totaldps)),'\n'],Totaldps')

figure(11)
idp = 3;
plotmethod ='norm';
myvals = 0;
ikaus = 4;
styles = {'-','--'};

% plot alpha 
Tsim = 0:0.1:data_rr.trS_mat(ikaus,idp,end,end)*1000;
plot(Tsim,Atrain(Tsim/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on

pois = [1,2];
kaus = [4,5];
colors = flare(7);
cm = colors(:,:);
colors = cm([1,4,6,3,5,2,7],:);

for ipoi=1:length(pois)
    for ikaus=1:length(kaus)
        yvals = squeeze(data_rr.rSVR_mat(kaus(ikaus),idp,pois(ipoi),:));
        ratio = min(yvals)/max(yvals);
        if abs(ratio)>1; flip = -1;else flip = 1; end
        if strcmpi(plotmethod,'flip')
            yvals = flip*yvals;
        elseif strcmpi(plotmethod,'norm')
            [~,idx_max] = max(abs(yvals));
            yvals = yvals/yvals(idx_max);
        else
            error('incorrect value for plotflag should be either flip or norm')
        end
        lPOI = sprintf('POI_{%i}',pois(ipoi));
        if pois(ipoi)==6
            lPOI = 'mPOI';
        elseif pois(ipoi) == 7
            lPOI = 'mPsO';
        end
        plot(squeeze(data_rr.trS_mat(kaus(ikaus),idp,pois(ipoi),:))*1000,yvals,'linestyle',styles{ikaus},'color',colors(pois(ipoi),:),...
            'linewidth',1,'HandleVisibility','off')
        myvals = max(max(yvals),myvals);
    end
end

for ipoi = 1:length(pois)
    plot(nan,nan,'color',colors(pois(ipoi),:),'linewidth',1,'DisplayName',sprintf('POI_{%i}',pois(ipoi)))
end
for ikaus = 1:length(kaus)
    plot(nan,nan,'color',[0.1,0.1,0.1],'linestyle',styles{ikaus},'linewidth',1,'DisplayName',sprintf('\\kappa = %0.2f',kAus(kaus(ikaus))))
end
hold off
xlim([0,25])
xlabel('time [ms]')
ylabel('norm. Signal [-]')
set(gca,'box','off')
title(num2str(Totaldps(idp)))
l = legend('show','box','off','NumColumns',1);
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
set(l,'position',[0.5390    0.5414    0.3816    0.3476])

%% CORTEX
clear al
folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1231\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_Amp_fkApdpset_NCfixed_nos4l_v2_12-30-20_2141';
data_rr_Amp05Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1224\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_kAus_nos4l_v2_12-29-20_0110';
data_rr_kAus05Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220304\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_kAus_nos4l_th10_v2_03-03-22_1754';
data_rr_kAus10Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220307\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_Amp_fkApdpset_NCfixed_th10_nos4l_v2_03-04-22_1753';
data_rr_Amp10Afun =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220304_aIfun\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_kAus_nos4l_th5_v2_03-08-22_1411';
data_rr_kAus05Atrain =  load(fullfile(folder_rawresult,filename_rawresult));

folder_rawresult = 'D:\no backup\EEGUS\HPC_files\20220310_aIfun\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_Amp_fkApdpset_NCfixed_th5_nos4l_v2_03-09-22_2102';
data_rr_Amp05Atrain =  load(fullfile(folder_rawresult,filename_rawresult));

Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI
Atrain = @(t) alphaTrainFun(t,'./Inputs/aTInput_220308.mat');
%%

data_rr = data_rr_kAus05Afun;

[M,N] = size(data_rr.Outall);
Amps = nan(M,N);
kaus = nan(M,N);
for i =1:M
    for j =1:N
        kaus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
        Amps(i,j) = data_rr.Outall(i,j).Out.Param.Aus;
    end
end
kAus = 1./kaus(:,1).*1000;
Amps = unique(Amps)*1e6;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kaus)),'\n'],kAus)
fprintf(['Amps: ', repmat('%5.2e ', 1, length(Amps)),'\n'],Amps')
Totaldps = data_rr.Input.Totaldps;
fprintf(['totaldps: ', repmat('%5.2e ', 1, length(Totaldps)),'\n'],Totaldps')

figure(10)
idp = 9;
plotmethod ='norm';
myvals = 0;
ikaus = 4;
styles = {'-','--'};

% plot alpha 
Tsim = 0:0.1:data_rr.trS_mat(ikaus,idp,end,end)*1000;
plot(Tsim,Ifun(Tsim/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on

pois = [1,2];
kaus = [4,5];
colors = flare(7);
cm = colors(:,:);
colors = cm([1,4,6,3,5,2,7],:);

for ipoi=1:length(pois)
    for ikaus=1:length(kaus)
        yvals = squeeze(data_rr.rSVR_mat(kaus(ikaus),idp,pois(ipoi),:));
        ratio = min(yvals)/max(yvals);
        if abs(ratio)>1; flip = -1;else flip = 1; end
        if strcmpi(plotmethod,'flip')
            yvals = flip*yvals;
        elseif strcmpi(plotmethod,'norm')
            [~,idx_max] = max(abs(yvals));
            yvals = yvals/yvals(idx_max);
        else
            error('incorrect value for plotflag should be either flip or norm')
        end
        lPOI = sprintf('POI_{%i}',pois(ipoi));
        if pois(ipoi)==6
            lPOI = 'mPOI';
        elseif pois(ipoi) == 7
            lPOI = 'mPsO';
        end
        plot(squeeze(data_rr.trS_mat(kaus(ikaus),idp,pois(ipoi),:))*1000,yvals,'linestyle',styles{ikaus},'color',colors(pois(ipoi),:),...
            'linewidth',1,'HandleVisibility','off')
        myvals = max(max(yvals),myvals);
    end
end

for ipoi = 1:length(pois)
    plot(nan,nan,'color',colors(pois(ipoi),:),'linewidth',1,'DisplayName',sprintf('POI_{%i}',pois(ipoi)))
end
for ikaus = 1:length(kaus)
    plot(nan,nan,'color',[0.1,0.1,0.1],'linestyle',styles{ikaus},'linewidth',1,'DisplayName',sprintf('\\kappa = %0.2f',kAus(kaus(ikaus))))
end
hold off
xlim([0,25])
xlabel('time [ms]')
ylabel('norm. Signal [-]')
set(gca,'box','off')
title(num2str(Totaldps(idp)))
l = legend('show','box','off','NumColumns',1);
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
set(l,'position',[0.5390    0.5414    0.3816    0.3476])
%%
data_rr = data_rr_kAus05Atrain;

[M,N] = size(data_rr.Outall);
Amps = nan(M,N);
kaus = nan(M,N);
for i =1:M
    for j =1:N
        kaus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
        Amps(i,j) = data_rr.Outall(i,j).Out.Param.Aus;
    end
end
kAus = 1./kaus(:,1).*1000;
Amps = unique(Amps)*1e6;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kaus)),'\n'],kAus)
fprintf(['Amps: ', repmat('%5.2e ', 1, length(Amps)),'\n'],Amps')

Totaldps = data_rr.Input.Totaldps;
fprintf(['totaldps: ', repmat('%5.2e ', 1, length(Totaldps)),'\n'],Totaldps')

figure(11)
idp = 3;
plotmethod ='norm';
myvals = 0;
ikaus = 4;
styles = {'-','--'};

% plot alpha 
Tsim = 0:0.1:data_rr.trS_mat(ikaus,idp,end,end)*1000;
plot(Tsim,Atrain(Tsim/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on

pois = [1,2];
kaus = [4,5];
colors = flare(7);
cm = colors(:,:);
colors = cm([1,4,6,3,5,2,7],:);

for ipoi=1:length(pois)
    for ikaus=1:length(kaus)
        yvals = squeeze(data_rr.rSVR_mat(kaus(ikaus),idp,pois(ipoi),:));
        ratio = min(yvals)/max(yvals);
        if abs(ratio)>1; flip = -1;else flip = 1; end
        if strcmpi(plotmethod,'flip')
            yvals = flip*yvals;
        elseif strcmpi(plotmethod,'norm')
            [~,idx_max] = max(abs(yvals));
            yvals = yvals/yvals(idx_max);
        else
            error('incorrect value for plotflag should be either flip or norm')
        end
        lPOI = sprintf('POI_{%i}',pois(ipoi));
        if pois(ipoi)==6
            lPOI = 'mPOI';
        elseif pois(ipoi) == 7
            lPOI = 'mPsO';
        end
        plot(squeeze(data_rr.trS_mat(kaus(ikaus),idp,pois(ipoi),:))*1000,yvals,'linestyle',styles{ikaus},'color',colors(pois(ipoi),:),...
            'linewidth',1,'HandleVisibility','off')
        myvals = max(max(yvals),myvals);
    end
end

for ipoi = 1:length(pois)
    plot(nan,nan,'color',colors(pois(ipoi),:),'linewidth',1,'DisplayName',sprintf('POI_{%i}',pois(ipoi)))
end
for ikaus = 1:length(kaus)
    plot(nan,nan,'color',[0.1,0.1,0.1],'linestyle',styles{ikaus},'linewidth',1,'DisplayName',sprintf('\\kappa = %0.2f',kAus(kaus(ikaus))))
end
hold off
xlim([0,25])
xlabel('time [ms]')
ylabel('norm. Signal [-]')
set(gca,'box','off')
title(num2str(Totaldps(idp)))
l = legend('show','box','off','NumColumns',1);
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
set(l,'position',[0.5390    0.5414    0.3816    0.3476])

%% LOAD DATA
clear all
data_coll_AMP05Afun = load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_a2piCondBug_v2.mat');
data_coll_AMP05Atrain = load('./Results/SCs_789_cEG_Amp_nos4l_fkApdpset_a2piCondBug_th5_aIfun_v2.mat');
data_coll_AMP10Afun = load('./Results/SCs_789_cEG_Amp_nos4l_fkApdpset_a2piCondBug_th10_v2.mat');


dps1 = data_coll_AMP05Afun.Inputall(1).Input.Totaldps;
dps2 = data_coll_AMP05Atrain.Inputall(1).Input.Totaldps;
dps3 = data_coll_AMP10Afun.Inputall(1).Input.Totaldps;


dps = unique([dps1,dps2,dps3]);

RequiredAmp1 = data_coll_AMP05Afun.requiredAmpAll;
RequiredAmp2 = data_coll_AMP05Atrain.requiredAmpAll;
RequiredAmp3 = data_coll_AMP10Afun.requiredAmpAll;


M = max([size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1)]);
RequiredAmp = nan(M,length(dps),3);


nLayers = [3,5,10];
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
ROIs = {'Deep','Cortex'};
[nL_mat,ROI_mat,Models_mat] = ndgrid(nLayers,ROIs,Models);
nL_mat = nL_mat(:)'; ROI_mat = ROI_mat(:)'; Models_mat = Models_mat(:)';

[~,idx_dps,idx_dps2] = intersect(dps,dps1);
for i=1:length(data_coll_AMP05Afun.Inputall)
    Settings = data_coll_AMP05Afun.Inputall(i).Input.Settings;
    info1(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info1(i).Model = data_coll_AMP05Afun.Inputall(i).Input.Model;
    info1(i).ROI = data_coll_AMP05Afun.Inputall(i).Input.ROI;
    
    
    idx = info1(i).nLayers==nL_mat & strcmpi(info1(i).Model,Models_mat) &...
        strcmpi(info1(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,1) = RequiredAmp1(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps2);
for i=1:length(data_coll_AMP05Atrain.Inputall)
    Settings = data_coll_AMP05Atrain.Inputall(i).Input.Settings;
    info2(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info2(i).Model = data_coll_AMP05Atrain.Inputall(i).Input.Model;
    info2(i).ROI = data_coll_AMP05Atrain.Inputall(i).Input.ROI;
    
    idx = info2(i).nLayers==nL_mat & strcmpi(info2(i).Model,Models_mat) & strcmpi(info2(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,2) = RequiredAmp2(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps3);
for i=1:length(data_coll_AMP10Afun.Inputall)
    Settings = data_coll_AMP10Afun.Inputall(i).Input.Settings;
    info3(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info3(i).Model = data_coll_AMP10Afun.Inputall(i).Input.Model;
    info3(i).ROI = data_coll_AMP10Afun.Inputall(i).Input.ROI;
    
    idx = info3(i).nLayers==nL_mat & strcmpi(info3(i).Model,Models_mat) & strcmp(info3(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,3) = RequiredAmp3(i,idx_dps2);
    
end
%% Lineplot Required amp vs # dipoles
Colors = thermal(4);

close all
types = {'-','-.',':','--'};
Markers = {'o','x','d','s'};
ytarget = 0.05;
dps_nr = [7,8,9];
YVAL = [];
Rsqval = [];
alphaval_ytarget = [];
Model_rval = [];
PSDval_ytarget = [];
bs = [];
bs2 = [];
sets = [1,2,3];

for isubplot = 1:length(Models)
    figure()
    for iM = 1:length(Models_mat)
        for iset = 1:length(sets)        
            
            idx_c = iset;
            idx_lt = find(strcmpi(ROI_mat{iM},ROIs));
            idx_m = find(strcmpi(ROI_mat{iM},ROIs));
            ltype = [Markers{isubplot},types{idx_lt}];
            yval = squeeze(RequiredAmp(iM,dps_nr,iset));
            yval = yval*1e6; % convert to um
            
            flag = strcmpi(Models_mat{iM},Models{isubplot}) && nL_mat(iM) == 3;
            
            if length(yval(~isnan(yval)))>1 && flag

                yintm = log10(yval'/10^mean(log10(yval)));
                %yintm = log10(yval'/yval(1));
                YVAL = [YVAL;yintm'];
                Model_rval = [Model_rval;[isubplot,iM]];
                b = [ones(numel(yval),1),log10(dps(dps_nr)')]\yintm;
                bs2 = [bs2,b];
                Rval = 1-sum((yintm'-[ones(numel(yval),1),log10(dps(dps_nr)')]*b).^2)/sum((yintm'-mean(yintm)).^2);
                Rsqval = [Rsqval;Rval];
                
                % find x for y =100nm
                yintm = log10(yval');
                X0 = [ones(numel(yval),1),log10(dps(dps_nr)')];
                b = X0\yintm;
                bs = [bs,b];
                xval_ytarget = (log10(ytarget)-b(1))/b(2);
                PSDval = calcPSD1MHz(xval_ytarget);
                alphaval_ytarget = [alphaval_ytarget;xval_ytarget];
                PSDval_ytarget = [PSDval_ytarget;PSDval];

                plot(dps(dps_nr)',yval,ltype(1),'color',Colors(idx_c,:),'HandleVisibility','off');
                hold on
                plot(dps(dps_nr)',10.^([ones(numel(yval),1),log10(dps(dps_nr)')]*b),ltype(2:end),'color',Colors(idx_c,:),'HandleVisibility','off');
            end          
            
        end
    end
    mytitle = split(Models{isubplot},'_');
    if strcmpi(mytitle{1},'mouse')
        mytitle{2} = 'scalp';
    end
    title(sprintf('%s_{%s}',mytitle{:}))
    set(gca,{'box','xscale','yscale'},{'off','log','log'})

    %ylabel('A_{max} [um]')
    set(findall(gcf,'type','axes'),'fontsize',11)
    set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
        {'centimeters',[1,1,1],[1,3,4.5,5.1],'centimeters',[1+4.5,3+5.1],'Painters'})
    
end
for i=1:length(dps_nr)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',sprintf('dps = %0.0e',dps(dps_nr(i))));
end
for i=1:length(ROIs)
    plot(nan,nan,[types{i}],'color','k','DisplayName',ROIs{i});
end
%for i=1:length(nLayers)
%    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',nLayers(i)));
%end
l = legend('show','box','off','NumColumns',1,'fontsize',11);


%-----------------------------------------------------------------------
% -----------------linear regresion----------------------------------
%-----------------------------------------------------------------------
XVAL = repmat(log10(dps(dps_nr)),size(YVAL,1),1);
b = [ones(numel(XVAL),1),XVAL(:)]\YVAL(:);

Rsq = 1-sum((YVAL(:)-[ones(numel(XVAL),1),XVAL(:)]*b).^2)/sum((YVAL(:)-mean(YVAL(:))).^2);
Rsq_min = min(Rsqval);
Rsq_max = max(Rsqval);
figure
for im = 1:length(Markers)
    xvals = XVAL((im-1)*6+1:im*6,:);
    yvals = YVAL((im-1)*6+1:im*6,:);
scatter(xvals(:),yvals(:),Markers{im});
hold on
end

plot(log10(dps(dps_nr)),[ones(3,1),log10(dps(dps_nr))']*b,'k-')
hold off

%% BARGRAPHS AMP
cmap = thermal(4);
cm  = nan(3);
cm(1,:) = cmap(2,:);
cm(3,:) = cmap(3,:);
cmap = thermal(9);
cm(2,:) = cmap(7,:);



dps_nr = [7,8,9];
sets = [1,2,3];

xticklabels = {};
for imodel=1:length(Models)
    myxticklabel = split(Models{imodel},'_');
    if strcmpi(myxticklabel{1},'mouse')
        myxticklabel{2} = 'scalp';
    end
    xticklabels{imodel} = sprintf('%s_{%s}',myxticklabel{:});
end

for isubplot = 1:length(ROIs)
    figure()
    
    idx_m = find(strcmpi(ROI_mat,ROIs{isubplot}));
    
    bardata = [RequiredAmp(idx_m,dps_nr,2)./RequiredAmp(idx_m,dps_nr,1),...
        RequiredAmp(idx_m,dps_nr,3)./RequiredAmp(idx_m,dps_nr,1)];
    bardata(any(isnan(bardata), 2),:) = [];
    hp = bar(bardata);
    %https://nl.mathworks.com/matlabcentral/answers/1670649-bar-plot-with-a-hatched-fill-pattern?s_tid=srchtitle


    
    

    %ylabel('A_{max} [um]')
    set(findall(gcf,'type','axes'),'fontsize',11)
    set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
        {'centimeters',[1,1,1],[1,3,8.5,5.1],'centimeters',[1+8.5,3+5.1],'Painters'})
    
    for i=1:3
        %hatchfill2(hp(i),'fill','HatchAngle',0,'hatchcolor',cm(i,:));
        hp(i).FaceColor = cm(i,:);
    end
    
    for i = 4:6
        hatchfill2(hp(i),'single','HatchAngle',45,'hatchcolor',cm(mod(i-1,3)+1,:));
        hp(i).FaceColor = 'none';
    end
    %title(ROIs{isubplot})
    ylim([0,1.2])
    set(gca,'box','off')
    %set(gca,{'xticklabel','xticklabelrotation'},{xticklabels,45})

%https://nl.mathworks.com/matlabcentral/answers/478956-getting-hatchfill-to-properly-display-a-patch-legend
%Make the legend
Lentry = split(sprintf('%0.0e,',dps(dps_nr)),',')
Legend  = [Lentry(1:end-1)',{'alpha-train','RMSE = 10%'}];
[legend_h,object_h,plot_h,text_str] = legendflex(hp,Legend,'fontSize',10,'box','off');
%object_h is the handle for the lines, patches and text objects
%hatch the legend patches to match the patches 

hHatch 	= hatchfill2(object_h(10),'single','HatchAngle',45,'hatchColor',[0.3,0.3,0.3],'hatchdensity',20);
%I couldn't use hatchfill2 to set the FaceAlpha's so I did it manually
for ii=1:3
set(get(object_h(ii+5),'children'),'FaceColor',cm(ii,:));  
end
set(get(object_h(4+5),'children'),'FaceColor',[0.5,0.5,0.5]);
ylabel("A'_{max} [-]")
end
%% kAUS figures
clear all
data_coll_kAus05Afun = load('./Results/SCs_all_cEG_kAus_nos4l_a2piCondBug_v2.mat');
data_coll_kAus05Atrain = load('./Results/SCs_789_cEG_kAus_nos4l_a2piCondBug_th5_aIfun_v2.mat');
data_coll_kAus10Afun = load('./Results/SCs_789_cEG_kAus_nos4l_a2piCondBug_th10_v2.mat');


dps1 = data_coll_kAus05Afun.Inputall(1).Input.Totaldps;
dps2 = data_coll_kAus05Atrain.Inputall(1).Input.Totaldps;
dps3 = data_coll_kAus10Afun.Inputall(1).Input.Totaldps;


dps = unique([dps1,dps2,dps3]);

RequiredAmp1 = 1000./data_coll_kAus05Afun.requiredAmpAll;
RequiredAmp2 = 1000./data_coll_kAus05Atrain.requiredAmpAll;
RequiredAmp3 = 1000./data_coll_kAus10Afun.requiredAmpAll;


M = max([size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1)]);
RequiredAmp = nan(M,length(dps),3);


nLayers = [3,5,10];
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
ROIs = {'Deep','Cortex'};
[nL_mat,ROI_mat,Models_mat] = ndgrid(nLayers,ROIs,Models);
nL_mat = nL_mat(:)'; ROI_mat = ROI_mat(:)'; Models_mat = Models_mat(:)';

[~,idx_dps,idx_dps2] = intersect(dps,dps1);
for i=1:length(data_coll_kAus05Afun.Inputall)
    Settings = data_coll_kAus05Afun.Inputall(i).Input.Settings;
    info1(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info1(i).Model = data_coll_kAus05Afun.Inputall(i).Input.Model;
    info1(i).ROI = data_coll_kAus05Afun.Inputall(i).Input.ROI;
    
    
    idx = info1(i).nLayers==nL_mat & strcmpi(info1(i).Model,Models_mat) &...
        strcmpi(info1(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,1) = RequiredAmp1(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps2);
for i=1:length(data_coll_kAus05Atrain.Inputall)
    Settings = data_coll_kAus05Atrain.Inputall(i).Input.Settings;
    info2(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info2(i).Model = data_coll_kAus05Atrain.Inputall(i).Input.Model;
    info2(i).ROI = data_coll_kAus05Atrain.Inputall(i).Input.ROI;
    
    idx = info2(i).nLayers==nL_mat & strcmpi(info2(i).Model,Models_mat) & strcmpi(info2(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,2) = RequiredAmp2(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps3);
for i=1:length(data_coll_kAus10Afun.Inputall)
    Settings = data_coll_kAus10Afun.Inputall(i).Input.Settings;
    info3(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info3(i).Model = data_coll_kAus10Afun.Inputall(i).Input.Model;
    info3(i).ROI = data_coll_kAus10Afun.Inputall(i).Input.ROI;
    
    idx = info3(i).nLayers==nL_mat & strcmpi(info3(i).Model,Models_mat) & strcmp(info3(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,3) = RequiredAmp3(i,idx_dps2);
    
end
%% Lineplot Required amp vs # dipoles
Colors = thermal(4);

close all
types = {'-','-.',':','--'};
Markers = {'o','x','d','s'};
ytarget = 0.05;
dps_nr = [7,8,9];
YVAL = [];
Rsqval = [];
alphaval_ytarget = [];
Model_rval = [];
PSDval_ytarget = [];
bs = [];
bs2 = [];
sets = [1,2,3];

for isubplot = 1:length(Models)
    figure()
    for iM = 1:length(Models_mat)
        for iset = 1:length(sets)        
            
            idx_c = iset;
            idx_lt = find(strcmpi(ROI_mat{iM},ROIs));
            idx_m = find(strcmpi(ROI_mat{iM},ROIs));
            ltype = [Markers{isubplot},types{idx_lt}];
            yval = squeeze(RequiredAmp(iM,dps_nr,iset));
            
            flag = strcmpi(Models_mat{iM},Models{isubplot}) && nL_mat(iM) == 3;
            
            if length(yval(~isnan(yval)))>1 && flag

                yintm = log10(yval'/10^mean(log10(yval)));
                %yintm = log10(yval'/yval(1));
                YVAL = [YVAL;yintm'];
                Model_rval = [Model_rval;[isubplot,iM]];
                b = [ones(numel(yval),1),log10(dps(dps_nr)')]\yintm;
                bs2 = [bs2,b];
                Rval = 1-sum((yintm'-[ones(numel(yval),1),log10(dps(dps_nr)')]*b).^2)/sum((yintm'-mean(yintm)).^2);
                Rsqval = [Rsqval;Rval];
                
                % find x for y =100nm
                yintm = log10(yval');
                X0 = [ones(numel(yval),1),log10(dps(dps_nr)')];
                b = X0\yintm;
                bs = [bs,b];
                xval_ytarget = (log10(ytarget)-b(1))/b(2);
                PSDval = calcPSD1MHz(xval_ytarget);
                alphaval_ytarget = [alphaval_ytarget;xval_ytarget];
                PSDval_ytarget = [PSDval_ytarget;PSDval];

                plot(dps(dps_nr)',yval,ltype(1),'color',Colors(idx_c,:),'HandleVisibility','off');
                hold on
                plot(dps(dps_nr)',10.^([ones(numel(yval),1),log10(dps(dps_nr)')]*b),ltype(2:end),'color',Colors(idx_c,:),'HandleVisibility','off');
            end          
            
        end
    end
    mytitle = split(Models{isubplot},'_');
    if strcmpi(mytitle{1},'mouse')
        mytitle{2} = 'scalp';
    end
    title(sprintf('%s_{%s}',mytitle{:}))
    set(gca,{'box','xscale','yscale'},{'off','log','log'})

    %ylabel('A_{max} [um]')
    set(findall(gcf,'type','axes'),'fontsize',11)
    set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
        {'centimeters',[1,1,1],[1,3,4.5,5.1],'centimeters',[1+4.5,3+5.1],'Painters'})
    
end
for i=1:length(dps_nr)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',sprintf('dps = %0.0e',dps(dps_nr(i))));
end
for i=1:length(ROIs)
    plot(nan,nan,[types{i}],'color','k','DisplayName',ROIs{i});
end
%for i=1:length(nLayers)
%    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',nLayers(i)));
%end
l = legend('show','box','off','NumColumns',1,'fontsize',11);


%-----------------------------------------------------------------------
% -----------------linear regresion----------------------------------
%-----------------------------------------------------------------------
XVAL = repmat(log10(dps(dps_nr)),size(YVAL,1),1);
b = [ones(numel(XVAL),1),XVAL(:)]\YVAL(:);

Rsq = 1-sum((YVAL(:)-[ones(numel(XVAL),1),XVAL(:)]*b).^2)/sum((YVAL(:)-mean(YVAL(:))).^2);
Rsq_min = min(Rsqval);
Rsq_max = max(Rsqval);
figure
for im = 1:length(Markers)
    xvals = XVAL((im-1)*6+1:im*6,:);
    yvals = YVAL((im-1)*6+1:im*6,:);
scatter(xvals(:),yvals(:),Markers{im});
hold on
end

plot(log10(dps(dps_nr)),[ones(3,1),log10(dps(dps_nr))']*b,'k-')
hold off

%% BARGRAPHS kAUS
cmap = thermal(4);
cm  = nan(3);
cm(1,:) = cmap(2,:);
cm(3,:) = cmap(3,:);
cmap = thermal(9);
cm(2,:) = cmap(7,:);

RequiredAmp(22,9,2) = 0.01; %simulation was limited to 0.01 but by chance it used the lower boundary for this sim -> no iterative proces done so value inaccurate therefore fix do opposed limits
dps_nr = [7,8,9];
sets = [1,2,3];

xticklabels = {};
for imodel=1:length(Models)
    myxticklabel = split(Models{imodel},'_');
    if strcmpi(myxticklabel{1},'mouse')
        myxticklabel{2} = 'scalp';
    end
    xticklabels{imodel} = sprintf('%s_{%s}',myxticklabel{:});
end

for isubplot = 1:length(ROIs)
    figure()
    
    idx_m = find(strcmpi(ROI_mat,ROIs{isubplot}));
    
    bardata = [RequiredAmp(idx_m,dps_nr,2)./RequiredAmp(idx_m,dps_nr,1),...
        RequiredAmp(idx_m,dps_nr,3)./RequiredAmp(idx_m,dps_nr,1)];
    bardata(any(isnan(bardata), 2),:) = [];
    hp = bar(bardata);
    %https://nl.mathworks.com/matlabcentral/answers/1670649-bar-plot-with-a-hatched-fill-pattern?s_tid=srchtitle


    
    

    %ylabel('A_{max} [um]')
    set(findall(gcf,'type','axes'),'fontsize',11)
    set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
        {'centimeters',[1,1,1],[1,3,8.5,5.1],'centimeters',[1+8.5,3+5.1],'Painters'})
    
    for i=1:3
        %hatchfill2(hp(i),'fill','HatchAngle',0,'hatchcolor',cm(i,:));
            hp(i).FaceColor = 'flat';
            hp(i).CData = cm(i,:).*ones(size(hp(i).CData(:,1),2),1);
    end
    
    for i = 4:6
        hatchfill2(hp(i),'single','HatchAngle',45,'hatchcolor',cm(mod(i-1,3)+1,:));
        hp(i).FaceColor = 'none';
    end
    %title(ROIs{isubplot})
    ylim([0,2.5])
    set(gca,'box','off')
    %set(gca,{'xticklabel','xticklabelrotation'},{xticklabels,45})

%https://nl.mathworks.com/matlabcentral/answers/478956-getting-hatchfill-to-properly-display-a-patch-legend
%Make the legend
Lentry = split(sprintf('%0.0e,',dps(dps_nr)),',')
Legend  = [Lentry(1:end-1)',{'alpha-train','RMSE = 10%'}];
[legend_h,object_h,plot_h,text_str] = legendflex(hp,Legend,'fontSize',10,'box','off');
%object_h is the handle for the lines, patches and text objects
%hatch the legend patches to match the patches 

hHatch 	= hatchfill2(object_h(10),'single','HatchAngle',45,'hatchColor',[0.3,0.3,0.3],'hatchdensity',20);
%I couldn't use hatchfill2 to set the FaceAlpha's so I did it manually
for ii=1:3
set(get(object_h(ii+5),'children'),'FaceColor',cm(ii,:));  
end
set(get(object_h(4+5),'children'),'FaceColor',[0.5,0.5,0.5]);
ylabel("\kappa' [-]")
end
%%
function PSD = calcPSD1MHz(alpha2)
%noise not included
wn = 14;
theta = 0.6;         %max is reached between wn and theta*wn
w2 = 10^3;                      %breakpoint higher powerlaw
alpha1 = 2;
%noise2 =@(x) 10.^(randn(1,length(x)));  %different noise methods
%were tested see script An_spec
noise3 = @(x) lognrnd(0,0.5,1,length(x));
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
fun = @(x) y(x).*noise3(x);
PSD = y(1e6);
end
