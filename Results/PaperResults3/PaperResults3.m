% Create results figure 2 for paper
% results from ampus study (change in spatial constant of ultrasonic field)
close all; clear all; clc
folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1231\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_Amp_fkApdpset_NCfixed_nos4l_v2_12-30-20_1825';
data_rr =  load(fullfile(folder_rawresult,filename_rawresult));


Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI

[M,N] = size(data_rr.Outall);
Amps = nan(M,N);
kaus = nan(M,N);
for i =1:M
    for j =1:N
        kaus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
        Amps(i,j) = data_rr.Outall(i,j).Out.Param.Aus;
    end
end
kaus = 1./kaus(:,1).*1000;
Amps = unique(Amps)*1e6;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kaus)),'\n'],kaus')
fprintf(['Amps: ', repmat('%5.2e ', 1, length(Amps)),'\n'],Amps')

%%
close all
iampus=1; idp=5; ipoi = 2;
plotmethod = 'flip'%'flip'; %'norm'
for iampus=3:6
    yvals = squeeze(data_rr.rSVR_mat(iampus,idp,ipoi,:));
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
    
    plot(squeeze(data_rr.trS_mat(iampus,idp,ipoi,:))*1000,yvals)
    hold on
end
hold off
%%
figure(10)
idp = 5;
plotmethod ='flip';
myvals = 0;
iampus = 4;
styles = {'--','-',':'};



pois = [1,2];
ampus = [5,6];
colors = flare(7);
cm = colors(:,:);
colors = cm([1,4,6,3,5,2,7],:);

% plot alpha 
Tsim = 0:0.1:data_rr.trS_mat(iampus,idp,end,end)*1000;
for ipoi=1:length(pois)
    for iampus=1:length(ampus)
        yvals = squeeze(data_rr.rSVR_mat(ampus(iampus),idp,pois(ipoi),:));
        myvals = max(max(yvals),myvals);
    end
end

plot(Tsim,myvals*Ifun(Tsim/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on

for ipoi=1:length(pois)
    for iampus=1:length(ampus)
        yvals = squeeze(data_rr.rSVR_mat(ampus(iampus),idp,pois(ipoi),:));
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
        plot(squeeze(data_rr.trS_mat(ampus(iampus),idp,pois(ipoi),:))*1000,yvals,'linestyle',styles{iampus},'color',colors(pois(ipoi),:),...
            'linewidth',1,'HandleVisibility','off')
        myvals = max(max(yvals),myvals);
    end
end

for ipoi = 1:length(pois)
    plot(nan,nan,'color',colors(pois(ipoi),:),'linewidth',1,'DisplayName',sprintf('POI_{%i}',pois(ipoi)))
end
for iampus = 1:length(ampus)
    plot(nan,nan,'color',[0.1,0.1,0.1],'linestyle',styles{iampus},'linewidth',1,'DisplayName',sprintf('A_{max} = %0.0f',Amps(ampus(iampus))))
end
hold off
xlim([0,25])
xlabel('time [ms]')
ylabel('Signal [uV]')
set(gca,'box','off')
l = legend('show','box','off','NumColumns',1);
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
set(l,'position',[0.5390    0.5414    0.3816    0.3476])
%%

colors = [27,158,119
217,95,2
117,112,179
231,41,138
102,166,30
230,171,2
166,118,29
102,102,102]/255;
colors = flare(7);
cm = colors;
colors = cm([1,4,6,3,5,2,7],:);
ampus = [4,5,6];
figure()
bardata = squeeze(data_rr.RMS_mat(ampus,idp,1:7));
b = bar(bardata);



xlim = get(gca,'xlim');

set(gca,'box','off')

set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,8.5,5.1],'centimeters',[1+8.5,3+5.1],'Painters'})

set(gca,'xticklabel',arrayfun(@(x) sprintf('%0.0f',x),Amps(ampus),'UniformOutput',false))
xlabel('A_{max}[\mum]')
ylabel('RMSE')
hold on
for i =1:length(b)
    b(i).FaceColor = 'flat';
    b(i).CData = colors(i,:).*ones(size(b(i).CData(:,1),2),1);
    bh(i) = bar(nan,nan,'FaceColor',colors(i,:));
end
set(gca,'xlim',xlim);
set(gca,'ylim',[0,0.81])

%% second row of figures
% Create results figure 2 for paper
% results from ampus study (change in spatial constant of ultrasonic field)
clear all;
folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1231\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_Amp_fkApdpset_NCfixed_nos4l_v2_12-30-20_2141';
data_rr =  load(fullfile(folder_rawresult,filename_rawresult));

folder_collective = 'D:\users\rschoeters\Documents\Imec USEEG\Matlab\Results';
filename_collective = 'SCs_all_cEG_Amp_nos4l_fkApdpset_a2piCondBug_v2';
data_col = load(fullfile(folder_collective,filename_collective));
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI

[M,N] = size(data_rr.Outall);
Amps = nan(M,N);
kaus = nan(M,N);
for i =1:M
    for j =1:N
        kaus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
        Amps(i,j) = data_rr.Outall(i,j).Out.Param.Aus;
    end
end
kaus = 1./kaus(:,1).*1000;
Amps = unique(Amps)*1e6;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kaus)),'\n'],kaus')
fprintf(['Amps: ', repmat('%5.2e ', 1, length(Amps)),'\n'],Amps')

%%
figure()
idp = 5;
plotmethod ='flip';
myvals = 0;
iampus = 4;
styles = {'--','-',':'};



pois = [1,2];
ampus = [5,6];
colors = flare(7);
cm = colors(:,:);
colors = cm([1,4,6,3,5,2,7],:);

% plot alpha 
Tsim = 0:0.1:data_rr.trS_mat(iampus,idp,end,end)*1000;
for ipoi=1:length(pois)
    for iampus=1:length(ampus)
        yvals = squeeze(data_rr.rSVR_mat(ampus(iampus),idp,pois(ipoi),:));
        myvals = max(max(yvals),myvals);
    end
end
if strcmpi(plotmethod,'norm')
    myvals = 1;
end
plot(Tsim,myvals*Ifun(Tsim/1000),'color',[0.9,0.9,0.9],'linewidth',3,'Handlevisibility','Off')
hold on

for ipoi=1:length(pois)
    for iampus=1:length(ampus)
        yvals = squeeze(data_rr.rSVR_mat(ampus(iampus),idp,pois(ipoi),:));
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
        plot(squeeze(data_rr.trS_mat(ampus(iampus),idp,pois(ipoi),:))*1000,yvals,'linestyle',styles{iampus},'color',colors(pois(ipoi),:),...
            'linewidth',1,'HandleVisibility','off')
        myvals = max(max(yvals),myvals);
    end
end

for ipoi = 1:length(pois)
    plot(nan,nan,'color',colors(pois(ipoi),:),'linewidth',1,'DisplayName',sprintf('POI_{%i}',pois(ipoi)))
end
for iampus = 1:length(ampus)
    plot(nan,nan,'color',[0.1,0.1,0.1],'linestyle',styles{iampus},'linewidth',1,'DisplayName',sprintf('A_{max} = %0.0f',Amps(ampus(iampus))))
end
hold off
xlim([0,25])
xlabel('time [ms]')
ylabel('Signal [uV]')
set(gca,'box','off')
l = legend('show','box','off','NumColumns',1);
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
set(l,'position',[0.5390    0.5414    0.3816    0.3476])
%%
colors = [27,158,119
217,95,2
117,112,179
231,41,138
102,166,30
230,171,2
166,118,29
102,102,102]/255;
colors = flare(7);
cm = colors;
colors = cm([1,4,6,3,5,2,7],:);
ampus = [4,5,6];
figure()
bardata = squeeze(data_rr.RMS_mat(ampus,idp,1:7));
b = bar(bardata);



xlim = get(gca,'xlim');

set(gca,'box','off')

set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,8.5,5.1],'centimeters',[1+8.5,3+5.1],'Painters'})

set(gca,'xticklabel',arrayfun(@(x) sprintf('%0.0f',x),Amps(ampus),'UniformOutput',false))
xlabel('A_{max}[\mum]')
ylabel('RMSE')
hold on
for i =1:length(b)
    b(i).FaceColor = 'flat';
    b(i).CData = colors(i,:).*ones(size(b(i).CData(:,1),2),1);
    bh(i) = bar(nan,nan,'FaceColor',colors(i,:));
end
set(gca,'xlim',xlim);
set(gca,'ylim',[0,0.81])

%%
clear all
folder_collective = 'D:\users\rschoeters\Documents\Imec USEEG\Matlab\Results';
filename_collective = 'SCs_all_cEG_Amp_nos4l_fkApdpset_a2piCondBug_v2';
data_col = load(fullfile(folder_collective,filename_collective));


Nin = length(data_col.Inputall);

colors = [27,158,119
217,95,2
117,112,179
231,41,138]/255;

%colors = thermal(5);
%colors = colors(1:2:end,:)

Ylims = [min(data_col.requiredAmpAll(:)),max(data_col.requiredAmpAll(:))];

Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
labels = {'scalp_{h}','cortical_{h}','air_{h}','scalp_{m}'};
styles = {'x-','o--','*:'};
ROIS = {'deep','cortex'};
NLAYERS = [3,5,10];
colors = [27,158,119
217,95,2
117,112,179
231,41,138]/255;

Rsq1 = nan(length(ROIS),length(Models));
Rsq2 = nan(length(ROIS),length(Models));
fignr = 100;
lsts = {'-','-.'};

for i = 1:length(ROIS)
figure(33+i)
data_sortedn = repmat({[]},1,4);
xdata_sortedn = repmat({[]},1,4);
for iin = 1:Nin

    myin = data_col.Inputall(iin).Input;
    Totaldps = myin.Totaldps;
    Model = myin.Model;
    ROI = myin.ROI;
    nLayers = myin.Settings{find(strcmpi('nLayers',myin.Settings))+1};
    iclr = find(strcmpi(Models,Model));
    ist = find([3,5,10]==nLayers);
    if strcmpi(ROI,ROIS{i})
        y = data_col.requiredAmpAll(iin,:)*1e6;
        plot(Totaldps,y,styles{ist}(1),'color',colors(iclr,:))
        hold on
        data_sortedn{iclr} = [data_sortedn{iclr},y];
        xdata_sortedn{iclr} = [xdata_sortedn{iclr},Totaldps];

    end
end
clearvars problem
bs = [];
for iclr = 1:length(Models)
    if ~isempty(data_sortedn{iclr})
    y = log10(data_sortedn{iclr});
    x = log10(xdata_sortedn{iclr});
    X = [ones(length(x),1),x'];
    b1 = X\y';
    bs = [bs,b1]
    ycalc = [ones(length(Totaldps),1),log10(Totaldps)']*b1;
    plot(Totaldps,10.^ycalc',lsts{i},'color',colors(iclr,:))

    
    Rsq1(i,iclr) = 1 - sum((y' - X*b1).^2)/sum((y - mean(y)).^2);

    
    
%     xdata = 10.^x; ydata = 10.^y;
%     [xData, yData] = prepareCurveData( xdata, ydata );

    % Set up fittype and options.
%     ft = fittype( 'a*x^b+c', 'independent', 'x', 'dependent', 'y' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.Lower = [-Inf -Inf 0];
%     opts.StartPoint = [0.0893911202344084 0.72000945452844 0.0557319445106051];
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     Rsq2(i,iclr) = gof.rsquare;
%     mycoef = coeffvalues(fitresult);
%     fun = @(x) mycoef(1)*x.^mycoef(2)+mycoef(3);
%     plot(Totaldps,fun(Totaldps),'-.','color',colors(iclr,:))
    end
end

hold off

set(gca,{'xscale','yscale'},{'log','log'})
set(gca,'box','off')
ylim(10.^[-3,4])
xlim(10.^[1.9,5.1])
set(gca,{'xtick'},{10.^[2,3,4,5]})
set(gca,{'ytick'},{10.^[-2,0,2,4]})
%set(gca,{'ytick','yticklabel'},{10.^[-2,0,2,4],{'0.01','1','100','10000'}})
title(ROIS{i})
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,8.5,5.1],'centimeters',[1+8.5,3+5.1],'Painters'})
%xlabel('# dipoles');
ylabel('A_{max} [\mum]')
set(gca,'position',[0.2341    0.2073    0.6709    0.6891])
set(get(gca,'ylabel'),'position', [25.6960    3.1623   -1.0000]);
end

%%
%Collect all info
clear all
intm1 = load('./Results/SCs_all_cEG_Amp_nos4l_fkApdpset_a2piCondBug_v2.mat');
intm2 = load('./Results/SCs_all_cEG_Amp_NCfixed2_nos4l_fkApdpset_v2.mat');
intm3 = load('./Results/SCs_all_cEG_Amp_NCfixed2_sOSC_nos4l_fkApdpset_v2.mat');
intm4 = load('./Results/SCs_all_cEG_Amp_NCfixed3_sOSC_nos4l_fkApdpset_v2.mat');

dps1 = intm1.Inputall(1).Input.Totaldps;
dps2 = intm2.Inputall(1).Input.Totaldps;
dps3 = intm3.Inputall(1).Input.Totaldps;
dps4 = intm4.Inputall(1).Input.Totaldps;

dps = unique([dps1,dps2,dps3,dps4]);

RequiredAmp1 = intm1.requiredAmpAll;
RequiredAmp2 = intm2.requiredAmpAll;
RequiredAmp3 = intm3.requiredAmpAll;
RequiredAmp4 = intm4.requiredAmpAll;

M = max([size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1),size(RequiredAmp1,1)]);
RequiredAmp = nan(M,length(dps),4);


nLayers = [3,5,10];
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
ROIs = {'Deep','Cortex'};
[nL_mat,ROI_mat,Models_mat] = ndgrid(nLayers,ROIs,Models);
nL_mat = nL_mat(:)'; ROI_mat = ROI_mat(:)'; Models_mat = Models_mat(:)';

[~,idx_dps,idx_dps2] = intersect(dps,dps1);
for i=1:length(intm1.Inputall)
    Settings = intm1.Inputall(i).Input.Settings;
    info1(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info1(i).Model = intm1.Inputall(i).Input.Model;
    info1(i).ROI = intm1.Inputall(i).Input.ROI;
    
    
    idx = info1(i).nLayers==nL_mat & strcmpi(info1(i).Model,Models_mat) &...
        strcmpi(info1(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,1) = RequiredAmp1(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps2);
for i=1:length(intm2.Inputall)
    Settings = intm2.Inputall(i).Input.Settings;
    info2(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info2(i).Model = intm2.Inputall(i).Input.Model;
    info2(i).ROI = intm2.Inputall(i).Input.ROI;
    
    idx = info2(i).nLayers==nL_mat & strcmpi(info2(i).Model,Models_mat) & strcmpi(info2(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,2) = RequiredAmp2(i,idx_dps2);
end

[~,idx_dps,idx_dps2] = intersect(dps,dps3);
for i=1:length(intm3.Inputall)
    Settings = intm3.Inputall(i).Input.Settings;
    info3(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info3(i).Model = intm3.Inputall(i).Input.Model;
    info3(i).ROI = intm3.Inputall(i).Input.ROI;
    
    idx = info3(i).nLayers==nL_mat & strcmpi(info3(i).Model,Models_mat) & strcmp(info3(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,3) = RequiredAmp3(i,idx_dps2);
    
end

[~,idx_dps,idx_dps2] = intersect(dps,dps4);
for i=1:length(intm4.Inputall)
    Settings = intm4.Inputall(i).Input.Settings;
    info4(i).nLayers = Settings{find(strcmpi(Settings,'nLayers'))+1};
    info4(i).Model = intm4.Inputall(i).Input.Model;
    info4(i).ROI = intm4.Inputall(i).Input.ROI;
    
    idx = info4(i).nLayers==nL_mat & strcmpi(info4(i).Model,Models_mat) & strcmpi(info4(i).ROI,ROI_mat);
    RequiredAmp(idx,idx_dps,4) = RequiredAmp4(i,idx_dps2);
end

row1 = {'2','3','3','4'};
row2 = {'1e-11','1e-14','1e-14','1e-17'};
xvals = [1,3,4];
labelArray = [row1(xvals);row2(xvals)];
labelArray = strjust(pad(labelArray),'center');
tickLabels = sprintf('  %s\\newline%s\n', labelArray{:});

types = {'-','-.',':','--'};
Markers = {'o','x','d','s'};

Colors = thermal(4);

close all
ytarget = 0.05;
dps_nr = [5,7,9];
YVAL = [];
Rsqval = [];
alphaval_ytarget = [];
Model_rval = [];
PSDval_ytarget = [];
bs = [];
bs2 = [];

for isubplot = 1:length(Models)
    figure()
    for iM = 1:length(Models_mat)
        for idps = 1:length(dps_nr)        
            
            idx_c = idps;
            idx_lt = find(strcmpi(ROI_mat{iM},ROIs));
            idx_m = find(strcmpi(ROI_mat{iM},ROIs));
            ltype = [Markers{isubplot},types{idx_lt}];
            yval = squeeze(RequiredAmp(iM,dps_nr(idps),xvals));
            yval = yval*1e6; % convert to um
            
            flag = strcmpi(Models_mat{iM},Models{isubplot}) && nL_mat(iM) == 3;
            
            if length(yval(~isnan(yval)))>1 && flag

                yintm = log10(yval'/10^mean(log10(yval)));
                %yintm = log10(yval'/yval(1));
                YVAL = [YVAL;yintm];
                Model_rval = [Model_rval;[isubplot,iM]];
                b = [ones(numel(yval),1),[2;3;4]]\yintm';
                bs2 = [bs2,b];
                Rval = 1-sum((yintm'-[ones(numel(yval),1),[2;3;4]]*b).^2)/sum((yintm'-mean(yintm)).^2);
                Rsqval = [Rsqval;Rval];
                
                % find x for y =100nm
                yintm = log10(yval');
                X0 = [ones(numel(yval),1),[2;3;4]];
                b = X0\yintm';
                bs = [bs,b];
                xval_ytarget = (log10(ytarget)-b(1))/b(2);
                PSDval = calcPSD1MHz(xval_ytarget);
                alphaval_ytarget = [alphaval_ytarget;xval_ytarget];
                PSDval_ytarget = [PSDval_ytarget;PSDval];

                plot(1:length(xvals),yval,ltype(1),'color',Colors(idx_c,:),'HandleVisibility','off');
                hold on
                plot(1:length(xvals),10.^([ones(numel(yval),1),[2;3;4]]*b),ltype(2:end),'color',Colors(idx_c,:),'HandleVisibility','off');
            end           
            
            ax = gca;
            ax.XLim = [0.8,length(xvals)+0.2];
            ax.XTick = 1:length(xvals);
            ax.XTickLabel = tickLabels;
            
        end
    end
    mytitle = split(Models{isubplot},'_');
    if strcmpi(mytitle{1},'mouse')
        mytitle{2} = 'scalp';
    end
    title(sprintf('%s_{%s}',mytitle{:}))
    set(gca,{'box','xscale','yscale'},{'off','lin','log'})
    set(gca,{'ytick'},{10.^[-4,-2,0,2,4]})
    ylim(10.^[-4,4])
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
%%
% linear regresion
XVAL = repmat([2,3,4],size(YVAL,1),1);
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

plot([2,3,4],[ones(3,1),[2,3,4]']*b,'k-')
hold off
xlim([1.8,4.2])

%%
ytarget = 0.05;
dps_nr = [3,5,7,9];
YVAL = [];
Rsqval = [];
alphaval_ytarget = [];
Model_rval = [];
PSDval_ytarget = [];
bs = [];
figure
for isubplot = 1:length(dps_nr)
    for iM = 1:length(Models_mat)        
        
    idx_c = find(strcmpi(Models_mat{iM},Models));
    idx_lt = find(strcmpi(ROI_mat{iM},ROIs));
    idx_m = find(nL_mat(iM)==nLayers);
    ltype = [Markers{idx_m},types{idx_lt}];
    yval = squeeze(RequiredAmp(iM,dps_nr(isubplot),xvals));
    yval = yval*1e6; % convert to um
    
    if length(yval(~isnan(yval)))>1
        yintm = log10(yval'/mean(yval));
        yintm = log10(yval'/10^mean(log10(yval)));
        %yintm = log10(yval'/yval(1));
        YVAL = [YVAL;yintm];
        Model_rval = [Model_rval;[isubplot,iM]];
        b = [ones(numel(yval),1),[1;2;3]]\yintm';
        Rval = 1-sum((yintm'-[ones(numel(yval),1),[1;2;3]]*b).^2)/sum((yintm'-mean(yintm)).^2);
        Rsqval = [Rsqval;Rval];
        
        % find x for y =100nm
        yintm = log10(yval');
        X0 = [ones(numel(yval),1),[2;3;4]];
        b = X0\yintm';
        bs = [bs,b];
        xval_ytarget = (log10(ytarget)-b(1))/b(2);
        PSDval = calcPSD1MHz(xval_ytarget);
        alphaval_ytarget = [alphaval_ytarget;xval_ytarget];
        PSDval_ytarget = [PSDval_ytarget;PSDval];
    end
    subplot(2,2,isubplot)
    p = plot(1:length(xvals),log10(yval),ltype,'color',Colors(idx_c,:),'markersize',10);
    if iM==1; hold on; end
    
    %set(gca,{'xscale','yscale'},{'log','log'})
    %set(gca,{'yscale'},{'log'});
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ax = gca;
    ax.XLim = [0,length(xvals)+1];
    ax.XTick = 1:length(xvals);
    ax.XTickLabel = tickLabels;
    
    end
    title(['nr. dipoles =',  num2str(dps(dps_nr(isubplot)))])
end
subplot(2,2,4)
for i=1:length(Models)
    plot(nan,nan,'-','color',Colors(i,:),'DisplayName',[Models{i}(1:5),'_{',Models{i}(7:end),'}']);
end
for i=1:length(ROIs)
    plot(nan,nan,types{i},'color','k','DisplayName',ROIs{i});
end
for i=1:length(nLayers)
    plot(nan,nan,Markers{i},'color','k','DisplayName',sprintf('# Layers = %i',nLayers(i)));
end


% linear regresion
XVAL = repmat([1,2,3],size(YVAL,1),1);
b = [ones(numel(XVAL),1),XVAL(:)]\YVAL(:)

Rsq = 1-sum((YVAL(:)-[ones(numel(XVAL),1),XVAL(:)]*b).^2)/sum((YVAL(:)-mean(YVAL(:))).^2)
Rsq_min = min(Rsqval)
Rsq_max = max(Rsqval)
figure
scatter(XVAL(:),YVAL(:));
hold on
plot([1,2,3],[ones(3,1),[1,2,3]']*b,'k-')
hold off
xlim([0,4])
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
