% Create results figure 2 for paper
% results from kAus study (change in spatial constant of ultrasonic field)
close all; clear all; clc
folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1224\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset1_kAus_nos4l_v2_12-28-20_1901';
data_rr =  load(fullfile(folder_rawresult,filename_rawresult));


Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI

[M,N] = size(data_rr.Outall);
kAus = nan(M,N);
for i =1:M
    for j =1:N
        kAus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
    end
end
kAus = 1./kAus(:,1).*1000;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kAus)),'\n'],kAus')

%%
close all
ikaus=1; idp=9; ipoi = 1;
plotmethod = 'norm'%'flip'; %'norm'
for ikaus=3:5
    yvals = squeeze(data_rr.rSVR_mat(ikaus,idp,ipoi,:));
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
    
    plot(squeeze(data_rr.trS_mat(ikaus,idp,ipoi,:))*1000,yvals)
    hold on
end
hold off
%%
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
figure()
bardata = squeeze(data_rr.RMS_mat(kaus,idp,1:7));
b = bar(kAus(kaus),bardata);
yticks = get(gca,'ytick');
yticks(end) = 1.1
set(gca,'XDir','reverse')
f = gcf;


y_break_start = 0.45;
y_break_end = 1.05;
y_break_mid=(y_break_end-y_break_start)./2+y_break_start;
dy = 0.01;

y = bardata(:);

y2=y;
y2(y2>=y_break_end)=y2(y2>=y_break_end)-y_break_end+y_break_start+dy;
y2 = reshape(y2,size(bardata));

figure()
b = bar(y2);
xlim = get(gca,'xlim');

% this can be vectorized
dx=(xlim(2)-xlim(1))./10;
yy=repmat([y_break_start+dy y_break_start+2*dy],1,6);
xx=xlim(1)+dx.*[0:11];
patch([xx(:);flipud(xx(:))], ...
    [yy(:);flipud(yy(:)-dy)], ...
    [.8 .8 .8])
        
set(gca,'xlim',xlim);
y_ticks_end = yticks(yticks>y_break_start)+(-y_break_end+y_break_start+dy);
y_ticks_end(y_ticks_end<=y_break_start) = [];
yticks = [yticks(yticks<=y_break_start),y_ticks_end];
set(gca,'ytick',yticks)

% map back
yticks2 = yticks;
yticks2(yticks>y_break_start)=yticks(yticks>y_break_start)-(-y_break_end+y_break_start+dy);
for i=1:length(yticks2)
   yticklabel{i}=sprintf('%0.1f',yticks2(i));
end
yticklabel{1} = 0;
set(gca,'yticklabel',yticklabel);
set(gca,'ylim',[0,yticks(end)*1.1])
set(gca,'box','off')

set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})

set(gca,'xticklabel',arrayfun(@(x) sprintf('%0.2f',x),kAus(kaus),'UniformOutput',false))
xlabel('\kappa[mm]')
ylabel('RMSE')
hold on
for i =1:length(b)
    b(i).FaceColor = 'flat';
    b(i).CData = colors(i,:).*ones(size(b(i).CData(:,1),2),1);
    bh(i) = bar(nan,nan,'FaceColor',colors(i,:));
end

legend(bh,horzcat(arrayfun(@(x) sprintf('POI_%i',x),1:5,'UniformOutput',false),{'mPOI','mPsO'}),...
    'box','off','location','northoutside','numColumns',7);
close(f)
%% create spheres
rng(1)
myin = data_rr.Outall(ikaus,idp).Out.Param;
sigma = myin.sigma;
SphereRes = 40;
ndps = 50;
[Options,RSphere] = getSettings(data_rr.Outall(ikaus,idp).Out.Param.SolutionType,1);
Options{find(strcmpi('tAir',Options))+1} = 0; % sphere output scalp
RPOI = RSphere+0.017;

[CSource,CSink] = GendpPos(myin.dpDistribution,myin.dpOrientation,myin.d,0.065,...
ndps,1,'nLayers',myin.nLayers);
OSCindices = data_rr.Outall(ikaus,idp).Out.OSCindices;
POIs = data_rr.Outall(ikaus,idp).Out.POIs;
CSource = [myin.posDp+myin.d*myin.OrienDipole;CSource];
CSink = [myin.posDp-myin.d*myin.OrienDipole;CSink];

Totaldps_plot = min(ndps,size(CSource,1));
idx = [1,sort(randperm(size(CSource,1)-1,Totaldps_plot))+1];
CSource_plot = CSource(idx,:);
CSink_plot = CSink(idx,:);

Iarray = [10,normrnd(10,3,1,ndps)];

varargin2 = horzcat(Options,{'POI',POIs,'DOI',1,'OSC',1:Totaldps_plot+1,'scale',0});
PotentialSphere_Multi(CSource_plot,CSink_plot,Iarray,sigma,RSphere,SphereRes,get(gcf,'number')+1,varargin2);
%%
colormap('inferno')
title('')
set(gca,{'color','box'},{'None','off'})
grid off
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
get(gcf,'children')
children = get(gcf,'children');
c = children(2);
l = children(1);
l.String = l.String(1:2);
l.String{2} = 'Vibrating dp';
l.ItemTokenSize(1) = 10
set(l,{'box','fontsize','position'},{'off',10,[0.0066,0.7994,0.3404,0.1565]})

set(gca,'position',[-0.0083    0.0659    0.7472    0.9066])
set(c,{'position','fontsize'},{[0.6987    0.0789    0.0489    0.8542],10})
c.Label.String = 'V_{scalp} [\muV]';
xl = get(gca,'xlabel'); set(xl,'position', [0.0297   -0.0488   -0.1056]);
yl = get(gca,'ylabel'); set(yl,'position', [-0.0108    0.0549   -0.1017]);
zl = get(gca,'zlabel'); set(zl,'position',  [0.0241   -0.0946    0.0174]);
numbers = findall(gca,'type','text');
set(numbers(1),'position',[0.0563 0.0485 0.0563]) %5
set(numbers(2),'position',[-0.0088 0.01 0.1133])
set(numbers(3),'position',[-0.0023 0.0942 0.0042])
set(numbers(4),'position',[-0.0088 -0.0067 0.1133])
set(numbers(5),'position',[0.0933 -0.0024 0.0042]) %1
numbers(2).String = '';
numbers(4).String = '2, 4';
%% second row of figures
% Create results figure 2 for paper
% results from kAus study (change in spatial constant of ultrasonic field)
clear all;
folder_rawresult = 'D:\no backup\EEGUS\HPC_files\1224\Crop';
filename_rawresult = 'ResultsCrop_cEG_Settingset4_kAus_nos4l_v2_12-29-20_0110';
data_rr =  load(fullfile(folder_rawresult,filename_rawresult));

folder_collective = 'D:\users\rschoeters\Documents\Imec USEEG\Matlab\Results';
filename_collective = 'SCs_all_cEG_kAus_nos4l_a2piCondBug_v2';
data_col = load(fullfile(folder_collective,filename_collective));
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI

[M,N] = size(data_rr.Outall);
kAus = nan(M,N);
for i =1:M
    for j =1:N
        kAus(i,j) = data_rr.Outall(i,j).Out.Param.k_Aus;
    end
end
kAus = 1./kAus(:,1).*1000;
fprintf(['spatial constants k: ', repmat('%5.2f ', 1, length(kAus)),'\n'],kAus')

%%
figure()
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
set(gca,'xlim',[0,25])
xlabel('time [ms]')
ylabel('norm. Signal [-]')
set(gca,'box','off')
%l = legend('show','box','off','NumColumns',1);
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
%set(l,'position',[0.5390    0.5414    0.3816    0.3476])
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
figure()
bardata = squeeze(data_rr.RMS_mat(kaus,idp,1:7));
b = bar(bardata);
yticks = get(gca,'ytick');
f = gcf;


y_break_start = 0.45;
y_break_end = 0.61;
y_break_mid=(y_break_end-y_break_start)./2+y_break_start;
dy = 0.01;

y = bardata(:);

y2=y;
y2(y2>=y_break_end)=y2(y2>=y_break_end)-y_break_end+y_break_start+dy;
y2 = reshape(y2,size(bardata));

figure()
b = bar(y2);
xlim = get(gca,'xlim');

% this can be vectorized
dx=(xlim(2)-xlim(1))./10;
yy=repmat([y_break_start+dy y_break_start+2*dy],1,6);
xx=xlim(1)+dx.*[0:11];
patch([xx(:);flipud(xx(:))], ...
    [yy(:);flipud(yy(:)-dy)], ...
    [.8 .8 .8])
        
set(gca,'xlim',xlim);
y_ticks_end = yticks(yticks>y_break_start+dy)+(-y_break_end+y_break_start+dy);
y_ticks_end(y_ticks_end<=y_break_start+dy) = [];
yticks = [yticks(yticks<=y_break_start),y_ticks_end];

yticks(2:2:end) = [];
set(gca,'ytick',yticks)

% map back
yticks2 = yticks;
yticks2(yticks>y_break_start)=yticks(yticks>y_break_start)-(-y_break_end+y_break_start+dy);
for i=1:length(yticks2)
   yticklabel{i}=sprintf('%0.1f',yticks2(i));
end
yticklabel{1} = 0;
set(gca,'yticklabel',yticklabel);
set(gca,'ylim',[0,yticks(end)*1.1])
set(gca,'box','off')

set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})

set(gca,'xticklabel',arrayfun(@(x) sprintf('%0.2f',x),kAus(kaus),'UniformOutput',false))
xlabel('\kappa[mm]')
ylabel('RMSE')
hold on
for i =1:length(b)
    b(i).FaceColor = 'flat';
    b(i).CData = colors(i,:).*ones(size(b(i).CData(:,1),2),1);
    bh(i) = bar(nan,nan,'FaceColor',colors(i,:));
end

legend(bh,horzcat(arrayfun(@(x) sprintf('POI_%i',x),1:5,'UniformOutput',false),{'mPOI','mPsO'}),...
    'box','off','location','northoutside','numColumns',7);
close(f)
%% create spheres
rng(1)
myin = data_rr.Outall(ikaus,idp).Out.Param;
sigma = myin.sigma;
SphereRes = 40;
ndps = 50;
[Options,RSphere] = getSettings(data_rr.Outall(ikaus,idp).Out.Param.SolutionType,1);
Options{find(strcmpi('tAir',Options))+1} = 0; % sphere output scalp
RPOI = RSphere+0.017;

[CSource,CSink] = GendpPos(myin.dpDistribution,myin.dpOrientation,myin.d,0.065,...
ndps,1,'nLayers',myin.nLayers);
OSCindices = data_rr.Outall(ikaus,idp).Out.OSCindices;
POIs = data_rr.Outall(ikaus,idp).Out.POIs;
CSource = [myin.posDp+myin.d*myin.OrienDipole;CSource];
CSink = [myin.posDp-myin.d*myin.OrienDipole;CSink];

Totaldps_plot = min(ndps,size(CSource,1));
idx = [1,sort(randperm(size(CSource,1)-1,Totaldps_plot))+1];
CSource_plot = CSource(idx,:);
CSink_plot = CSink(idx,:);

Iarray = [10,normrnd(10,3,1,ndps)];

varargin2 = horzcat(Options,{'POI',POIs,'DOI',1,'OSC',1:Totaldps_plot+1,'scale',0});
PotentialSphere_Multi(CSource_plot,CSink_plot,Iarray,sigma,RSphere,SphereRes,get(gcf,'number')+1,varargin2);
%%
colormap('inferno')
title('')
set(gca,{'color','box'},{'None','off'})
grid off
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
get(gcf,'children')
children = get(gcf,'children');
c = children(2);
l = children(1);
l.String = l.String(1:2);
l.String{2} = 'Vibrating dp';
l.ItemTokenSize(1) = 10
set(l,{'box','fontsize','position'},{'off',10,[0.0066,0.7994,0.3404,0.1565]})

set(gca,'position',[-0.0083    0.0659    0.7472    0.9066])
set(c,{'position','fontsize'},{[0.6987    0.0789    0.0489    0.8542],10})
c.Label.String = 'V_{scalp} [\muV]';
xl = get(gca,'xlabel'); set(xl,'position', [0.0297   -0.0488   -0.1056]);
yl = get(gca,'ylabel'); set(yl,'position', [-0.0108    0.0549   -0.1017]);
zl = get(gca,'zlabel'); set(zl,'position',  [0.0241   -0.0946    0.0174]);
numbers = findall(gca,'type','text');
set(numbers(1),'position',[0.0563 0.0485 0.0563]) %5
set(numbers(2),'position',[-0.0088 0.01 0.1133])
set(numbers(3),'position',[-0.0023 0.0942 0.0042])
set(numbers(4),'position',[-0.0088 -0.0067 0.1133])
set(numbers(5),'position',[0.0933 -0.0024 0.0042]) %1
numbers(2).String = '';
numbers(4).String = '2, 4';
%%
clear all
close all
folder_collective = 'D:\users\rschoeters\Documents\Imec USEEG\Matlab\Results';
filename_collective = 'SCs_all_cEG_kAus_nos4l_a2piCondBug_v2';
data_col = load(fullfile(folder_collective,filename_collective));
distPOIdps = load(fullfile(folder_collective,'mindistPOIdps.mat'));

Nin = length(data_col.Inputall);

colors = [27,158,119
217,95,2
117,112,179
231,41,138]/255;

%colors = thermal(5);
%colors = colors(1:2:end,:)

Ylims = [min(data_col.requiredAmpAll(:)),max(data_col.requiredAmpAll(:))];

%%
close all
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
labels = {'scalp_{h}','cortical_{h}','air_{h}','scalp_{m}'};
styles = {'x-','o--','*:'};
ROIS = {'deep','cortex'};
NLAYERS = [3,5,10];
colors = [27,158,119
217,95,2
117,112,179
231,41,138]/255;
for iROI=1:length(ROIS)
figure()
for iin = 1:Nin

    myin = data_col.Inputall(iin).Input;
    Totaldps = myin.Totaldps;
    Model = myin.Model;
    ROI = myin.ROI;
    nLayers = myin.Settings{find(strcmpi('nLayers',myin.Settings))+1};
    iclr = find(strcmpi(Models,Model));
    ist = find(NLAYERS==nLayers);
    if strcmpi(ROI,ROIS{iROI})
        plot(Totaldps,1000./data_col.requiredAmpAll(iin,:),styles{ist},'color',colors(iclr,:),'HandleVisibility','off')
        hold on
    end
end

for iclr = 1:length(Models)
    
    plot(nan,nan,'color',colors(iclr,:),'DisplayName',labels{iclr})
end
for ist =1:length(NLAYERS)
    plot(nan,nan,styles{ist},'color',[0.1,0.1,0.1],'DisplayName',sprintf('n = %2.0i',NLAYERS(ist)))
end

%l = legend('show');
xlim([90,1.1e5])
ylim([0.01,100])
set(gca,'xtick',[100,1000,1e4,1e5])
set(gca,'yticklabel',{'0.01','1','100'})
hold off
title(ROIS{iROI})
xlabel('# dipoles')
ylabel('\kappa [mm]')
set(findall(gcf,'type','axes'),'fontsize',10)
set(gca,{'xscale','yscale','box','color'},{'log','log','off','None'})
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.55,5.3],'centimeters',[1+7.55,3+5.3],'Painters'})
end

%%
poidist = [repmat(distPOIdps.mindist_all(1:6,:),3,1);distPOIdps.mindist_all(7:end,:)];
YVALS = nan(2,24,9);
for iROI=1:length(ROIS)
figure()
for iin = 1:Nin

    myin = data_col.Inputall(iin).Input;
    Totaldps = myin.Totaldps;
    Model = myin.Model;
    ROI = myin.ROI;
    nLayers = myin.Settings{find(strcmpi('nLayers',myin.Settings))+1};
    iclr = find(strcmpi(Models,Model));
    ist = find([3,5,10]==nLayers);
    if strcmpi(ROI,ROIS{iROI})
        yval = 1000./data_col.requiredAmpAll(iin,:)./(poidist(iin,:)*1000);
        YVALS(iROI,iin,:) = yval
        plot(Totaldps,yval,styles{ist}(1),'color',colors(iclr,:))
        hold on
    end
end
hold off
xlim([0.9*1e3,1.1e5])
ylim([0,1])
set(gca,'xtick',[1e3,1e4,1e5])
hold off
xlabel('# dipoles')
ylabel('rel. \kappa [-]')
set(findall(gcf,'type','axes'),'fontsize',10)
set(gca,{'xscale','yscale','box','color'},{'log','lin','off','None'})
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.55,5.3],'centimeters',[1+7.55,3+5.3],'Painters'})
end
YVALS = squeeze(nanmean(YVALS,2))

%%
close all
colors = thermal(4)
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
styles = {'x-','o--','*:'};
ROIS = {'deep','cortex'};
poidist = [repmat(distPOIdps.mindist_all(1:6,:),3,1);distPOIdps.mindist_all(7:end,:)];
for i = 1:length(Models)
figure()
for iin = 1:Nin

    myin = data_col.Inputall(iin).Input;
    Totaldps = myin.Totaldps;
    Model = myin.Model;
    ROI = myin.ROI;
    nLayers = myin.Settings{find(strcmpi('nLayers',myin.Settings))+1};
    ist = find(strcmpi(ROIS,ROI));
    iclr = find([3,5,10]==nLayers);
    if strcmpi(Model,Models{i})
        plot(Totaldps,1000./data_col.requiredAmpAll(iin,:),styles{ist},'color',colors(iclr,:))
        hold on
        clr2 = hsv2rgb(min(rgb2hsv(colors(iclr,:)).*[1,10,1],[1,0.1,1]));
        plot(Totaldps,(poidist(iin,:)*1000),styles{ist}(2:end),'color',clr2,'linewidth',1.5)
    end
end
hold off
set(gca,{'xscale','yscale'},{'log','log'})
xlim([1e2,1e5])
end

