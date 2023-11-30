%s=rng;
close all
clear all
clc
addpath(genpath('Functions'))
s = load('randomseed');
s=s.s;
%% test function and parallelcomp
%rng(s)
% to test parallelcomp set i = 1;
for i=0
    
    
tstart = tic;
[RMSpc,SNRpc] = investBiologicalNoise('posDp',[0,0,0.065;0,0.05,0],'Totaldps',100,'dps_run',10^2,...
    'Aus',500e-9,'PLOT',1,'ParallelCompute',i,'randomseed',s);
Time.PC(i+1) = toc(tstart);
    
%%
Y=[[RMSpc.RMS]',[RMSpc.RMSptherm]',[RMSpc.RMSDOI]',[RMSpc.RMSDOIvr]',[RMSpc.RMSDOIvrpOSCvr]',[RMSpc.RMSDOIvrpstatnoise]'];
%Y=[[RMSpc.RMS]'];
figure
subplot(1,2,1)
bar(Y)
set(gca,'xticklabel',{RMSpc.info})
legend('RMS','RMSptherm','RMSDOI','RMSDOIvr','RMSDOIvrpOSCvr','RMSDOIvrpstatnoise')

YSNR = [SNRpc.DOIvr_NoiseAll(2,:)',SNRpc.DOIvr_NoiseAllptherm(2,:)',SNRpc.DOIvr_NoiseOSCvr(2,:)',SNRpc.DOIvr_NoiseStatic(2,:)']
title('RMS')
subplot(1,2,2)
bar(YSNR)
set(gca,'xticklabel',{RMSpc.info})
title('SNR')
legend('SNR_{all}','SNR_{all+therm}','SNR_{OSCvr}','SNR_{static}')
end

% mean I =10µA sigma = 3µA
continue_flag = input('continue runing script? \n');
if ~continue_flag
    return
end
%% validate calcErrorGrid (new) vs calcStrength Curve

fus = 1e6; fusMat = repmat(fus,1,13); fusMat([2,3]) = [1e5,1e7];
meanI = 10; meanIMat = repmat(meanI,1,13); meanIMat([5,6,7]) = [50,50,100];
stdI = 3; stdIMat = repmat(stdI,1,13); stdIMat([4,5,6,7,8]) = [1,10,25,log(10),log(10)];
ThermalNoiseAmp = 0; ThermalNoiseMat = repmat(ThermalNoiseAmp,1,13); ThermalNoiseMat(13) = 10e3;
dpDistribution = 'fibonaccisingle'; dpDistributionMat = repmat(cellstr(dpDistribution),1,13); 
dpDistributionMat([9,10,11,12]) = [{'fibonaccimultiple'},{'fibonaccimultiple'},{'fibonaccimultiple'},{'fibonaccimultiple'}];
dpOrientation = 'radial'; dpOrientationMat = repmat(cellstr(dpOrientation),1,13); dpOrientationMat([10,12]) = {'random'};
dpI_space = 'normal'; dpISpaceMat = repmat(cellstr(dpI_space),1,13); dpISpaceMat([7,8]) = {'lognormal'};

isettings = 9;
Amp = logspace(0,4,3)*1e-9;
Totaldps = round(logspace(1,3,3)); 
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
  Totaldps(i) = round(Totaldps(i),orders(i));
end
Totaldps = [10,100,300];
dps_run = min(1e2,Totaldps);
Settings = {'fus',fusMat(isettings),'meanI',meanIMat(isettings),'stdI',stdIMat(isettings),...
    'ThermalNoiseAmp',ThermalNoiseMat(isettings),'dpDistribution',dpDistributionMat{isettings},...
    'dpOrientation',dpOrientationMat{isettings},'dpI_space',dpISpaceMat{isettings}};
Fignr = get(gcf,'number')+1;
[requiredAmp,requiredAmplower,itermat] = ...
    calcStrengthCurve2(Amp,Totaldps,dps_run,Settings,Fignr,['Settingset: ',num2str(isettings)]);
%% test if same result with calcErrorGrid
% all functionalities are encorporated in the calcErrorGrid function:
% Validate
Input.Settings = Settings; Input.Amp = Amp; Input.Totaldps = Totaldps; Input.dps_run = dps_run;
Input.evalMode = 'normal'; Input.vibrField = 'none'; Input.ROI = 'cortex'; Input.Model = 'human_scalp';
Input.SettingsStr = ['Settingset: ',num2str(isettings)]; Input.calcStrengthCurve = 1;
[RMSall,SNRall,Outall,TSTOPall,Totaldps,requiredAmp2,requiredAmplower2,itermat2] = calcErrorGrid(Input);


%% run study for 13 setting sets. Look to table in raportations first months 50819-weeklyupdate week15
% see end first 4 months
DEBUG = 1;

fus = 1e6; fusMat = repmat(fus,1,13); fusMat([2,3]) = [1e5,1e7];
meanI = 10; meanIMat = repmat(meanI,1,13); meanIMat([5,6,7]) = [50,50,100];
stdI = 3; stdIMat = repmat(stdI,1,13); stdIMat([4,5,6,7,8]) = [1,10,25,log(10),log(10)];
ThermalNoiseAmp = 0; ThermalNoiseMat = repmat(ThermalNoiseAmp,1,13); ThermalNoiseMat(13) = 10e3;
dpDistribution = 'fibonaccisingle'; dpDistributionMat = repmat(cellstr(dpDistribution),1,13); 
dpDistributionMat([9,10,11,12]) = [{'fibonaccimultiple'},{'fibonaccimultiple'},{'fibonaccimultiple'},{'fibonaccimultiple'}];
dpOrientation = 'radial'; dpOrientationMat = repmat(cellstr(dpOrientation),1,13); dpOrientationMat([10,12]) = {'random'};
dpI_space = 'normal'; dpISpaceMat = repmat(cellstr(dpI_space),1,13); dpISpaceMat([7,8]) = {'lognormal'};

for isettings=1:13
    if DEBUG
        Amp = logspace(0,4,3)*1e-9;
        Totaldps = round(logspace(1,3,3));
        
    else
        
        Amp = logspace(0,5,6).*1e-9;
        Totaldps = round(logspace(0,4,5));%round(logspace(0,5+double(isettings>=9 & isettings<=12),6+2*+double(isettings>=9 & isettings<=12)));
    end
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
  Totaldps(i) = round(Totaldps(i),orders(i));
end
if isettings == 3
    dps_run = min(1e2,Totaldps);
else
dps_run = min(1e3,Totaldps);
end


Settings = {'fus',fusMat(isettings),'meanI',meanIMat(isettings),'stdI',stdIMat(isettings),...
    'ThermalNoiseAmp',ThermalNoiseMat(isettings),'dpDistribution',dpDistributionMat{isettings},...
    'dpOrientation',dpOrientationMat{isettings},'dpI_space',dpISpaceMat{isettings}};
Fignr = get(gcf,'number')+1;
Out(isettings).Settings = Settings;
[Out(isettings).requiredAmp,Out(isettings).requiredAmplower,Out(isettings).itermat] = calcStrengthCurve2(Amp,Totaldps,dps_run,Settings,Fignr,['Settingset: ',num2str(isettings)]);

end

%%
function [requiredAmp,requiredAmplower,itermat] = calcStrengthCurve2(Amp,Totaldps,dps_run,Settings,Fignr,SettingsStr)
[X,Y] = meshgrid(Totaldps,Amp);
Y = Y*1e9;
txtbxstr = gentxtbxstr(Settings);
Models={'Human scalp'};%,'Human air','Human sub scalp','Human cortical','Mouse scalp','Mouse air'};
Solutiontypes = {'4SphereS8.7R25','4SphereS8.7R25','4SphereS8.7R25','4SphereS8.7R25','Mouse4Sphere~fair0','Mouse4Sphere~fair1'};
RPOI = [0.082,0.087,0.075,0.07,0.0059,0.0069];
maxIter = 10;
Thresholds = [0.1,0.05];
maxosc = Amp(end);
totaliters = length(Thresholds)*length(Totaldps)*maxIter;
try
s = load('randomseed.mat');
s = s.s;
catch
    warning('new random seed')
    s = rng;
end
for imodel = 1:length(Models)
for iAmp = 1:length(Amp)
    for idps = 1:length(Totaldps)
        rng(s)
        input_sets = horzcat({'Totaldps',Totaldps(idps),'dps_run',dps_run(idps),...
            'Aus',Amp(iAmp),'ParallelCompute',dps_run(idps)>100,'SolutionType',Solutiontypes{imodel},...
            'POIs',RPOI(imodel)*[0,0,1],'scale',0,'PLOT',0,'showSphere',0,'randomseed',s},Settings);
        [RMS,SNR] = investBiologicalNoise(input_sets{:});
        RMSnMAT(imodel,iAmp,idps) = RMS(1).RMS; RMSDOIvr(imodel,iAmp,idps) = RMS(1).RMSDOIvr;
        SNRMAT(imodel,iAmp,idps) = SNR.DOIvr_NoiseAll(2,1);
    end
end
%% Plot results
figure(Fignr+floor((imodel-1)/3))
subplot(min(3,length(Models)),3,mod(imodel-1,3)*3+1)
surf(X,Y,squeeze(RMSnMAT(imodel,:,:)),'EdgeColor','none','FaceColor','interp')
set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('# dipoles')
xlim([min(Totaldps),max(Totaldps)])
ylabel({['\bf', Models{imodel}],'','Vibration amplitude [nm]'})
colormap('jet')
colorbar
view(0,90)
title('RMS error')
subplot(min(3,length(Models)),3,2+mod(imodel-1,3)*3)
surf(X,Y,squeeze(RMSDOIvr(imodel,:,:)),'EdgeColor','none','FaceColor','interp')
set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('# dipoles')
xlim([min(Totaldps),max(Totaldps)])
ylabel('Vibration amplitude [nm]')
colormap('jet')
colorbar
view(0,90)
title('RMS error')
subplot(min(3,length(Models)),3,3+mod(imodel-1,3)*3)
surf(X,Y,squeeze(SNRMAT(imodel,:,:)),'EdgeColor','none','FaceColor','interp')
set(gca,'xscale','log'); set(gca,'yscale','log');
xlabel('# dipoles')
ylabel('Vibration amplitude [nm]')
minlim = Totaldps(Totaldps>1);
xlim([minlim(1),max(Totaldps)])
colormap('jet')
colorbar
view(0,90)
title('SNR')
set(gcf,'position',[-1919,41,1920,963])
if ~mod(imodel,3) || imodel == length(Models)
mtit(SettingsStr,'xoff',-0.01,'yoff',.05,'fontsize',14)
end
pause(0.01);
%%
    figure(Fignr+floor((length(Models)-1)/3)+1)
    sp{imodel} = subplot(ceil(length(Models)/2),1+double(length(Models)>1),imodel);
    h = waitbar(0,'start');
    hiter = 0;
for ithresh=1:length(Thresholds)
    RMSOI = Thresholds(ithresh);
for idps = 1:length(Totaldps)
    superval = Amp(RMSnMAT(imodel,:,idps)<=RMSOI);
    iter = 0;
    findhigherampiter = 0;
    if isempty(superval)
        superval = maxosc;
        newval = superval;
        findhigheramp = 1;        
    else
        superval = superval(1);
        findhigheramp = 0;
    end
    lowerval = Amp(RMSnMAT(imodel,:,idps)>=RMSOI);
    if isempty(lowerval)
        lowerval = superval;
    else
        lowerval = lowerval(end);
    end
    flag = abs(superval-lowerval)>1e-9;
    while flag && iter<maxIter && findhigherampiter<2
        rng(s)
        if findhigheramp
            newval = newval*10;
            iter = 0;
            findhigherampiter = findhigherampiter +1;
        else
            newval = 10^((log10(superval)+log10(lowerval))/2);
        end
        input_set = horzcat({'POIs',RPOI(imodel)*[0,0,1],'Totaldps',Totaldps(idps),...
            'dps_run',dps_run(idps),'Aus',newval,'ParallelCompute',dps_run(idps)>100,...
            'SolutionType',Solutiontypes{imodel},'scale',0,'scale',0,'PLOT',0,'showSphere',0,'randomseed',s},Settings);
        [RMS,SNR] = investBiologicalNoise(input_set{:});
        RMSn = RMS(1).RMS;
        if RMSn<=RMSOI
            superval = newval;
            findhigheramp = 0;
        else
            lowerval = newval;
        end
        flag = abs(superval-lowerval)>1e-9;
        waitbar((hiter*maxIter+iter)/totaliters,h, [num2str((hiter*maxIter+iter)/totaliters*100),'%'])
        iter = iter+1;
    end
    hiter = hiter+1;
    requiredAmp(imodel,ithresh,idps) = superval;
    requiredAmplower(imodel,ithresh,idps) = lowerval;
    itermat(imodel,ithresh,idps) = iter;
end
plot(Totaldps,squeeze(requiredAmp(imodel,ithresh,:)*1e9))
hold on
end
hold off
set(gca,'xscale','log');
xlabel('number of dipoles')
ylabel('threshold amplitude [nm]')
title(Models{imodel})
if mod(imodel+1,2) || imodel == length(Models)
legend(arrayfun(@(x) ['RMS threshold = ',num2str(x*100),'%'],Thresholds,'UniformOutput',false))
end
close(h)
end
set(gcf,'position',[-1919,41,1920,963])
pause(0.1)
pos = get(sp{1}, 'position');
dim = pos.*[1 1+0.5*pos(4)/pos(2) 0.5 0.5];
annotation(gcf, 'textbox', dim, 'String', txtbxstr,'FitBoxToText','on')
mtit(SettingsStr,'xoff',-0.001,'yoff',.05,'fontsize',14)
end
function out = gentxtbxstr(Settings)
out = {['fus = ',num2str(Settings{2}*1e-6),' MHz']
    ['mean I = ',num2str(Settings{4}),' µA']
    ['std I = ',num2str(Settings{6}),' µA']
    ['Thermal Noise = ',num2str(Settings{8}),' µV']
    ['dpDistribution = ',Settings{10}]
    ['dpOrientation = ',Settings{12}]
    ['dpI_{space} = ',Settings{14}]};    
    
end