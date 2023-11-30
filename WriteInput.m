% Generate input files
FileName = './Inputs/InputBiologicalNoisev3';
fus = 1e6; fusMat = repmat(fus,1,13); fusMat([2,3]) = [1e5,1e7];
meanI = 10; meanIMat = repmat(meanI,1,13); meanIMat([5,6,7]) = [50,50,100];
stdI = 3; stdIMat = repmat(stdI,1,13); stdIMat([4,5,6,7,8]) = [1,10,25,log(3.5),log(3.5)];
ThermalNoiseAmp = 0; ThermalNoiseMat = repmat(ThermalNoiseAmp,1,13); ThermalNoiseMat(13) = 10e-2; %µV
dpDistribution = 'fibonaccisingle'; dpDistributionMat = repmat(cellstr(dpDistribution),1,13); 
dpDistributionMat([9,10,11,12]) = [{'fibonaccimultiple'},{'fibonaccimultiple'},{'random1'},{'random1'}];
dpOrientation = 'radial'; dpOrientationMat = repmat(cellstr(dpOrientation),1,13); dpOrientationMat([10,12]) = {'random'};
dpI_space = 'normal'; dpISpaceMat = repmat(cellstr(dpI_space),1,13); dpISpaceMat([7,8]) = {'lognormal'};
Models={'Human scalp','Human cortical','Mouse scalp'};
Solutiontypes = {'4SphereS8.7R~f','4SphereS8.7R~f'};
RPOI = [0.082,0.07,0.0059];
for isettings=1:13
    Input = struct();
Amp = logspace(0,5,6).*1e-9;
Totaldps = round(logspace(0,5,11)); %round(logspace(0,5+double(isettings>=9 & isettings<=12),11+2*+double(isettings>=9 & isettings<=12)));
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
  Totaldps(i) = round(Totaldps(i),orders(i));
end
if isettings ==3
    dps_run = min(1e2,Totaldps);
else
dps_run = min(1e3,Totaldps);
end
Settings = {'fus',fusMat(isettings),'meanI',meanIMat(isettings),'stdI',stdIMat(isettings),...
    'ThermalNoiseAmp',ThermalNoiseMat(isettings),'dpDistribution',dpDistributionMat{isettings},...
    'dpOrientation',dpOrientationMat{isettings},'dpI_space',dpISpaceMat{isettings}};
Fignr = 1+10*(isettings-1);
SettingsStr = ['Settingset: ',num2str(isettings)];
Input.Amp = Amp; Input.Totaldps = Totaldps; Input.dps_run = dps_run; Input.Settings = Settings; Input.Fignr = Fignr;
Input.SettingsStr = SettingsStr; Input.HPC_flag = 1; Input.Models = Models;
Input.Solutiontypes = Solutiontypes; Input.RPOI = RPOI;
savenames{isettings,1} = ['Inputv3_',num2str(isettings)];
save(['./Inputs/',savenames{isettings,1},'.mat'],'Input');
end

% Convert cell to a table and use first row as variable names
T = cell2table(savenames,'VariableNames',{'Input'});
 
% Write the table to a CSV file
writetable(T,[FileName,'.csv'])
%% single model run
FileName = './Inputs/InputBiologicalNoiseSinglev2';
fus = 1e6; fusMat = repmat(fus,1,15); fusMat(1:3) = 1e7;
meanI = 10; meanIMat = repmat(meanI,1,15); 
stdI = 3; stdIMat = repmat(stdI,1,15); 
ThermalNoiseAmp = 0; ThermalNoiseMat = repmat(ThermalNoiseAmp,1,15); 
dpDistribution = 'fibonaccisingle'; dpDistributionMat = repmat(cellstr(dpDistribution),1,15); 
dpDistributionMat([4:end]) = [repmat({'fibonaccimultiple'},1,6),repmat({'random1'},1,6)];
dpOrientation = 'radial'; dpOrientationMat = repmat(cellstr(dpOrientation),1,15); dpOrientationMat([7:9,13:15]) = {'random'};
dpI_space = 'normal'; dpISpaceMat = repmat(cellstr(dpI_space),1,15);
Models=repmat({'Human scalp','Human cortical','Mouse scalp'},1,5);
Solutiontypes = repmat({'4SphereS8.7R25','4SphereS8.7R25','Mouse4Sphere~fair0'},1,5);
RPOI = repmat([0.082,0.07,0.0059],1,5);
%%
for isettings=1:15
    Input = struct();
Amp = logspace(0,5,6).*1e-9;
Totaldps = round(logspace(0,5+double(isettings>=4 & isettings<=15),11+2*+double(isettings>=4 & isettings<=15)));
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
  Totaldps(i) = round(Totaldps(i),orders(i));
end
if any(isettings ==[1:3])
    dps_run = min(1e2,Totaldps);
else
dps_run = min(1e3,Totaldps);
end
Settings = {'fus',fusMat(isettings),'meanI',meanIMat(isettings),'stdI',stdIMat(isettings),...
    'ThermalNoiseAmp',ThermalNoiseMat(isettings),'dpDistribution',dpDistributionMat{isettings},...
    'dpOrientation',dpOrientationMat{isettings},'dpI_space',dpISpaceMat{isettings}};
Fignr = 1+10*(isettings-1);
SettingsStr = ['Settingset: ',num2str(isettings)];
Input.Amp = Amp; Input.Totaldps = Totaldps; Input.dps_run = dps_run; Input.Settings = Settings; Input.Fignr = Fignr;
Input.SettingsStr = SettingsStr; Input.HPC_flag = 1; Input.Models = Models{isettings};
Input.Solutiontypes = Solutiontypes{isettings}; Input.RPOI = RPOI(isettings);
savenames{isettings,1} = ['Inputsingle',num2str(isettings)];
save(['./Inputs/',savenames{isettings,1},'.mat'],'Input');
end

% Convert cell to a table and use first row as variable names
T = cell2table(savenames,'VariableNames',{'Input'});
 
% Write the table to a CSV file
writetable(T,[FileName,'.csv'])

%% write inputs for study with pressure fields (single model run)
locHdeep = 'D:\no backup\Sim4Life\3spherical\CenterFocus\Wicasim3_sHHM\';
locHcort = 'D:\no backup\Sim4Life\3spherical\CortexFocus\Results_SEFT\';
fnHdeep = 'ExportedData_SHHM_AW164_gs0375mm';
fnHcort = 'ExportedData_sHHM_1SEFT_fp0_55_0_gsrs0375mm';
locMdeep = 'D:\no backup\Sim4Life\3spherical\CenterFocus\Absorptioncoefficient effect\';
fnMdeep = 'ExportedData_MM_changeACSkull.mat';
locMcort = 'D:\no backup\Sim4Life\3spherical\CortexFocus\Special_focus_SEFT_Mouse.smash_Results\';
fnMcort = 'ExportedData_sMM_1SEFT_fp0_4_0';



FileName = './Inputs/Input_Mosc_sM_v3';
tnoI = 18; %total number of Inputs
fus = 1e6; fusMat = repmat(fus,1,tnoI);
meanI = 10; meanIMat = repmat(meanI,1,tnoI); 
stdI = 3; stdIMat = repmat(stdI,1,tnoI);  
dpDistributionMat = repmat(horzcat(repmat({'fibonaccisingle'},1,2),repmat({'fibonaccimultiple'},1,2),...
    repmat({'s4l_hotspot'},1,2)),1,3);
dpOrientation = 'radial'; dpOrientationMat = repmat(cellstr(dpOrientation),1,tnoI);
dpI_space = 'normal'; dpISpaceMat = repmat(cellstr(dpI_space),1,tnoI);
Models = horzcat(repmat({'human_scalp'},1,6),repmat({'human_cortical'},1,6),...
    repmat({'mouse_fair0'},1,6));
ROIs = repmat({'Deep','Cortex'},1,tnoI/2);
S4l_locations = horzcat(repmat(horzcat(cellstr(locHdeep),cellstr(locHcort)),1,6),...
    repmat(horzcat(cellstr(locMdeep),cellstr(locMcort)),1,3));
S4l_Filenames = horzcat(repmat(horzcat(cellstr(fnHdeep),cellstr(fnHcort)),1,6),...
    repmat(horzcat(cellstr(fnMdeep),cellstr(fnMcort)),1,3));
S4l_resolutions = repmat([1,3,1]*10^(-3),tnoI,1);
%%
for isettings=1:tnoI
    Input = struct();
    Amp = logspace(2,4,3).*1e-9;
    Totaldps = round(logspace(1,5,9));
    orders = -floor(log10(Totaldps));
    for i=1:length(Totaldps)
        Totaldps(i) = round(Totaldps(i),orders(i));
    end
    dps_run = 100*ones(size(Totaldps));
    
Settings = {'fus',fusMat(isettings),'meanI',meanIMat(isettings),'stdI',stdIMat(isettings),...
    'dpDistribution',dpDistributionMat{isettings},...
    'dpOrientation',dpOrientationMat{isettings},'dpI_space',dpISpaceMat{isettings}};
Fignr = 1;
SettingsStr = ['Settingset: ',num2str(isettings)];
Input.Model = Models{isettings};
Input.Amp = Amp; Input.Totaldps = Totaldps; Input.dps_run = dps_run; 
Input.Settings = Settings; Input.Fignr = Fignr;
Input.SettingsStr = SettingsStr; 
Input.HPC_flag = 1; Input.Models = Models{isettings};
Input.S4l_flag = 1;
Input.Location_S4l = S4l_locations{isettings};
Input.Filename_S4l = S4l_Filenames{isettings};
Input.ROI = ROIs{isettings};
Input.S4l_res = S4l_resolutions(isettings,:);
savenames{isettings,1} = [FileName,'_',num2str(isettings)];
save([savenames{isettings,1},'.mat'],'Input');

end

% Convert cell to a table and use first row as variable names
T = cell2table(savenames,'VariableNames',{'Input'});
 
% Write the table to a CSV file
writetable(T,[FileName,'.csv'])


%% write inputs for study with ratio differences (single model run)
FileName = './Inputs/Input_Mosc_sM_nos4l_v1';
tnoI = 18; %total number of Inputs
fus = 1e6; fusMat = repmat(fus,1,tnoI);
meanI = 10; meanIMat = repmat(meanI,1,tnoI); 
stdI = 3; stdIMat = repmat(stdI,1,tnoI);  
dpDistributionMat = repmat(horzcat(repmat({'fibonaccisingle'},1,2),repmat({'fibonaccimultiple'},1,2),...
    repmat({'s4l_hotspot'},1,2)),1,3);
dpOrientation = 'radial'; dpOrientationMat = repmat(cellstr(dpOrientation),1,tnoI);
dpI_space = 'normal'; dpISpaceMat = repmat(cellstr(dpI_space),1,tnoI);
Models = horzcat(repmat({'human_scalp'},1,6),repmat({'human_cortical'},1,6),...
    repmat({'mouse_fair0'},1,6));
ROIs = repmat({'Deep','Cortex'},1,tnoI/2);
Aus_way = repmat({'kexp'},1,tnoI);
%%
for isettings=1:tnoI
    Input = struct();
    kAus = [70,100,700];
    orders = -floor(log10(kAus));
    for i=1:length(kAus)
        kAus(i) = round(kAus(i),orders(i));
    end
    Amp = 1e-6*ones(size(kAus));
    Totaldps = round(logspace(1,5,9));
    orders = -floor(log10(Totaldps));
    for i=1:length(Totaldps)
        Totaldps(i) = round(Totaldps(i),orders(i));
    end
    dps_run = 100*ones(size(Totaldps));
    
Settings = {'fus',fusMat(isettings),'meanI',meanIMat(isettings),'stdI',stdIMat(isettings),...
    'dpDistribution',dpDistributionMat{isettings},...
    'dpOrientation',dpOrientationMat{isettings},'dpI_space',dpISpaceMat{isettings}};
Fignr = 1;
SettingsStr = ['Settingset: ',num2str(isettings)];
Input.Model = Models{isettings};
Input.Aus_way = Aus_way{isettings};
Input.kAus = kAus;
Input.Amp = Amp; Input.Totaldps = Totaldps; Input.dps_run = dps_run; 
Input.Settings = Settings; Input.Fignr = Fignr;
Input.SettingsStr = SettingsStr; 
Input.HPC_flag = 1; Input.Models = Models{isettings};
Input.S4l_flag = 0;
Input.ratioAmpflag = 1;
Input.ROI = ROIs{isettings};
savenames{isettings,1} = [FileName,'_',num2str(isettings)];
save([savenames{isettings,1},'.mat'],'Input');

end

% Convert cell to a table and use first row as variable names
T = cell2table(savenames,'VariableNames',{'Input'});
 
% Write the table to a CSV file
writetable(T,[FileName,'.csv'])

