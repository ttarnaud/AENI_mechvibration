% Input file for calcErrorGrid
% Necessary inputs fields are:
%   Amp:        different vibration amplitudes
%   Totaldps:   number of dipoles included
%   dp_run:     number of dipoles/run (min(100,totaldps))
%   evalMode:   debug, hpc, local, normal
%   vibrField:  s4l, art(ificial), none
%   sim_method: 'seq','single_sim_totaldps' UPDATE 4/12/2020
%    if art: 
%        Aus_way: spatial distribution vibration amplitudes
%        ('kexp','kscale','normal')
%        kAus: scale or spatial constant, (needs same length Amp!)
%        UPDATE 4/12
%        {if you want to investigate effect kAus make Amp equal length kAus
%        with all values in Amp the same + kAus_con='amp'. If investigate Amp while kAus is
%        fixed for a set of dipoles make kAus same length as totaldps + kAus_con='totaldps'
%             kAus_con: 'none' 'amp' or 'totaldps
%             to which parameter is kAus connected: none = single kAus input
%             Totaldps: kAus is fixed for given nr of dipoles: kAus needs same length as Totaldps
%             Amp: kAus is connected to Amp,kAus same length as Amp  } 
%    if s4l:
%        Filename_S4l
%        Location_S4l
%   ROI:        Region of interest; cortex, deep
%   Model:      human_scalp, human_cortex, human_air, mouse_fair0
%   Settings:   parameters to be changed in invest biological noise
%      eg dpOrientation, dpI_time, dpI_space, dpDistribution, fus, meanI,
%      nLayers, resUS...
%
%Optional:
%   calcStrengthCurve:     flag if calculation of Strength Curve is needed
%                          (iteratief proces to determine necessary amp in
%                          order to have RMS below 10%)
%   maxTime;               maxtime of simulations
%   
%%
FileName = './Inputs/aIfun/Input_cEG';
Ifun = func2str(@(t) alphaTrainFun(t,'aTInput_220308.mat'));
idxdps = 9;
threshold = 0.05;
sets = [1:3:22]; %1:tnoI
suffix = sprintf('kAus_dps%i%i_nos4l_th%i_v2',idxdps(1),idxdps(end),threshold*100);


%look at vibrational interference, Multilayer ({3-}5-10, DOI = Deep,
%Cortex, POI = Scalp, Cortex, Air, Mouse

% variable components
nLayers = [3,5,10];
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
ROIs = {'Deep','Cortex'};
[nL_mat,ROI_mat,Models_mat] = ndgrid(nLayers,ROIs,Models);
nL_mat = nL_mat(:)'; ROI_mat = ROI_mat(:)'; Models_mat = Models_mat(:)';

tnoI = numel(nL_mat); %total number of inputs

%Inputs
kAus = 1./[0.03,0.02,0.01,0.005,0.001];
Amp = repmat(100*1e-9,1,numel(kAus));
Totaldps = round(logspace(1,5,9));
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
    Totaldps(i) = round(Totaldps(i),orders(i));
end
Totaldps = Totaldps(idxdps)
dps_run = min(40,Totaldps);


Aus_way = {'kexp'}; Aus_way_mat = repmat(Aus_way,1,tnoI);
evalMode = {'hpc'}; evalMode_mat = repmat(evalMode,1,tnoI);
vibrField = {'art'}; vibrField_mat = repmat(vibrField,1,tnoI);
sim_method = {'single_sim_totaldps'}; sim_method_mat = repmat(sim_method,1,tnoI);
kAus_con = {'Amp'}; kAus_con_mat = repmat(kAus_con,1,tnoI);

calcStrengthCurve = 1;
maxTime = 23.9*3600;





% Settings
fus = 1e6; fus_mat = repmat(fus,1,tnoI);
meanI = 10; meanI_mat = repmat(meanI,1,tnoI); 
stdI = 3; stdI_mat = repmat(stdI,1,tnoI);  

dpDistribution_mat = repmat({'fibonaccimultiple'},1,tnoI);
dpOrientation = 'radial'; dpOrientation_mat = repmat(cellstr(dpOrientation),1,tnoI);
dpI_space = 'normal'; dpISpace_mat = repmat(cellstr(dpI_space),1,tnoI);
dpI_time = 'artificialPSD_BC'; dpITime_mat = repmat(cellstr(dpI_time),1,tnoI);



for iset=sets%1:tnoI
    Input = struct();   
    
    Settings = {'fus',fus_mat(iset),'meanI',meanI_mat(iset),'stdI',stdI_mat(iset),...
        'dpDistribution',dpDistribution_mat{iset},'dpOrientation',dpOrientation_mat{iset},...
        'dpI_space',dpISpace_mat{iset},'dpI_time',dpITime_mat{iset},...
        'nLayers',nL_mat(iset),'Ifun',str2func(Ifun)};

SettingsStr = ['Settingset: ',num2str(iset)];
Input.Model = Models_mat{iset};
Input.ROI = ROI_mat{iset};
Input.evalMode = evalMode_mat{iset};
Input.sim_method = sim_method_mat{iset};
Input.vibrField = vibrField_mat{iset};
Input.Aus_way = Aus_way_mat{iset};
Input.kAus = kAus;
Input.kAus_con = kAus_con_mat{iset};
Input.Amp = Amp; Input.Totaldps = Totaldps; Input.dps_run = dps_run; 
Input.Settings = Settings; 
Input.SettingsStr = SettingsStr;
Input.calcStrengthCurve = calcStrengthCurve;
Input.maxTime = maxTime;
Input.threshold = threshold;

savenames{iset,1} = sprintf('%s_%s_%s',FileName,suffix,num2str(iset));
save([savenames{iset,1},'.mat'],'Input');


% create table with info
if any(strcmpi('Ifun',Settings))
    Settings{find(strcmpi('Ifun',Settings))+1} = func2str(Settings{find(strcmpi('Ifun',Settings))+1});
end
Set_struct = cell2struct(Settings(2:2:end),Settings(1:2:end-1),2);
intm_struct = catstruct(Input,Set_struct);
Tableinfo_struct(iset) = intm_struct;

end

% Convert cell to a table and use first row as variable names
T = cell2table(savenames,'VariableNames',{'Input'});
 
% Write the table to a CSV file
writetable(T,[FileName,'_',suffix,sprintf('_sets%i-%i',sets(1),sets(end)),'.csv'])

% convert struct to table and save 
Tableinfo = struct2table(Tableinfo_struct);
Tableinfo = movevars(Tableinfo, 'nLayers', 'Before', 'evalMode');
writetable(Tableinfo, [FileName,'_param_',suffix,sprintf('_sets%i-%i',sets(1),sets(end)),'.csv'])




%%
FileName = './Inputs/aIfun/Input_cEG';
Ifun = func2str(@(t) alphaTrainFun(t,'aTInput_220308.mat'));
sets = [1:3:22]; %1:tnoI
idxdps = 7:8;
dpI_time = 'artificialPSD_fixed';
threshold = 0.05;
suffix = sprintf('Amp_fkApdpset_NCfixed_dps%i%i_th%i_nos4l_v2',idxdps(1),idxdps(end),threshold*100);

%look at vibrational interference, Multilayer ({3-}5-10, DOI = Deep,
%Cortex, POI = Scalp, Cortex, Air, Mouse

% variable components
nLayers = [3,5,10];
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
ROIs = {'Deep','Cortex'};
[nL_mat,ROI_mat,Models_mat] = ndgrid(nLayers,ROIs,Models);
nL_mat = nL_mat(:)'; ROI_mat = ROI_mat(:)'; Models_mat = Models_mat(:)';

tnoI = numel(nL_mat); %total number of inputs

%Inputs
intm = load('./Results/SCs_789_cEG_kAus_nos4l_a2piCondBug_th5_aIfun_v2.mat');

kAus_MAT = intm.requiredAmpAll;
%kAus_MAT = kAus_MAT(:,idxdps);
kAus_MAT = repmat(kAus_MAT(:,end)*2,1,length(idxdps));
Amp = logspace(0,5,6).*1e-9;
Totaldps = round(logspace(1,5,9));
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
    Totaldps(i) = round(Totaldps(i),orders(i));
end
Totaldps = Totaldps(idxdps)
dps_run = min(40,Totaldps); %50*500001*5*8 = 953 MB victiny 96GB/36cores = 2.666GB/core


Aus_way = {'kexp'}; Aus_way_mat = repmat(Aus_way,1,tnoI);
evalMode = {'hpc'}; evalMode_mat = repmat(evalMode,1,tnoI);
vibrField = {'art'}; vibrField_mat = repmat(vibrField,1,tnoI);
sim_method = {'single_sim_totaldps'}; sim_method_mat = repmat(sim_method,1,tnoI);
kAus_con = {'Totaldps'}; kAus_con_mat = repmat(kAus_con,1,tnoI);

calcStrengthCurve = 1;
maxTime = 23.9*3600;





% Settings
fus = 1e6; fus_mat = repmat(fus,1,tnoI);
meanI = 10; meanI_mat = repmat(meanI,1,tnoI); 
stdI = 3; stdI_mat = repmat(stdI,1,tnoI);  

dpDistribution_mat = repmat({'fibonaccimultiple'},1,tnoI);
dpOrientation = 'radial'; dpOrientation_mat = repmat(cellstr(dpOrientation),1,tnoI);
dpI_space = 'normal'; dpISpace_mat = repmat(cellstr(dpI_space),1,tnoI);
dpITime_mat = repmat(cellstr(dpI_time),1,tnoI);



for i=1:length(sets)%1:tnoI only with air layer need to be recalc
    iset = sets(i); 
    Input = struct();   
    
    Settings = {'fus',fus_mat(iset),'meanI',meanI_mat(iset),'stdI',stdI_mat(iset),...
        'dpDistribution',dpDistribution_mat{iset},'dpOrientation',dpOrientation_mat{iset},...
        'dpI_space',dpISpace_mat{iset},'dpI_time',dpITime_mat{iset},...
        'nLayers',nL_mat(iset),'Ifun',str2func(Ifun)};

SettingsStr = ['Settingset: ',num2str(iset)];
Input.Model = Models_mat{iset};
Input.ROI = ROI_mat{iset};
Input.evalMode = evalMode_mat{iset};
Input.sim_method = sim_method_mat{iset};
Input.vibrField = vibrField_mat{iset};
Input.Aus_way = Aus_way_mat{iset};
if size(kAus_MAT,1)==tnoI
Input.kAus = kAus_MAT(iset,:);
else
    Input.kAus = kAus_MAT(i,:);
end
Input.kAus_con = kAus_con_mat{iset};
Input.Amp = Amp; Input.Totaldps = Totaldps; Input.dps_run = dps_run; 
Input.Settings = Settings; 
Input.SettingsStr = SettingsStr;
Input.calcStrengthCurve = calcStrengthCurve;
Input.maxTime = maxTime;
Input.threshold = threshold;

savenames{iset,1} = sprintf('%s_%s_%s',FileName,suffix,num2str(iset));
save([savenames{iset,1},'.mat'],'Input');


% create table with info 
Set_struct = cell2struct(Settings(2:2:end),Settings(1:2:end-1),2);
intm_struct = catstruct(Input,Set_struct);
Tableinfo_struct(iset) = intm_struct;

end

% Convert cell to a table and use first row as variable names
T = cell2table(savenames,'VariableNames',{'Input'});
 
% Write the table to a CSV file
writetable(T,[FileName,'_',suffix,sprintf('_sets%i-%i',sets(1),sets(end)),'.csv'])

% convert struct to table and save 
Tableinfo = struct2table(Tableinfo_struct);
Tableinfo = movevars(Tableinfo, 'nLayers', 'Before', 'evalMode');
writetable(Tableinfo, [FileName,'_param_',suffix,sprintf('_sets%i-%i',sets(1),sets(end)),'.csv'])


%% Single Oscillator: vibrField = none; adjusted AMP
FileName = './Inputs/Input_cEG';
sets = [1:3:22];
idxdps = [7,9];
dpI_time = 'artificialPSD_fixed3';
threshold = 0.1;
suffix = sprintf('Amp_fkApdpset_sOsc_NCfixed3_dps%i%i_th%i_nos4l_v3',idxdps(1),idxdps(end),threshold*100);

%look at vibrational interference, Multilayer ({3-}5-10, DOI = Deep,
%Cortex, POI = Scalp, Cortex, Air, Mouse

% variable components
nLayers = [3,5,10];
Models = {'human_scalp','human_cortical','human_air','mouse_fair0'};
ROIs = {'Deep','Cortex'};
[nL_mat,ROI_mat,Models_mat] = ndgrid(nLayers,ROIs,Models);
nL_mat = nL_mat(:)'; ROI_mat = ROI_mat(:)'; Models_mat = Models_mat(:)';

tnoI = numel(nL_mat); %total number of inputs

%Inputs
intm = load('./Results/SCs_all_cEG_kAus_nos4l_a2piCondBug_v2.mat');
requiredAmpAll = load('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Results\SCs_all_cEG_Amp_nos4l_fkApdpset_v2.mat', 'requiredAmpAll');
% narrow down Ampgrid
requiredAmpAll = requiredAmpAll.requiredAmpAll(:,idxdps);
min_rAA = floor(log10(min(requiredAmpAll,[],2)))-1;
max_rAA = ceil(log10(max(requiredAmpAll,[],2)))+1;
if strcmpi(dpI_time,'artificialPSD_fixed3')
    max_rAA = min_rAA;
    min_rAA = min_rAA-3;
elseif strcmpi(dpI_time,'artificialPSD_fixed2')
    max_rAA = max_rAA-1;
    min_rAA = min_rAA-1;
end
Amp = arrayfun(@(a,b)10.^(a:b),min_rAA,max_rAA,'UniformOutput',false);
kAus_MAT = intm.requiredAmpAll;
%kAus_MAT = kAus_MAT(:,idxdps);
kAus_MAT = repmat(kAus_MAT(:,end)*2,1,length(idxdps));
%Amp = logspace(0,5,6).*1e-9;
Totaldps = round(logspace(1,5,9));
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
    Totaldps(i) = round(Totaldps(i),orders(i));
end
Totaldps = Totaldps(idxdps)
dps_run = min(20,Totaldps); %50*500001*5*8 = 953 MB victiny 96GB/36cores = 2.666GB/core


Aus_way = {'kexp'}; Aus_way_mat = repmat(Aus_way,1,tnoI);
evalMode = {'pcwica'}; evalMode_mat = repmat(evalMode,1,tnoI);
vibrField = {'none'}; vibrField_mat = repmat(vibrField,1,tnoI);
sim_method = {'single_sim_totaldps'}; sim_method_mat = repmat(sim_method,1,tnoI);
kAus_con = {'Totaldps'}; kAus_con_mat = repmat(kAus_con,1,tnoI);

calcStrengthCurve = 1;
maxTime = 23.9*3600;




% Settings
fus = 1e6; fus_mat = repmat(fus,1,tnoI);
meanI = 10; meanI_mat = repmat(meanI,1,tnoI); 
stdI = 3; stdI_mat = repmat(stdI,1,tnoI);  

dpDistribution_mat = repmat({'fibonaccimultiple'},1,tnoI);
dpOrientation = 'radial'; dpOrientation_mat = repmat(cellstr(dpOrientation),1,tnoI);
dpI_space = 'normal'; dpISpace_mat = repmat(cellstr(dpI_space),1,tnoI);
dpITime_mat = repmat(cellstr(dpI_time),1,tnoI);



for iset=sets
    Input = struct();   
    
    Settings = {'fus',fus_mat(iset),'meanI',meanI_mat(iset),'stdI',stdI_mat(iset),...
        'dpDistribution',dpDistribution_mat{iset},'dpOrientation',dpOrientation_mat{iset},...
        'dpI_space',dpISpace_mat{iset},'dpI_time',dpITime_mat{iset},...
        'nLayers',nL_mat(iset)};

SettingsStr = ['Settingset: ',num2str(iset)];
Input.Model = Models_mat{iset};
Input.ROI = ROI_mat{iset};
Input.evalMode = evalMode_mat{iset};
Input.sim_method = sim_method_mat{iset};
Input.vibrField = vibrField_mat{iset};
Input.Aus_way = Aus_way_mat{iset};
Input.kAus = kAus_MAT(iset,:);
Input.kAus_con = kAus_con_mat{iset};
Input.Amp = Amp{iset}; Input.Totaldps = Totaldps; Input.dps_run = dps_run; 
Input.Settings = Settings; 
Input.SettingsStr = SettingsStr;
Input.calcStrengthCurve = calcStrengthCurve;
Input.maxTime = maxTime;
Input.threshold = threshold;

savenames{iset,1} = sprintf('%s_%s_%s',FileName,suffix,num2str(iset));
save([savenames{iset,1},'.mat'],'Input');


% create table with info 
Set_struct = cell2struct(Settings(2:2:end),Settings(1:2:end-1),2);
intm_struct = catstruct(Input,Set_struct);
Tableinfo_struct(iset) = intm_struct;

end

% Convert cell to a table and use first row as variable names
T = cell2table(savenames,'VariableNames',{'Input'});
 
% Write the table to a CSV file
writetable(T,[FileName,'_',suffix,sprintf('_sets%i-%i',sets(1),sets(end)),'.csv'])

% convert struct to table and save 
Tableinfo = struct2table(Tableinfo_struct);
Tableinfo = movevars(Tableinfo, 'nLayers', 'Before', 'evalMode');
writetable(Tableinfo, [FileName,'_param_',suffix,sprintf('_sets%i-%i',sets(1),sets(end)),'.csv'])



