%processing the results of full pressure fields on EM simulations
clear all
close all

Old_dir = cd('D:\no backup\EEGUS\ResultsHPC');
%list of all mat files in current folder
lall = dir('*.mat');
cd(Old_dir);
l = lall([lall.datenum]>737835 & [lall.bytes]>1e5); % only results after 14/02 and large enough file (to small probably error)
idx = find(contains({l(:).name},'Settingset1_'));

load(fullfile(l(idx).folder,l(idx).name));

Amp = Input.Amp;
Totaldps = Input.Totaldps;
lAmp  = length(Amp);
lTdps = length(Totaldps);

Totaldps_plot = 10;
PLOT_flag = 0;
PLOT_s4l = 0;
Display = 0;
%% reobtain data
dps_run = Input.dps_run; 
Settings = Input.Settings;  SettingsStr = Input.SettingsStr;
Fignr = Input.Fignr;
HPC_flag = Input.HPC_flag;
S4l_flag = Input.S4l_flag; S4l_fn = Input.Filename_S4l; S4l_loc = Input.Location_S4l;
dpDistribution = Settings{find(strcmpi(Settings,'dpDistribution'))+1};
dpOrientation = Settings{find(strcmpi(Settings,'dpOrientation'))+1};
sigma = 0.33;
SphereRes = 30;
d = 500*10^-6; %distance between current source and current sink
nLayers = 10;

Model = Input.Model;
switch lower(Model)
    case 'human_scalp'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.082;
        RBrain = 0.07;
        dRBrain = 0.005;
    case 'human_cortical'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.07;
        RBrain = 0.07;
        dRBrain = 0.005;
    case 'mouse_fair0'
        SolutionType = 'Mouse4Sphere~fair0';
        RPOI = 0.0059;
        RBrain = 0.0046;
        dRBrain = 0.0005;
    otherwise 
        error('false input model')        
end  

[Options,RSphere] = getSettings(SolutionType,1);
if strcmpi(SolutionType,'Mouse4Sphere~fair0') || strcmpi(SolutionType,'Mouse4Sphere~fair1')
dRSphere = 0.0005;
else
dRSphere = 0.005;
end


[X,Y] = meshgrid(Totaldps,Amp);
Y = Y*1e9;
% Load pressure field
S4l_input = struct();
[S4l_input.veloc, S4l_input.press, S4l_input.loc, S4l_input.fus, maxv,maxv_pos,maxv_dir] = ...
    calcvelocity(S4l_loc,S4l_fn,RBrain-dRBrain,1,'snapshot',0,'rotate',0','totalperiods',50,'inittime',48);
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
POIs = [1,0,0;0,0,1;0,1,0;maxv_pos+ts(idxt)*maxv_dir;1,1,1];
POIs = POIs./vecnorm(POIs,2,2);
POIs = RPOI.*POIs;
Input_AdaptVamp.flag = 1;
Input_AdaptVamp.vmax = maxv;
posDp = maxv_pos;
OrienDipole = maxv_dir;
[CSource,CSink] = GendpPos(dpDistribution,dpOrientation,d,RSphere-dRSphere,Totaldps_plot-size(posDp,1),Display,'nLayers',nLayers);
CSource = vertcat(posDp+d/2*OrienDipole,CSource);
CSink = vertcat(posDp-d/2*OrienDipole,CSink);
%% start plotting
varargin2 = horzcat(Options,{'POI',POIs,'DOI',1,'OSC',1:Totaldps_plot,'scale',0});
PotentialSphere_Multi(CSource,CSink,10,sigma,RSphere,SphereRes,1,varargin2);
disp('succes')

for iAmp = 1:lAmp
    for idps = 1:lTdps
        
    break
    
    end    
end
