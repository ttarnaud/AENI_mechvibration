%% run multi comparisons multiply oscillators for efficiency database
% nnot finished probably
% !!!!!!!!!!!!IN CONSTRUCTION!!!!!!!!!!!!!!!!!!!!

if HPC_flag
    try
        c = parcluster('local');
        pool = parpool(c.NumWorkers);
    end
    S4l_loc = '';
else
close all
clear all
clc
addpath(genpath('./Functions'));
addpath(genpath('./DataBases'));
addpath('D:\users\rschoeters\Documents\Imec USEEG\Sim4life')
HPC_flag = 0;
end
try
s = load('randomseed.mat');
s = s.s;
catch
    warning('new random seed')
    s = rng;
end

Models = {'sHHM','sHHM_cortex'};

for imodel = 1:length(Models)
model = 'sHHM';
switch model
    case 'sHHM_old'
        location = 'D:\no backup\Sim4Life\3spherical\CenterFocus\SimplifiedHumanHeadModel_3spheres.smash_Results\';
        Filename = 'ExportedDataOF_sHHM180_focusedSEFT_1MHz_aw80_gs075mm';
        Input.RBrain = 0.07;
        SolutionType = '3SphereS8.2R25';
        orienDipole = [0,1,0];
        posDp = [0,0,0];
        dpDistri = 'fibonaccimultiple';
        plane = 'XY';
    case 'sHHM'
        location = 'D:\no backup\Sim4Life\3spherical\CenterFocus\Wicasim3_sHHM\';
        Filename = 'ExportedData_SHHM_AW164_gs0375mm';
        Input.RBrain = 0.07;
        SolutionType = '3SphereS8.2R25';
        orienDipole = [];
        posDp = [];
        dpDistri = 's4l_hotspot';
        plane = 'YZ';
    case 'sHHM_11s'
        location = 'D:\no backup\Sim4Life\3spherical\CenterFocus\Wicasim3_sHHM\';
        Filename = 'ExportedData_SHHM_AW164_cropx_11slices';
        Input.RBrain = 0.07;
        SolutionType = '3SphereS8.2R25';
        
        posDp = [0,0,0;0,0.010,0.010;0,-0.010,-0.05];
        orienDipole = [0,0,1;0,0.010,0.010;0,-0.010,-0.05];
        orienDipole = orienDipole./vecnorm(orienDipole,2,2);
        dpDistri = 'fibonaccisingle';
        plane = 'YZ';
    case 'RRM180'
        location = 'D:\no backup\Sim4Life\3spherical\CenterFocus\3sphericalMouse.smash_Results\';
        Filename = 'ExportedDataRR_M180_focusedSEFT_1MHz_aw10';
        Input.RBrain=0.0046;
        SolutionType = 'Mouse4Sphere~fair0';
        orienDipole = [1,0,0];
        posDp = [0,0,0];
        dpDistri = 'fibonaccimultiple';
        plane = 'XY';
end

Plot_flag = 1;
[veloc, press, loc, fus] = calcvelocity(location,Filename,Input.RBrain,1,...
    'Database',1,'Plane',plane,'rotate',0','totalperiods',50,'inittime',48);
Input.veloc = veloc; Input.fus = fus; Input.loc = loc; Input.loc.unit = 'mm';
Input.S4l_flag = 1;
Input.resolution=0.005;


nr_dps = logspace(0,5,11);
dps_run = nr_dps;
dps_run(dps_run>100) = 100;
Distributions = {'fibonaccisingle','fibonaccimultiple','s4l_hotspot'};
for idp_distri = 1:length(dp_distri)
    for inr_dp = 1:length(nr_dps)
        %% test normal
        rng(s)
        Input_iBN = {'posDp',posDp,'POIs',[0,0,1;0,1,0;1,0,0;1/sqrt(3),1/sqrt(3),1/sqrt(3)],'OrienDipole',orienDipole,...
            'Plot',0,'totaldps',nr_dps(inr_dp),'dps_run',dps_run(inr_dp),'showSphere','','Solutiontype',SolutionType,...
            'dpDistribution',Distributions{idp_distri},'Sim4life',Input};
        [RMSpc,SNRpc,Time,OUT] = investBiologicalNoise(horzcat(Input_iBN,{'eval_method','normal'}));
        
        rng(s)
        [RMSpc_DB,SNRpc_DB,Time_DB,OUT_DB] = investBiologicalNoise(horzcat(Input_iBN,{'eval_method','database'}));
    end
end
end