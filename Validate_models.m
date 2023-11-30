% test if same output is given
close all
clear all
clc
addpath(genpath('./Functions'));
addpath(genpath('./Inputs'));
addpath('D:\users\rschoeters\Documents\Imec USEEG\Sim4life')
model = 'RRM180';
SolutionTypes = {'3SphereS8.2R25','3SphereS8.2R~f','4SphereS8.7R25',...
    '4SphereS8.7R~f','Mouse4Sphere~fair0','Mouse4Sphere~fair1'};
ParallelCompute=[0,1]



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
    case 'RRM180'
        location = 'D:\no backup\Sim4Life\3spherical\CenterFocus\3sphericalMouse.smash_Results\';
        Filename = 'ExportedDataRR_M180_focusedSEFT_1MHz_aw10';
        %Input.RBrain=0.0046;
        Input.RBrain=0.07;
        SolutionType = 'Mouse4Sphere~fair0';
        orienDipole = [1,0,0];
        posDp = [0,0,0];
        dpDistri = 'fibonaccimultiple';
        plane = 'XY';
end

Plot_flag = 1;
[veloc, press, loc, fus] = calcvelocity(location,Filename,1,'Plane',plane,'rotate',0','totalperiods',50,'inittime',48);
Input.veloc = veloc; Input.fus = fus; Input.loc = loc; Input.loc.unit = 'mm';
Input.S4l_flag = 1;
Input.resolution=0.005;

%% test
s=rng
for isoltype = 6
rng(s)
[RMS,SNR] = investBiologicalNoise('posDp',posDp,'POIs',[0,0,1;0,1,0;1,0,0;1/sqrt(3),1/sqrt(3),1/sqrt(3)],'OrienDipole',orienDipole,...
    'Plot',0,'totaldps',3,'dps_run',3,'showSphere','wpoi','Solutiontype',SolutionTypes{isoltype},...
    'dpDistribution',dpDistri,'Sim4life',Input,'Parallelcompute',ParallelCompute(1));
rng(s)
[RMSpc,SNRpc] = investBiologicalNoise('posDp',posDp,'POIs',[0,0,1;0,1,0;1,0,0;1/sqrt(3),1/sqrt(3),1/sqrt(3)],'OrienDipole',orienDipole,...
    'Plot',0,'totaldps',3,'dps_run',3,'showSphere','wpoi','Solutiontype',SolutionTypes{isoltype},...
    'dpDistribution',dpDistri,'Sim4life',Input,'Parallelcompute',ParallelCompute(2));

 1
 close all
 
    
end
