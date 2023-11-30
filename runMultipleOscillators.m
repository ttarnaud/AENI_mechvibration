% Run Multiple oscillators
close all
clear all
clc
addpath(genpath('./Functions'));
addpath(genpath('./DataBases'));
addpath('D:\users\rschoeters\Documents\Imec USEEG\Sim4life')
try
s = load('randomseed.mat');
s = s.s;
catch
    warning('new random seed')
    s = rng;
end
model = 'sHHM_11s';
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

%% test normal
rng(s)
Input_iBN = {'posDp',posDp,'POIs',[0,0,1;0,1,0;1,0,0;1/sqrt(3),1/sqrt(3),1/sqrt(3)],'OrienDipole',orienDipole,...
    'Plot',0,'totaldps',10,'dps_run',10,'showSphere','','Solutiontype',SolutionType,...
    'dpDistribution',dpDistri,'Sim4life',Input};
[RMSpc,SNRpc,Time,OUT] = investBiologicalNoise(horzcat(Input_iBN,{'eval_method','normal'}));
%  [RMSpc,SNRpc] = investBiologicalNoise('posDp',[0.001,0.001,0.001;0,0,0;0.01,0.01,0.01],'POIs',[1,0,0;0,1,0;0,0,1],'OrienDipole',[-1,0,0],...
%     'Plot',1,'totaldps',10,'dps_run',10,'showSphere','','Solutiontype',SolutionType,...
%     'Dirus',[1,1,1],'Phaseus',[pi,pi/2,pi/3],'tend',100e-3);
%%  Plot metrics
%metric fields can be changed. check fieldnames
Y=[[RMSpc.RMS]',[RMSpc.RMSDOI]'];
%Y=[[RMSpc.RMSnorm]'];
figure
subplot(1,3,1)
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
legend('Normalized signal','wrt reconstruction of oscillating dipole only')
title('RMS')
subplot(1,3,2)
Y=[[RMSpc.Q2]'];
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
title('Q^2')

Y = [[SNRpc.DOI_noiseall(2,:)]',[SNRpc.DOI_noiseosc(2,:)]',[SNRpc.DOI_noisestat(2,:)]'];

subplot(1,3,3)
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
legend('SNR_{all}','SNR_{osc}','SNR_{stat}')
title('SNR')
pause(1)
%% test database
rng(s)
[RMSpc_DB,SNRpc_DB,Time_DB,OUT_DB] = investBiologicalNoise(horzcat(Input_iBN,{'eval_method','database'}));
%%
Y=[[RMSpc_DB.RMS]',[RMSpc_DB.RMSDOI]'];
%Y=[[RMSpc.RMSnorm]'];
figure
subplot(1,3,1)
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
legend('Normalized signal','wrt reconstruction of oscillating dipole only')
title('RMS')
subplot(1,3,2)
Y=[[RMSpc_DB.Q2]'];
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
title('Q^2')

Y = [[SNRpc_DB.DOI_noiseall(2,:)]',[SNRpc_DB.DOI_noiseosc(2,:)]',[SNRpc_DB.DOI_noisestat(2,:)]'];

subplot(1,3,3)
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
legend('SNR_{all}','SNR_{osc}','SNR_{stat}')
title('SNR')

%% test database_f
rng(s)
[RMSpc_DB_f,SNRpc_DB_f,Time_DB_f,OUT_DB_f] = investBiologicalNoise(horzcat(Input_iBN,{'eval_method','database_f'}));
%%
Y=[[RMSpc_DB_f.RMS]',[RMSpc_DB_f.RMSDOI]'];
%Y=[[RMSpc.RMSnorm]'];
figure
subplot(1,3,1)
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
legend('Normalized signal','wrt reconstruction of oscillating dipole only')
title('RMS')
subplot(1,3,2)
Y=[[RMSpc_DB_f.Q2]'];
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
title('Q^2')

Y = [[SNRpc_DB_f.DOI_noiseall(2,:)]',[SNRpc_DB_f.DOI_noiseosc(2,:)]',[SNRpc_DB_f.DOI_noisestat(2,:)]'];

subplot(1,3,3)
bar(Y)
set(gca,'xticklabel',{'POI_1','POI_2','POI_3','POI_4','mean all','mean dominat'})
legend('SNR_{all}','SNR_{osc}','SNR_{stat}')
title('SNR')

rmpath(genpath('./Functions'));
rmpath(genpath('./DataBases'));