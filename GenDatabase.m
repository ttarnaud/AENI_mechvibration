% A database/Lookup table is generated for a certain model(MODEL). the
% lookup table contains dp_pos along z-axis and N uniformly spread vibration directions.
%the signal is measured at defined POIs. Use simDipole_DB_f to calculate
%potential of given problem with database option
%Current Lookup table method has however still flaws, only linear
%polarisation and still significant errors: see raportations 
addpath(genpath('./Functions'))
Model = 'human_scalp';
fus = 1e6;

% Select model
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

% dipole positions
dp_pos = [0,0,1].*[0:0.5:65]'*1e-3;

% POIs
dtheta = 0.5/65;
theta = 0:dtheta:pi;
POIs = [zeros(size(theta')),sin(theta'),cos(theta')]*RPOI;

% Vibrations
N = 2000;
indices = [0:N-1]+0.5;
theta_vdir = acos(1-2*indices'/N);
phi_vdir = pi*(1+5^(1/2))*indices';
vdir = [cos(phi_vdir).*sin(theta_vdir),sin(phi_vdir).*sin(theta_vdir),cos(theta_vdir)];




Input = {'OrienDipole',[0,0,1],'Totaldps',1,'dps_run',1,'Aus',1e-6,...
    'ParallelCompute',0,'SolutionType',SolutionType,'POIs',POIs,...
    'scale',0,'PLOT',0,'genDatabase',1,'meanI',1,'fus',fus,'Tend',1/fus,'Ifun',@(t) 1+0*t};

Pflag = 0;
for ipos=1:size(dp_pos,1)
    DataBase = zeros(size(vdir,1),size(POIs,1),21);
    selpos = dp_pos(ipos,:);
    parfor ivib = 1:size(vdir,1)
        [RMS,SNR,TSTOP,Out] = investBiologicalNoise(...
            horzcat({'posDp',selpos,'Dirus',vdir(ivib,:)},Input));
        DataBase(ivib,:,:) = Out.VR*1e-6;% scale back to V       
    end
    Pval = ipos/size(dp_pos,1)*100;
    if Pval>Pflag
        Pflag = Pflag + 5;
        
        fprintf('%s\nComplete: %6.2f%%\n',datestr(now),Pval)
    end
    Filename = sprintf('Database%s_%0.2fMHz_pos%0.2fmm.mat',Model,fus/1e6,selpos(1,3)*1000);    
    save(['.\DataBases\DataBases_v3\',Filename],'DataBase','theta_vdir','phi_vdir','POIs','Model','fus','selpos')
end

%% adjust Database

% adjustement of first generated Database: rescaling was necessary. Double
% chekc if already corrected in part above
 Model = 'human_scalp';
fus = 1e6;
dp_pos = 0:0.5:65;
wus = 2*pi*fus;
dt = 1/fus/20;
load_location = '.\DataBases\DataBases_v1\';
save_location = '.\DataBases\DataBases_v2\';
for idp = 1:length(dp_pos)

    Filename = sprintf('Database%s_%0.2fMHz_pos%0.2fmm.mat',Model,fus/1e6,dp_pos(idp));
    Database = load([load_location,Filename]);
    DataBase = Database.DataBase*1e-7; %correct for µV and 10µA;
    
    statval = mean(DataBase(:,:,1:end-1),3);
    amp = (max(DataBase,[],3)-min(DataBase,[],3))./2;
    
    shiftdatabase = round(DataBase-statval,16);
    A1 = squeeze(shiftdatabase(:,:,1));
    A2 = squeeze(shiftdatabase(:,:,2));
    A1 = A1+eps(A2);
    num = sin(wus*dt);
    denum = A2./A1-cos(wus*dt);
    phase = round(atan2(num,denum),6);
    
    phi_vdir = Database.phi_vdir;
    theta_vdir = Database.theta_vdir;
    selpos = Database.selpos;
    POIs = Database.POIs;
    
    save([save_location,Filename],'DataBase','Model','fus','statval','amp','phase','phi_vdir','theta_vdir','selpos','POIs')
    pause(0.1);
    
    
end

%% add v2

   projectdir = 'D:\users\rschoeters\Documents\Imec USEEG\Matlab\DataBases\DataBases_v2';
   dinfo = dir( fullfile(projectdir, '*.mat') );
   oldnames = {dinfo.name};
   unwanted = cellfun(@isempty, regexp(oldnames, 'mat$') );
   oldnames(unwanted) = [];   
   oldnames_intm = regexprep(oldnames, '_v2.','.');
   newnames = regexprep(oldnames_intm, '_scalp',  '$0_v2');
   for K = 1 : length(oldnames)
      movefile( fullfile(projectdir, oldnames{K}), fullfile(projectdir, newnames{K}) );
   end