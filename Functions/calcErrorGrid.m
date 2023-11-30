function [RMSall,SNRall,Outall,TSTOPall,Totaldps,requiredAmp,requiredAmplower,itermat,cSC_mode] = calcErrorGrid(Input)
%Calculate StrengthCurve vs #dipoles with multiple oscilatores
%Output is saved when HPC mode
startTimer = tic;
disp('start function')
Input_str = '';        %add name of input to savename
    
Threshold = [0.05];
vibrField = 'none';     % vibrational field distribution 'none','Art',S4L
calcStrengthCurve = 0;

pp_flag = 1;

% settings strength curve
flaglim = 1e-10;
maxIter = 10;%5;
maxTime = 24*3600;
maxTime_overall = maxTime;

S4l_flag = 0;
HPC_flag = 0;
lim_numWorkers_flag = 0;

sim_method = 'seq';     %'seq', 'single_sim_totaldps', 'validate'
evalMode = 'hpc';        %'hpc','local','debug'
PLOT_flag = 0;          %plot results
PLOT_s4l = 0;
Outall = struct();
RMSall = struct();
SNRall = struct();

requiredAmp = [];
requiredAmplower = [];
itermat = [];
cSC_mode = 'amp';
disp('start loading input')


if ischar(Input)
    Input_str = Input;
    fprintf('\nInput: %s\n',Input)
    try
    Input_temp = load(['./Inputs/',Input,'.mat']);
    catch
        Input_temp = load([Input,'.mat']);
    end
    Input = Input_temp.Input;
end
  

Amp = Input.Amp; 
Totaldps = Input.Totaldps; dps_run = Input.dps_run; 
Settings = Input.Settings;  SettingsStr = Input.SettingsStr;
evalMode = Input.evalMode;  vibrField = Input.vibrField;
ROI = Input.ROI; Model = Input.Model; 

if isfield(Input,'calcStrengthCurve')
    calcStrengthCurve = Input.calcStrengthCurve;
end
if isfield(Input,'sim_method')
    sim_method = Input.sim_method;
end
if isfield(Input,'maxTime')
    maxTime = Input.maxTime;
end
if isfield(Input,'pp_flag')
    pp_flag = Input.pp_flag;
end
if isfield(Input,'threshold')
    Threshold = Input.threshold;
end

if strcmpi(vibrField,'s4l')    
    %Field according to S4L
    S4l_flag = 1;
    S4l_fn = Input.Filename_S4l; S4l_loc = Input.Location_S4l;   
    
end




% Change of settings based on Debug of local run flag
switch lower(evalMode)
    case 'debug'
        disp('DEBUG modus: limited number of dipoles !!!!')
        disp('')
        %idx_dpitime = find(strcmpi(Settings,'dpI_time'))+1;
        %Settings{idx_dpitime} = 'artificialPSD_BC_val';
        %idx_dpispace = find(strcmpi(Settings,'dpI_space'))+1;
        %Settings{idx_dpispace} = 'val';
        %maxIter = 5;
        %idxkAus = length(Totaldps)
        Totaldps = [40,80]%,10000]%min(20,Totaldps);
        dps_run = [40,40]%,100]%min(20,Totaldps);
        %kAus_end = length(Amp);
        %Amp = Amp(1:2:end);
        
        %dps_run = min(10,Totaldps);
        %Amp_Input = Amp_Input([1,end]);
        S4l_flag = 0;
        PLOT_flag = 0;
        PLOT_s4l = 0;
        save_flag = 1;
        showSphere = '';
        calcStrengthCurve = 1;
        %Amp = 1e-6;
        
    case 'local'
        disp('LOCAL RUN!!!!!')
        disp('___________________________')
        Amp = 10e-6;%logspace(2,5,4).*1e-9;
        Totaldps = 10%round(logspace(1,4,4));%round(logspace(0,5+double(isettings>=9 & isettings<=12),6+2*+double(isettings>=9 & isettings<=12)));
        orders = -floor(log10(Totaldps));
        for i=1:length(Totaldps)
            Totaldps(i) = round(Totaldps(i),orders(i));
        end
        dps_run = min(1e3,Totaldps);
        if S4l_flag
            S4l_fn = 'ExportedData_sHHM_1SEFT_fp0_55_0_cropx19s';
        end
        save_flag = 0;
        showSphere = '';
        
    case 'hpc'
        HPC_flag = 1;
        save_flag = 1;
        if isfield(Input,'maxTime')
            maxTime_overall = Input.maxTime;
        else
            maxTime_overall = 72*3600;
        end
        showSphere = '';
    case 'pcwica'
        
        if isdeployed
            lim_numWorkers_flag = 1;
            NumWorkers = 3;
        end
        HPC_flag = 1;
        save_flag = 1;
        if isfield(Input,'maxTime')
            maxTime_overall = Input.maxTime;
        else
            maxTime_overall = 72*3600;
        end
        showSphere = '';
            
    case 'normal'
        save_flag = 0;
        showSphere = '';
end

% which model to use: 3 predefined models (more possible but this way
% better for consistency) %change made here also adjust in Process_mosc/createfig
switch lower(Model)
    case 'human_scalp'
        SolutionType = '4SphereS8.7R25'; %could use 3spherical model but for comparative reasons not applied
        RPOI = 0.082;
        RBrain = 0.07;
        dRBrain = 0.005;
    case 'human_cortical'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.07;
        RBrain = 0.07;
        dRBrain = 0.005;
    case 'human_air'
        SolutionType = '4sphereS8.7R25';
        RPOI = 0.087;
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

%region of interest deep (centre or at cortex)
switch lower(ROI)
    case 'deep'
        adj_drSphere = 0;
        if any(strcmpi(vibrField,{'none','Art'}))
            Input_dppos = [0,0,0];
            Input_OrienDipole = [0,0,1];         
        end
    case 'cortex'
        adj_drSphere = 1;
        if any(strcmpi(vibrField,{'none','Art'}))
            Input_dppos = [0,0,RBrain-dRBrain];
            Input_OrienDipole = [0,0,1];
        end
    otherwise
        error('false input')
end

%if HPC_flag or isdeployed load parallelpool
if HPC_flag || isdeployed
    try
        c = parcluster('local');
        if lim_numWorkers_flag
            pool = parpool(NumWorkers);
        else
            pool = parpool(c.NumWorkers);
        end
    end
    if S4l_flag
        S4l_loc = './PressureFields/';
    end
elseif ~isdeployed
    addpath(genpath('./Functions'));
    %addpath('D:\users\rschoeters\Documents\Imec USEEG\Sim4life')
end

disp(['Startingmodel: ',Model]);
%fix rng for consistency
try
s = load('randomseed.mat');
s = s.s;
catch
    warning('new random seed')
    s = rng;
    save('randomseed.mat','s');
end



switch lower(vibrField)
    case 's4l'
        % load pressure field en covert to vibration field
        S4l_input = struct();
        [S4l_input.veloc, S4l_input.press, S4l_input.loc, S4l_input.fus, maxv,maxv_pos,maxv_dir] = ...
            calcvelocity(S4l_loc,S4l_fn,RBrain-dRBrain,PLOT_s4l,'snapshot',0,'rotate',0');
        S4l_input.loc.unit = 'mm';
        if strcmpi(S4l_input.loc.unit,'mm')
            maxv_pos = maxv_pos*1e-3;
        end
        S4l_input.resolution=Input.S4l_res; %resolution used in case of hotspot dipole placement
        S4l_input.RBrain = RBrain; 
        % place POI in direction of vibration
        t1 = fzero(@(t) vecnorm(maxv_pos+t'*maxv_dir,2,2)-RBrain,[0,1]);
        t2 = fzero(@(t) vecnorm(maxv_pos+t'*maxv_dir,2,2)-RBrain,[-1,0]);
        ts = [t1,t2];
        [~,idxt] = min(abs(ts));
        POIs = [1,0,0;0,0,1;0,1,0;maxv_pos+ts(idxt)*maxv_dir;1,1,1];
        POIs = POIs./vecnorm(POIs,2,2);
        POIs = RPOI.*POIs;                  %POIs at correct radius
        Input_AdaptVamp.flag = 1;           %Adapt vibrational amplitude for sweep
        Input_AdaptVamp.vmax = maxv;
        posDp = maxv_pos;                   % diple of interest =  biggest vibrator
        OrienDipole = maxv_dir;
        Dirus_way = 'input';                %S4L flag is on will not be evaluated 
        Dirus = maxv_dir;                   %orientation of idpole and vibration are alligned => best result
        Aus_way = 'normal';                 %S4L flag is on will not be evaluated 
        %savename_add = '_s4l';
        kAus = [];
    case 'art'        
        Input_AdaptVamp.flag = 0;           %normally not evaluated
        posDp = Input_dppos;
        OrienDipole = Input_OrienDipole;    %dipole oreintation and direction are the same
        Dirus = OrienDipole;
        Dirus_way = 'input';               % vibration of all other dipoles is radial
        Aus_way = Input.Aus_way;
        kAus = Input.kAus;
        kAus_con = Input.kAus_con;          % to which parameter is kAus connected: none = single kAus input
                                            % Totaldps: kAus is fixed for given nr of dipoles: kAus needs same length as Totaldps
                                            % Amp: kAus is connected to Amp,kAus same length as Amp   
                                           
        if exist('idxkAus','var')
            kAus = kAus(idxkAus);
        end
        %POIs % place POI as projecton of DOI in direction of vibration on
        %Radius of interest
        t1 = fzero(@(t) vecnorm(posDp+t'*OrienDipole,2,2)-RBrain,[0,1]);
        t2 = fzero(@(t) vecnorm(posDp+t'*OrienDipole,2,2)-RBrain,[-1,0]);
        ts = [t1,t2];
        [~,idxt] = min(abs(ts));
        POIs = [1,0,0;0,0,1;0,1,0;posDp+ts(idxt)*OrienDipole;1,1,1];
        POIs = POIs./vecnorm(POIs,2,2);
        POIs = RPOI.*POIs;
        % all dipoles are oscillating strength based on Aus_way and kAus
        all_osc = 1;
        Settings = horzcat(Settings,{'all_osc',all_osc});
        %savename_add = '_nos4l';
    case 'none'
        Input_AdaptVamp.flag = 0;           %normally not evaluated
        posDp = Input_dppos;
        OrienDipole = Input_OrienDipole;    %dipole oreintation and direction are the same
        Dirus = OrienDipole;
        Dirus_way = 'input';               % vibration of all other dipoles is radial
        Aus_way = 'normal';
        %POIs % place POI as projecton of DOI in direction of vibration on
        %Radius of interest
        t1 = fzero(@(t) vecnorm(posDp+t'*OrienDipole,2,2)-RBrain,[0,1]);
        t2 = fzero(@(t) vecnorm(posDp+t'*OrienDipole,2,2)-RBrain,[-1,0]);
        ts = [t1,t2];
        [~,idxt] = min(abs(ts));
        POIs = [1,0,0;0,0,1;0,1,0;posDp+ts(idxt)*OrienDipole;1,1,1];
        POIs = POIs./vecnorm(POIs,2,2);
        POIs = RPOI.*POIs;
        %savename_add = '';
        kAus = [];
        idx_kAus = [];
        kAus_con = [];
        
end
S4l_input.S4l_flag = S4l_flag;
%add variables so saved with them
Input.Dirus_way = Dirus_way;
Input.OrienDipole = OrienDipole;
Input.Dirus = Dirus;
Input.posDp = posDp;
Input.POIs = POIs;


%% start simulations
switch lower(sim_method)
    
    case 'seq'
        disp('start sim_method  seq')
        for iAmp = 1:length(Amp)
            for idps = 1:length(Totaldps)
                rng(s)
                fprintf('\n%s\n Amp: %i/%i, idps: %i/%i\n',datestr(now),iAmp,length(Amp),idps,length(Totaldps))
                Amp_input = Amp(iAmp);
                switch lower(vibrField)
                    case 's4l'
                        Input_AdaptVamp.Aus = Amp_input;
                        k_Aus = [];
                    case 'art'
                        switch lower(kAus_con)
                            case 'none'
                                k_Aus = kAus(1);
                            case 'amp'                                
                                k_Aus = kAus(iAmp);
                            case 'totaldps'
                                k_Aus = kAus(idps);
                            otherwise
                                error('false kAus_con input')
                        end
                    case 'none'
                        k_Aus = [];
                end
                input_set = horzcat({'posDp',posDp,'OrienDipole',OrienDipole,...
                    'Totaldps',Totaldps(idps),'dps_run',dps_run(idps),'Aus',Amp_input,'Aus_way',Aus_way,...
                    'k_Aus',k_Aus,'Dirus',Dirus,'Dirus_way',Dirus_way,'ParallelCompute',1,'SolutionType',SolutionType,...
                    'POIs',POIs,'scale',0,'showSphere',showSphere,'PLOT',PLOT_flag,'adj_drSphere',adj_drSphere,'HPC_flag',HPC_flag,...
                    'Adapt_vAmp',Input_AdaptVamp,'Sim4life',S4l_input,'randomseed',s},Settings);
                [RMS,SNR,TSTOP,Out] = investBiologicalNoise(input_set{:});
                RMSall(iAmp,idps).RMS = RMS;
                SNRall(iAmp,idps).SNR = SNR;
                Outall(iAmp,idps).Out = Out;
                TSTOPall(iAmp,idps).TSTOP = TSTOP;
            end
        end
        timeGrid = toc(startTimer);
        if calcStrengthCurve            
            for iov = 1:length(Amp)
                for idps = 1:length(Totaldps)
                    for iPOIs = 1:length(RMSall(iov,idps).RMS)
                        % Decompose RMS structure
                        RMS_mat(iov,idps,iPOIs) = RMSall(iov,idps).RMS(iPOIs).RMS;
                    end
                end
            end
            
            input_set_init = horzcat({'posDp',posDp,'OrienDipole',OrienDipole,...
                'Aus_way',Aus_way,'Dirus',Dirus,'Dirus_way',Dirus_way,'ParallelCompute',1,'SolutionType',SolutionType,...
                'POIs',POIs,'scale',0,'PLOT',PLOT_flag,'adj_drSphere',adj_drSphere,'HPC_flag',HPC_flag,...
                'Adapt_vAmp',Input_AdaptVamp,'Sim4life',S4l_input,'randomseed',s},Settings);
            [requiredAmp,requiredAmplower,itermat,cSC_mode] = cSC_cEG(Input,...
                Amp,kAus,kAus_con,vibrField,flaglim,evalMode,maxTime,maxTime_overall,timeGrid,...
                Totaldps,dps_run,RMS_mat,Threshold,maxIter,input_set_init,SettingsStr,Input_str,s);
        end
    case 'single_sim_totaldps'        
         disp('start sim_method  single_sim_totaldps')
        for idps = 1:length(Totaldps)
            % run in normal case
            startTimer_dprun = tic;
            %scale runtime to number of dipoles and dtime left, but at
            %least 30 min if possible (min(max,...))
            maxTime_dprun = floor((maxTime-toc(startTimer)).*Totaldps(idps)./sum(Totaldps(idps:end)));
            maxTime_dprun = min(max(maxTime_dprun,30*60),(maxTime-toc(startTimer)));
            s_intm = seconds(maxTime_dprun); s_intm.Format = 'hh:mm:ss.SSS';            
            iAmp = ceil(length(Amp)/2);
            rng(s)
            fprintf('\n%s\n Amp: %i/%i, idps: %i/%i\n',datestr(now),iAmp,length(Amp),idps,length(Totaldps))
            fprintf('Time this run: %s\n',s_intm)
            Amp_input = Amp(iAmp);
            switch lower(vibrField)
                case 's4l'
                    Input_AdaptVamp.Aus = Amp_input;
                    k_Aus = [];
                case 'art'
                    switch lower(kAus_con)
                        case 'none'
                            k_Aus = kAus(1);
                        case 'amp'
                            k_Aus = kAus(iAmp);
                        case 'totaldps'
                            k_Aus = kAus(idps);
                        otherwise
                            error('false kAus_con input')
                    end
                case 'none'
                    k_Aus = [];
            end
            input_set = horzcat({'coll_singlePeriod_flag',1,'posDp',posDp,'OrienDipole',OrienDipole,...
                'Totaldps',Totaldps(idps),'dps_run',dps_run(idps),'Aus',Amp_input,'Aus_way',Aus_way,...
                'k_Aus',k_Aus,'Dirus',Dirus,'Dirus_way',Dirus_way,'ParallelCompute',pp_flag,'SolutionType',SolutionType,...
                'POIs',POIs,'scale',0,'showSphere',showSphere,'PLOT',PLOT_flag,'adj_drSphere',adj_drSphere,'HPC_flag',HPC_flag,...
                'Adapt_vAmp',Input_AdaptVamp,'Sim4life',S4l_input,'randomseed',s},Settings);
            [RMS,SNR,TSTOP,Out] = investBiologicalNoise(input_set{:});
            Input_VRsP.VRMATsP = permute(Out.VRMATsP,[3,2,1]);  Out.VRMATsP = [];
            Input_VRsP.VRMATsPstat = permute(Out.VRMATsPstat,[3,2,1]); Out.VRMATsPstat = [];            
            Input_VRsP.Aus = Out.Aus; Out.Aus = [];
            idx_zero = Input_VRsP.Aus==0;
            min_Aus = min(Input_VRsP.Aus(~idx_zero));
            Input_VRsP.Aus(idx_zero) = min_Aus; 
            if strcmpi(vibrField,'s4l')
                Input_VRsP.vampx = Out.vampx; Out.vampx = [];
                Input_VRsP.vampy = Out.vampx; Out.vampy = [];
                Input_VRsP.vampz = Out.vampx; Out.vampz = [];
            end
            RMSall(iAmp,idps).RMS = RMS;
            SNRall(iAmp,idps).SNR = SNR;
            Outall(iAmp,idps).Out = Out;
            TSTOPall(iAmp,idps).TSTOP = TSTOP;
            
            iAmps = 1:length(Amp); iAmps(iAmp==iAmps) = [];
            for iAmp = iAmps
                rng(s)
                fprintf('\n%s\n Amp: %i/%i, idps: %i/%i\n',datestr(now),iAmp,length(Amp),idps,length(Totaldps))
                Amp_input = Amp(iAmp);
                switch lower(vibrField)
                    case 's4l'
                        Input_AdaptVamp.Aus = Amp_input;
                        k_Aus = [];
                    case 'art'
                        switch lower(kAus_con)
                            case 'none'
                                k_Aus = kAus(1);
                            case 'amp'
                                k_Aus = kAus(iAmp);
                            case 'totaldps'
                                k_Aus = kAus(idps);
                            otherwise
                                error('false kAus_con input')
                        end
                    case 'none'
                        k_Aus = [];
                end
                input_set = horzcat({'rel_singlePeriod_flag',1,'Input_VRsP',Input_VRsP,'posDp',posDp,'OrienDipole',OrienDipole,...
                    'Totaldps',Totaldps(idps),'dps_run',dps_run(idps),'Aus',Amp_input,'Aus_way',Aus_way,...
                    'k_Aus',k_Aus,'Dirus',Dirus,'Dirus_way',Dirus_way,'ParallelCompute',pp_flag,'SolutionType',SolutionType,...
                    'POIs',POIs,'scale',0,'showSphere',showSphere,'PLOT',PLOT_flag,'adj_drSphere',adj_drSphere,'HPC_flag',HPC_flag,...
                    'Adapt_vAmp',Input_AdaptVamp,'Sim4life',S4l_input,'randomseed',s},Settings);
                [RMS,SNR,TSTOP,Out] = investBiologicalNoise(input_set{:});
                RMSall(iAmp,idps).RMS = RMS;
                SNRall(iAmp,idps).SNR = SNR;
                Outall(iAmp,idps).Out = Out;
                TSTOPall(iAmp,idps).TSTOP = TSTOP;
            end
            
            timeGrid_dprun = toc(startTimer_dprun);
            
            if calcStrengthCurve              
                for iov = 1:length(Amp)                    
                    for iPOIs = 1:length(RMSall(iov,idps).RMS)
                        % Decompose RMS structure
                        RMS_mat(iov,idps,iPOIs) = RMSall(iov,idps).RMS(iPOIs).RMS;
                    end                    
                end                
                if ~any(Amp-Amp(1)) && strcmpi(vibrField,'art') 
                    if length(kAus)>1
                        cSC_mode = 'kAus';
                    else
                        warning('calc strength mode is kAus but only one value?')
                    end
                end
                if isfield(Input,'cSC_mode')
                    cSC_mode = Input.cSC_mode;
                end
                switch lower(cSC_mode)
                    case 'amp'
                        
                        if strcmpi(vibrField,'art')
                            
                            switch lower(kAus_con)
                                case {'none','amp'}
                                    idx_kAus = 1;
                                    if length(kAus)>1 && any(kAus-kAus(1))
                                        warning('Amp optimisation with mulitple kAus?')
                                    end
                                    
                                case 'totaldps'
                                    idx_kAus = idps;
                                otherwise
                                    error('false kAus_con input')
                            end
                            
                        end
                    case'kaus'
                        idx_kAus = 1:length(kAus);
                    otherwise
                        error('incorrect cSC_mode')
                end
                
                
                input_set_init = horzcat({'rel_singlePeriod_flag',1,'Input_VRsP',Input_VRsP,'posDp',posDp,'OrienDipole',OrienDipole,...
                    'Aus_way',Aus_way,'Dirus',Dirus,'Dirus_way',Dirus_way,'ParallelCompute',pp_flag,'SolutionType',SolutionType,...
                    'POIs',POIs,'scale',0,'PLOT',PLOT_flag,'adj_drSphere',adj_drSphere,'HPC_flag',HPC_flag,...
                    'Adapt_vAmp',Input_AdaptVamp,'Sim4life',S4l_input,'randomseed',s},Settings);
                [requiredAmp_idps,requiredAmplower_idps,itermat_idps,cSC_mode] = cSC_cEG(Input,...
                    Amp,kAus(idx_kAus),kAus_con,vibrField,flaglim,evalMode,maxTime_dprun,maxTime_dprun,timeGrid_dprun,...
                    Totaldps(idps),dps_run(idps),RMS_mat(:,idps,:),Threshold,maxIter,input_set_init,SettingsStr,Input_str,s);
                 requiredAmp(idps) = requiredAmp_idps;
                 requiredAmplower(idps) = requiredAmplower_idps;
                 itermat(idps) = itermat_idps;
            end
        end
    case 'validate'
    otherwise
        error('wrong sim_method input' )
end
TimeSimulations = toc(startTimer);

%% decompose to reduce size of final file
To = whos('Outall','RMSall','SNRall','TSTOPall');
dataTo = sum([To(:).bytes]);
% if (dataTo/1024^3)>2
%     Decomp_struct_flag = 1;
% else
%     Decomp_struct_flag = 0;
% end
Decomp_struct_flag = 1;
if isfield(Input,'save_flag'); save_flag = Input.save_flag; end

Itgen_mat = zeros(length(Amp),length(Totaldps),max(ceil(Totaldps./dps_run)));
VRgen_mat = zeros(length(Amp),length(Totaldps),max(ceil(Totaldps./dps_run)));
RMS_info = {};
RMS_DB_info = {};
liPOIs = size(POIs,1);
if Decomp_struct_flag
    for iAmp = 1:length(Amp)
        for idps = 1:length(Totaldps)
            for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                % Decompose RMS structure
                RMS_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMS;
                Q2_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2;
                RMSptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSptherm;
                Q2ptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2ptherm;                
                RMSDOI_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSDOI;                
                Q2DOI_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2DOI;
                RMSDOIvr_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvr;                
                Q2DOIvr_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvr;                
                RMSDOIvrptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvrpthn;                
                Q2DOIvrptherm_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvrpthn;
                RMSDOIvrpOSCvr_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvrpOSCvr;                
                Q2DOIvrpOSCvr_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvrpOSCvr;
                RMSDOIvrpstatnoise_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).RMSDOIvrpstatnoise;                
                Q2DOIvrpstatnoise_mat(iAmp,idps,iPOIs) = RMSall(iAmp,idps).RMS(iPOIs).Q2DOIvrpstatnoise;
                
                RMS_info(iAmp,idps,iPOIs) = {RMSall(iAmp,idps).RMS(iPOIs).info};
                
                try
                %Decompose Out reconsSignal
                rSVR_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVR{iPOIs};
                rSVRDOI_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVRDOI{iPOIs};
                rSVRDOIvr_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVRDOIvibr{iPOIs};
                rSVRDOIptherm_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVRptherm{iPOIs};
                rSVRdvpov_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVRdvpov{iPOIs};
                rSVRdvpsn_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.reconsSignalVRdvpsn{iPOIs};
                
                trS_mat(iAmp,idps,iPOIs,:) = Outall(iAmp,idps).Out.timereconS{iPOIs};
                catch
                    fprintf('\ncould not decompose reconsSignal poi: %i',iPOIs)
                end
            end
            % Decompose SNR structure
            SNR_DOIvr_noiseall_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOIvr_NoiseAll;
            SNR_DOIvr_noiseallptherm_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOIvr_NoiseAllptherm;
            SNR_DOIvr_noiseoscvr_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOIvr_NoiseOSCvr;
            SNR_DOIvr_noisestat_mat(iAmp,idps,:,:) = SNRall(iAmp,idps).SNR.DOIvr_NoiseStatic;
            SNR_info(iAmp,idps,:) = SNRall(iAmp,idps).SNR.info(:);
            
            
                
            
            
            % Decompose Out
            VR_mat(iAmp,idps,:,:) = Outall(iAmp,idps).Out.VR;
            VRptherm_mat(iAmp,idps,:,:) = Outall(iAmp,idps).Out.VRptherm;
            VRDOI_mat(iAmp,idps,:,:) = Outall(iAmp,idps).Out.VRDOI;
            VRDOIstat_mat(iAmp,idps,:,:) = Outall(iAmp,idps).Out.VRDOIstat;
            VROSC_mat(iAmp,idps,:,:) = Outall(iAmp,idps).Out.VROSC;
            VROSCstat_mat(iAmp,idps,:,:) = Outall(iAmp,idps).Out.VROSCstat;
            VRstatnoise_mat(iAmp,idps,:,:) = Outall(iAmp,idps).Out.VRstatnoise;
            
            
            
            
            
            
            %Decompose TStOP
            
            
            StartSim_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP.startSim;
            if ~strcmpi(sim_method,'single_sim_totaldps')
            Itgen_mat(iAmp,idps,1:length(TSTOPall(iAmp,idps).TSTOP.Itgen)) = TSTOPall(iAmp,idps).TSTOP.Itgen;
            end
            VRgen_mat(iAmp,idps,1:length(TSTOPall(iAmp,idps).TSTOP.VRgen)) = TSTOPall(iAmp,idps).TSTOP.VRgen;
            endsim_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP.endsim;
            signalreconstruction_mat(iAmp,idps,:) = TSTOPall(iAmp,idps).TSTOP.signalreconstruction;
            
            
            Outall(iAmp,idps).Out.VR = [];
            Outall(iAmp,idps).Out.VRptherm = [];
            Outall(iAmp,idps).Out.VRDOI = [];
            Outall(iAmp,idps).Out.VRDOIstat = [];
            Outall(iAmp,idps).Out.VROSC = [];
            Outall(iAmp,idps).Out.VROSCstat = [];
            Outall(iAmp,idps).Out.VRstatnoise = [];   
            
            Outall(iAmp,idps).Out.reconsSignalVR = [];            
            Outall(iAmp,idps).Out.reconsSignalVRDOI = [];
            Outall(iAmp,idps).Out.reconsSignalVRDOIvibr = [];
            Outall(iAmp,idps).Out.reconsSignalVRptherm = [];
            Outall(iAmp,idps).Out.reconsSignalVRdvpov = [];
            Outall(iAmp,idps).Out.reconsSignalVRdvpsn = [];
            Outall(iAmp,idps).Out.timereconS = [];
            
            
            

        end
    end
    disp('succes');
    Ta = whos('*_mat');
    dataTa = sum([Ta(:).bytes]);
    namestosave = horzcat({Ta(:).name},{'SNR_info','TSTOPall','POIs','RMS_info','Input','Outall','Decomp_struct_flag','TimeSimulations'});
    if calcStrengthCurve
        namestosave = horzcat(namestosave,{'requiredAmp','requiredAmplower','itermat','cSC_mode'});
    end
    disp('');
    disp(['ended, start saving']);
    if save_flag
        idxdoublepoint = find(SettingsStr==':');
        if isempty(idxdoublepoint)
            save_str_int = SettingsStr;
        else
            save_str_int = SettingsStr([1:idxdoublepoint-1,idxdoublepoint+2:end]);
        end
%         savename = sprintf('Results_cEG%s_%s_%s_%s.mat',savename_add,save_str_int,...
%             Input_str(regexp(Input_str,'_sM')+3:end),datestr(now,'mm-dd-yy_HHMM'));
savename = sprintf('Results_cEG_%s_%s_%s.mat',save_str_int,...
            Input_str(regexp(Input_str,'_cEG')+5:end-2),datestr(now,'mm-dd-yy_HHMM'));
        save(savename,namestosave{:});
    end
else
disp(' ');
disp(['ended, start saving']);

if save_flag
    idxdoublepoint = find(SettingsStr==':');
    if isempty(idxdoublepoint)
        save_str_int = SettingsStr;
    else
        save_str_int = SettingsStr([1:idxdoublepoint-1,idxdoublepoint+2:end]);
    end
    namestosave = {'RMSall','SNRall','Outall','TSTOPall','Input','POIs','Decomp_struct_flag','TimeSimulations'};
    if calcStrengthCurve
        namestosave = horzcat(namestosave,{'requiredAmp','requiredAmplower','itermat','cSC_mode'});
    end
    t = whos('Outall');
%     savename = sprintf('Results_cEG%s_%s_%s_%s.mat',savename_add,save_str_int,...
%             Input_str(regexp(Input_str,'_sM')+3:end),datestr(now,'mm-dd-yy_HHMM'));
savename = sprintf('Results_cEG_%s_%s_%s.mat',save_str_int,...
            Input_str(regexp(Input_str,'_cEG')+5:end-2),datestr(now,'mm-dd-yy_HHMM'));
    if (t.bytes/1024^3)>2        
        save(savename,namestosave{:},'-v7.3');
    else
        save(savename,namestosave{:});
    end
fprintf('save_succes');
end
end
if HPC_flag && isdeployed
    pause(60)
    exit
end
end
