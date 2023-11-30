function [PLOT,PLOT_RMS,Display,plotIarray_flag,plotUSwaves_flag,HPC_flag,showSphere,...
    ParallelCompute_flag,ThermalNoiseAmp,scale_flag,AdaptVamp_flag,input_posDp_flag,...
    sel_oscil_from_gendPos_flag,coll_singlePeriod_flag,rel_singlePeriod_flag,adj_drSphere_flag,...
    all_osc_flag,S4l_flag,genDatabase_flag,noRecon_flag,VR_given_flag,eval_method,vib_interp_type,interp3_method,Dirus_way,...
    Aus_way,k_Aus,nLayers,sigma,d,dI,meanI,stdI,Ifun,fus,Aus,Phaseus,Dirus,fbandwidth,resUS,Tend,...
    SolutionType,POIs,OrienDipole,Totaldps,dps_run,dpDistribution,dpOrientation,...
    dpI_time,dpI_space,vmax_AVa,Aus_AVa,posDp,ROI_OSC,Input_VR,Input_VRsP,Input_AdaptVamp,S4l_Input,mrunpBatch,randomseed,pulsed,Param] = declare_parameters(varargin)
%% Initialise default variables
tf_str = {'Off','On'};

PLOT = 1;                        %Display plots generated in reconstruction fase (calFourierGetsignal) +
                                 %Generate plot of oscillation seen from first 3 POIs
                                 %Generate plot of the RMS difference
                                 %Generate large 4*3 plot conatining frequency
                                 %Domains near 0 and fus
PLOT_RMS = 0;                               
                        
Display = 0;                     %Display progress
plotIarray_flag = 0;             %Create a figure with selection of applied currents
plotUSwaves_flag = 0;            %Plot vibration of a select group of dipoles (can be  xyz or ellipsoidal
HPC_flag = 0;                    %HPC flag mtit cannot be used
showSphere='';                   % '','nopoi','wpoi' show before each iteration position of dipoles with or without POIs
ThermalNoiseAmp = 10^-2;         %Amplitude of thermalNois ein uV
scale_flag = 1;                  %Rescale POIs to be on outer sphere
AdaptVamp_flag = 0;              %Set max vibration,flag only necessary when S4L input is used if 0 q,plitude fro, S4L si, is used
input_posDp_flag = 0;            %Dipole potions given as input if not DOI is declared in cortex 0,0,Rsphere-dRSphere
sel_oscil_from_gendPos_flag = 0; %Assign OSC_flag to dipoles that are in certain region definded by ROI_OSC
ROI_OSC = [];                    %Region of Interest (region that contains oscillating dipoles)
adj_drSphere_flag = 0;           %Adjust drSphere = max radius where dipoles are located, if 1: all dipoles are located deeper than dipole of interest
all_osc_flag = 0;                %If  this flag = 1 all dipoles generated below will oscilate as well. By default only dipoles (posDp) that are given as input will be oscillating, except when ROI is used.  
S4l_flag = 0;                    %S4l_flag, a velocity field created in S4L determines the oscillating directions and amplitudes
S4l_Input = [];                  %Structure containing all S4l info when S4l flag =1;
genDatabase_flag = 0;            %Stop after VR is simulated => no reconstruction, eg in case of database simulation
noRecon_flag = 0;                %Same purpose as genDatabase_flag (could be merged?)
VR_given_flag = 0;               %VR already determined only reconstruction necessary
Input_VR = [];                   %structure containing previously determined potential field (VR = total,VRDOI = DOI only, VRosc= all oscilaters;VRnoise)
coll_singlePeriod_flag = 0;      % database like method where one period of each dipole is saved
rel_singlePeriod_flag = 0;       % database like method where one period of each dipole is used (released vs collect)
Input_VRsP = [];                 % input containing single Period data
eval_method = 'normal';          %'normal','database' using database ,'database_f' using faster version of database,'validate' both database and normal to compare
vib_interp_type = 'linear';      %How is vibration amplitude interpolated when using a database
interp3_method = 'linear';       %Phase and vibration determination via interpolation from Velocity field(S4L)
Dirus_way = 'input';             %'dipole_dir','radial','input', Direction of oscilation
Aus_way = 'normal';              %'kexp' exponential decreasing amplitude from DOI starting on.
                                 %'kscale' all other dipoles are a factor differ a factor k_AUS
                                 %'normal' Aus is not altered = equal to input
                                 %Aus DOI = AUS(1) 
                                 
k_Aus = 1;                       %[m^-1] case kexp: for x=k [m^-1] Aus(i) = 0.36*Aus(DOI)
                                 %[-] cae kscale: Aus(i) = k_Aus*Aus(DOI)
                                 
mrunpBatch = 5;                  % default multiplication number if batch run stelected                                                    


% Electrical properties
Npercolumn = 1e4;                %number on neurons per cortical column
meanI = 1*1e-3*Npercolumn;       %µA strength of single neuron is 1nA? see raportations
sigma = 0.33;                    %S/m conductance brain tissue
d = 500*10^-6; %distance between current source and current sink
dI = meanI*d;                    %dipole moment
stdI = 3;                        % standard deviation of mean current
%Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+1/10*double(t(1)~=t0).*gaussmf(t,[0.0005,t0+0.0005]);
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0; %alpha function
Tau = 0.005;                     % [s]
AlphaDelay = 0;                 %delay in alpha function [s]
Ifun = @(t) Alphafun(t,Tau,AlphaDelay); %default function applied as current on DOI
%Ifun = @(t) sin(2*pi*80.*t);

% US options
fus = 1e6;                     %[Hz] frequency of US wave
Aus = 10e-9;                    %[m] Amplitude of displacement
Phaseus = 0;                    %[rad] phase on vibration
Dirus = [0,0,1];                % vibration direction
fbandwidth = 1e3;               %bandwidth used for reconstruction
Input_AdaptVamp = [];           %Adapt max vibration amplitude when given by S4L
vmax_AVa = [];                  %max velocity in S4L sim
Aus_AVa = [];                   %max displacement amp wanted

% Simulation settings
resUS = 20;                     %resolution of US wave: # points per period
Tend = 0.025;                   %end time [ms]
pulsed = struct();
pulsed.flag = 0;
pulsed.prp = 0;
pulsed.dc = 1;


% Geometric  properties
SolutionType = '4SphereS8.7R25';%Select which headmodel is used see getSettings

% POIs by default (on unit sphere, will be scaled with scale flag) second
% position is of centre to 0,0,1
POIs = [0,0,1;cos(0)*sin(pi/180),sin(0)*sin(pi/180),cos(pi/180);0,1,0;0,1/sqrt(2),1/sqrt(2)];

% Dipoles
OrienDipole = [0,0,1];              %orientation of input dipoles posDp
Totaldps = 10;                     % Total number of dipoles used                  
dps_run = 10;                      % Total number of dipoles in each run (best not bigger than 1e3: matrix to heavy)
dpDistribution = 'fibonaccisingle'; % 'random1','random2','random3','concentricsingle'
                                    %'fibonaccisingle',fibonaccimultiple','S4l_hotspots'
                                    %distribution ot Totaldps-#posDp
                                    %dipoles see GendpPos more info
nLayers = 10;                       %# Of layers in fibonaccimultiple
dpOrientation = 'radial';           %'radial','random' orientation
dpI_time = 'artificialPSD';         %'randomdelay','randomi','randominewgen','randomirandomdelay','sync','artificialPSD','artificialPSD_WC','artificialPSD_BC'
                                    %time course of currents applied to
                                    %dipoles see getIarray
dpI_space = 'normal';               %'same','normal','lognormal'%
                                    %distribution of Max amplitudes of current in dipoles
posDp = [];                         %dipoles position given as input
                                    

ParallelCompute_flag = Totaldps>100;        %Use parallel computation: SimDipoleOSCPC is used (in SimDipoleOSC)
randomseed = rng('shuffle','twister');
%% change default variables
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'Plot'))
            PLOT = varargin{find(strcmpi(varargin,'Plot'))+1};
            % Plot settings if 1  change defaults flags below
            if PLOT
                plotIarray_flag = 1;
                showSphere = 'wpoi';
                Display = 1;
                plotUSwaves_flag = 1;
            end
        end
        if any(strcmpi(varargin,'Plot_RMS'))
            PLOT_RMS = varargin{find(strcmpi(varargin,'Plot_RMS'))+1};
        end
        

        if any(strcmpi(varargin,'Aus'))
            Aus = varargin{find(strcmpi(varargin,'Aus'))+1};
        end
        if any(strcmpi(varargin,'Aus_way'))
            Aus_way = varargin{find(strcmpi(varargin,'Aus_way'))+1};
        end
        if any(strcmpi(varargin,'all_osc'))
            all_osc_flag = varargin{find(strcmpi(varargin,'all_osc'))+1};
        end
        if any(strcmpi(varargin,'Adapt_vAmp'))
            Input_AdaptVamp = varargin{find(strcmpi(varargin,'Adapt_vAmp'))+1};
            AdaptVamp_flag = Input_AdaptVamp.flag;
            if AdaptVamp_flag
            vmax_AVa = Input_AdaptVamp.vmax;
            Aus_AVa = Input_AdaptVamp.Aus;
            end
        end
        if any(strcmpi(varargin,'adj_drSphere'))
            adj_drSphere_flag = varargin{find(strcmpi(varargin,'adj_drSphere'))+1};
        end
        if any(strcmpi(varargin,'coll_singlePeriod_flag'))
            coll_singlePeriod_flag = varargin{find(strcmpi(varargin,'coll_singlePeriod_flag'))+1};
        end
        if any(strcmpi(varargin,'Dirus'))
            Dirus = varargin{find(strcmpi(varargin,'Dirus'))+1};
        end
        if any(strcmpi(varargin,'Dirus_way'))
            Dirus_way = varargin{find(strcmpi(varargin,'Dirus_way'))+1};
        end
        if any(strcmpi(varargin,'dpDistribution'))
            dpDistribution = varargin{find(strcmpi(varargin,'dpDistribution'))+1};
        end
        if any(strcmpi(varargin,'dpI_space'))
            dpI_space = varargin{find(strcmpi(varargin,'dpI_space'))+1};
        end
        if any(strcmpi(varargin,'dpI_time'))
            dpI_time = varargin{find(strcmpi(varargin,'dpI_time'))+1};
        end
        if any(strcmpi(varargin,'dpOrientation'))
            dpOrientation = varargin{find(strcmpi(varargin,'dpOrientation'))+1};
        end
        if any(strcmpi(varargin,'dps_run'))
            dps_run = varargin{find(strcmpi(varargin,'dps_run'))+1};
            dps_run = min(dps_run,100); %if higher than 100 difficulties with Iarray matrix (too much Memory)
        end
        if any(strcmpi(varargin,'eval_method'))
            eval_method = varargin{find(strcmpi(varargin,'eval_method'))+1};
        end
        if any(strcmpi(varargin,'fus'))
            fus = varargin{find(strcmpi(varargin,'fus'))+1};
        end
        if any(strcmpi(varargin,'genDatabase'))
            genDatabase_flag = varargin{find(strcmpi(varargin,'genDatabase'))+1};
        end
        if any(strcmpi(varargin,'HPC_flag'))
            HPC_flag = varargin{find(strcmpi(varargin,'HPC_flag'))+1};
        end
        if any(strcmpi(varargin,'Ifun'))
            Ifun = varargin{find(strcmpi(varargin,'Ifun'))+1};
        end
        if any(strcmpi(varargin,'k_Aus'))
            k_Aus = varargin{find(strcmpi(varargin,'k_Aus'))+1};
        end
        if any(strcmpi(varargin,'meanI'))
            meanI = varargin{find(strcmpi(varargin,'meanI'))+1};
        end
        if any(strcmpi(varargin,'mrunpBatch'))
            mrunpBatch = varargin{find(strcmpi(varargin,'mrunpBatch'))+1};
        end
        if any(strcmpi(varargin,'nLayers'))
            nLayers = varargin{find(strcmpi(varargin,'nLayers'))+1};
        end
        if any(strcmpi(varargin,'noRecon'))
            noRecon_flag = varargin{find(strcmpi(varargin,'noRecon'))+1};
        end
        if any(strcmpi(varargin,'OrienDipole'))
            OrienDipole = varargin{find(strcmpi(varargin,'OrienDipole'))+1};
        end
        if any(strcmpi(varargin,'ParallelCompute'))
            ParallelCompute_flag = varargin{find(strcmpi(varargin,'ParallelCompute'))+1};
            inputPC_flag = 1;
        else
            inputPC_flag = 0;
        end
        if any(strcmpi(varargin,'Phaseus'))
            Phaseus = varargin{find(strcmpi(varargin,'Phaseus'))+1};
        end     
        if any(strcmpi(varargin,'POIs'))
            POIs = varargin{find(strcmpi(varargin,'POIs'))+1};
        end
        if any(strcmpi(varargin,'posDp'))
            posDp = varargin{find(strcmpi(varargin,'posDp'))+1};
            if ~isempty(posDp)
                input_posDp_flag = 1;
            end
        end
        if any(strcmpi(varargin,'pulsed'))
            pulsed = varargin{find(strcmpi(varargin,'pulsed'))+1};
        end
        if any(strcmpi(varargin,'randomseed'))
            randomseed = varargin{find(strcmpi(varargin,'randomseed'))+1};
        end
        if any(strcmpi(varargin,'rel_singlePeriod_flag'))
            rel_singlePeriod_flag = varargin{find(strcmpi(varargin,'rel_singlePeriod_flag'))+1};
            if rel_singlePeriod_flag
                Input_VRsP = varargin{find(strcmpi(varargin,'Input_VRsP'))+1};
            end
            if isempty(Input_VRsP)
                error('release singlePeriod_flag but no input VRsP given')
            end
        end
        if any(strcmpi(varargin,'resUS'))
            resUS = varargin{find(strcmpi(varargin,'resUS'))+1};
        end
        if any(strcmpi(varargin,'ROI_OSC'))
            ROI_OSC = varargin{find(strcmpi(varargin,'ROI_OSC'))+1};
            sel_oscil_from_gendPos_flag = 1;
        end
        if any(strcmpi(varargin,'scale'))
            scale_flag = varargin{find(strcmpi(varargin,'scale'))+1};
        end
        if any(strcmpi(varargin,'showSphere'))
            showSphere = varargin{find(strcmpi(varargin,'showSphere'))+1};
        end
        if any(strcmpi(varargin,'Sim4life'))
            S4l_Input = varargin{find(strcmpi(varargin,'Sim4life'))+1};
            S4l_flag = S4l_Input.S4l_flag;
        end
        if any(strcmpi(varargin,'SolutionType'))
            SolutionType = varargin{find(strcmpi(varargin,'SolutionType'))+1};
        end
        if any(strcmpi(varargin,'stdI'))
            stdI = varargin{find(strcmpi(varargin,'stdI'))+1};
        end
        if any(strcmpi(varargin,'Tend'))
            Tend = varargin{find(strcmpi(varargin,'Tend'))+1};
        end
        if any(strcmpi(varargin,'ThermalNoiseAmp'))
            ThermalNoiseAmp = varargin{find(strcmpi(varargin,'ThermalNoiseAmp'))+1};
        end
        if any(strcmpi(varargin,'Totaldps'))
            Totaldps = varargin{find(strcmpi(varargin,'Totaldps'))+1};
            if ~inputPC_flag
                ParallelCompute_flag = Totaldps>100;
            end            
        end
        if any(strcmpi(varargin,'vib_interp_type'))
            vib_interp_type = varargin{find(strcmpi(varargin,'vib_interp_type'))+1};
        end
        if any(strcmpi(varargin,'VR_Input'))
            Input_VR = varargin{find(strcmpi(varargin,'VR_Input'))+1};
            VR_given_flag = Input_VR.flag;
        end
    end
end

if Display
fprintf('\nModel parameters\n\nFlags:\n------------\n')
fprintf('\tPlot: %s\n',tf_str{PLOT+1});
fprintf('\tPlot Iarray: %s\n',tf_str{plotIarray_flag+1});
fprintf('\tPlot US waves_flag: %s\n',tf_str{plotUSwaves_flag+1});
fprintf('\tHPC flag: %s\n',tf_str{HPC_flag+1});
fprintf('\tParralelcomputing flag: %s\n',tf_str{ParallelCompute_flag+1});
fprintf('\tShow sphere: %s\n',showSphere);
fprintf('\tscale flag: %s\n',tf_str{scale_flag+1});
fprintf('\tAdapt vamp s4l input: %s\n',tf_str{AdaptVamp_flag+1});
fprintf('\tposDp as input: %s\n',tf_str{input_posDp_flag+1});
fprintf('\tselect oscillators form ROI: %s\n',tf_str{sel_oscil_from_gendPos_flag+1});
fprintf('\tadj max radius of dipole locations: %s\n',tf_str{adj_drSphere_flag+1});
fprintf('\tall dipoles oscillate: %s\n',tf_str{all_osc_flag+1});
fprintf('\tS4l  input: %s\n',tf_str{S4l_Input.S4l_flag+1});
fprintf('\tgenerate Database : %s\n',tf_str{genDatabase_flag+1});
fprintf('\tno reconstruciton of DOI signal: %s\n',tf_str{noRecon_flag+1});
fprintf('\tVR already calculated: %s\n',tf_str{VR_given_flag+1});
fprintf('\t collect single Period of each dipole at each POI: %s\n',tf_str{coll_singlePeriod_flag+1})
fprintf('\t Release single Period of each dipole at each POI: %s\n',tf_str{rel_singlePeriod_flag+1})
fprintf('\t multiplicator batch run: %i\n',mrunpBatch)
fprintf('\nParameter deriviation settings\n---------\n')
fprintf('\tShow sphere: %s\n',showSphere);
fprintf('\tevaluation method: %s\n',eval_method);
fprintf('\tinterpolation of vibration when database used: %s\n',vib_interp_type);
fprintf('\tinterpolation of vibration when S4l fields: %s\n',interp3_method);
fprintf('\tinterpolation of vibration when database used: %s\n',vib_interp_type);
fprintf('\tmethod how Dirus is derived: %s\n',Dirus_way);
fprintf('\tmethod how Aus is derived: %s\n',Aus_way);
fprintf('\tk_Aus: %5.2f\n',k_Aus);
fprintf('\tmaximal current strength DOI: %5.2f\n',meanI);
fprintf('\tconductivity brain matter: %5.2f\n',sigma);
fprintf('\tstandard deviation mean current amplitude: %5.2f\n',stdI);
fprintf('\ttimconstant alpha function: %5.2f\n',Tau);
fprintf('\tDelay alpha function: %5.2f\n',AlphaDelay);
fprintf('\tultrasound frequency: %5.2f\n',fus);
fprintf('\tUltrasound amplitude: %5.2f\n',Aus);
fprintf('\tPhaseUS: %5.2f\n',Phaseus);
strforplot = sprintf('%5.2e, ',Dirus);
fprintf('\tDirectionUS: [%s]\n',strforplot(1:end-2));
fprintf('\tfrequency range of interest: %5.2f\n',fbandwidth);
fprintf('\tresolution US wave: %5.2f\n',resUS);
fprintf('\tend of simulation: %5.4f\n',Tend);
fprintf('\tSimulation method: %s\n',SolutionType);
fprintf('\tend of simulation: %5.4f\n',Tend);
if size(POIs,1)>4
    fprintf('\t only first 4 POIs printed\n')
end    
for ipoi = 1:min(size(POIs,1),4)
    strforplot = sprintf('%5.2e, ',POIs(ipoi,:));
    fprintf('\tPOI %i: [%s]\n', ipoi,strforplot(1:end-2))
end
if size(OrienDipole,1)>4
    fprintf('\t only first 4 Dipole orientations printed\n')
end    
for ipoi = 1:min(size(OrienDipole,1),4)
    strforplot = sprintf('%5.2e, ',OrienDipole(ipoi,:));
    fprintf('\tdipole orientation %i: [%s]\n', ipoi,strforplot(1:end-2))
end
if size(posDp,1)>4
    fprintf('\t only first 4 Dipole orientations printed\n')
end    
for ipoi = 1:min(size(posDp,1),4)
    strforplot = sprintf('%5.2e, ',posDp(ipoi,:));
    fprintf('\tdipole positions %i: [%s]\n', ipoi,strforplot(1:end-2))
end

 fprintf('\ttotal number of dipoles: %5.2f\n',Totaldps); 
 fprintf('\tdipole per run: %5.2f\n',dps_run);
 fprintf('\tdpdistribution: %s\n',dpDistribution);
 if strcmpi(dpDistribution,'fibonaccimultiple')
      fprintf('\tnubmer of layers : %i\n',nLayers);
 end
 fprintf('\tdipole orientation: %s\n',dpOrientation);
 fprintf('\tdipole current course: %s\n',dpI_time);
 fprintf('\tamplitude distribution over space: %s\n',dpI_space);

end
disp('');
%store important parameters in Param
Param.ThermalNoiseAmp = ThermalNoiseAmp; Param.scale_flag = scale_flag; Param.AdaptVamp_flag = AdaptVamp_flag; Param.input_posDp_flag = input_posDp_flag;
Param.adj_drSphere_flag = adj_drSphere_flag; Param.all_osc_flag = all_osc_flag; Param.S4l_flag = S4l_flag; Param.genDatabase_flag = genDatabase_flag;
Param.coll_singlePeriod_flag = coll_singlePeriod_flag; Param.rel_singlePeriod_flag = rel_singlePeriod_flag;
Param.eval_method = eval_method; Param.vib_interp_type = vib_interp_type; Param.interp3_method = interp3_method; Param.Dirus_way = Dirus_way;
Param.Aus_way = Aus_way; Param.k_Aus = k_Aus; Param.nLayers = nLayers; Param.sigma = sigma; Param.d = d; Param.dI = dI; Param.meanI = meanI;
Param.stdI = stdI; Param.Ifun = Ifun; Param.fus = fus; Param.Aus = Aus; Param.Phaseus = Phaseus; Param.Dirus = Dirus; Param.fbandwidth = fbandwidth;
Param.resUS = resUS; Param.Tend = Tend; Param.SolutionType = SolutionType; Param.POIs = POIs; Param.OrienDipole = OrienDipole;
Param.Totaldps = Totaldps; Param.dps_run = dps_run; Param.dpDistribution = dpDistribution; Param.dpOrientation = dpOrientation;
Param.dpI_time = dpI_time; Param.dpI_space = dpI_space; Param.vmax_AVa = vmax_AVa; Param.Aus_AVa = Aus_AVa; Param.posDp = posDp;
Param.ROI_OSC = ROI_OSC; Param.Input_VR = Input_VR; Param.Input_AdaptVamp = Input_AdaptVamp; Param.randomseed = randomseed; 
Param.pulsed = pulsed;

if coll_singlePeriod_flag && ~strcmpi(eval_method,'normal')
    warning('collect single period only possible with normal run')
end

end