function createfig(Input,Outall,RMSall,SNRall,POIs,Debug_flag,PSM_flag,s,nos4l_flag,Param,date_sel)
% create figures of results in calcErrorGrid. Outputs all figures. 
%Inputs: Input: the input structure used in calcErrorGrid. contains some
%               chosen parameter values.
%        Outall: a structure with results of the calcErrorGrid
%        RMSall: structure with root mean square error results
%        SNRall: structure with Signal to noise ratios
%        POIs: Point of interest positions
%        Debug_flag: debug flag => reduces plots and computation time to
%        debug code
%        PSM_flag: create figure of potential at sphere boundary with
%        multiple dipoles
%        s: randomseed for consistency between results
%        nos4l_flag: results not based on s4l result
%        Param: prameters used in investbiological noise
%        Date_sel: date of selected results => code differs for later
%        results

figpos = [89, 83, 1567, 1015];
DOIindice = 1;               % dipole of interest is by default at 1          
Amp = Input.Amp;             % Amplitudes run
Totaldps = Input.Totaldps;   % number of dipoles included
lAmp  = length(Amp);         % length of Amp vector
lTdps = length(Totaldps);    % length of dipoles vector
Tend = 0.025;                % end simulation time (default) changed if parameter struct included
resUS = 20;                  % resolution US wave # points per period
Alphafun = @(t,Tau,t0) double((t-t0)>=0).*((t-t0)./Tau).*exp(1-(t-t0)./Tau)+0;
Tau = 0.005;         
AlphaDelay = 0;
Ifun = @(t) Alphafun(t,Tau,AlphaDelay);  % function of interest
meanI = 10;                  % mean strength of current through dipoles
Window = 1000;               % window in reverse fourier transform

iAmp_sp = 1; %index Amp used for surface sphere plot
idps_sp = 2; %index dps used for surface sphere plot

%% reobtain data
dps_run = Input.dps_run;
Settings = Input.Settings;  SettingsStr = Input.SettingsStr;
%Fignr = Input.Fignr;
%HPC_flag = Input.HPC_flag;
if ~nos4l_flag
S4l_flag = Input.S4l_flag; S4l_fn = Input.Filename_S4l; S4l_loc = Input.Location_S4l;
end
dpDistribution = Settings{find(strcmpi(Settings,'dpDistribution'))+1};
dpOrientation = Settings{find(strcmpi(Settings,'dpOrientation'))+1};
fus = Settings{find(strcmpi(Settings,'fus'))+1};
sigma = 0.33;       % conductivity of brain tissue
SphereRes = 30;     % resolution of sphere
d = 500*10^-6; %distance between current source and current sink
nLayers = 10;  % number of layers in multi fibonacci

if ~isempty(Param)
    Window = Param.fbandwidth;
    Tend = Param.Tend;
    resUS = Param.resUS;
    meanI = Param.meanI;
    Ifun = Param.Ifun;
    sigma = Param.sigma;
    d = Param.d;
    nLayers = Param.nLayers;
    
end

if ~nos4l_flag
    % when s4l region is decided by position highest vibration amplitude
    % files saved in folders centerfocus or cortex focus (s4l_loc)
if contains(S4l_loc,'CenterFocus')
    title_suffix = 'ROI: Deep';
elseif contains(S4l_loc,'CortexFocus')
    title_suffix = 'ROI: Cortex';
end
end
if isfield(Input,'ROI')
    title_suffix = ['ROI: ',Input.ROI];
else
    title_suffix = '';
end
Model = Input.Model;
switch lower(Model)
    case 'human_scalp'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.082;
        RBrain = 0.07;
        dRBrain = 0.005;
        title_prefix = 'Measure: Human Scalp,  ';
    case 'human_cortical'
        SolutionType = '4SphereS8.7R25';
        RPOI = 0.07;
        RBrain = 0.07;
        dRBrain = 0.005;
        title_prefix = 'Measure: Human Cortical,  ';
    case 'human_air'
        SolutionType = '4sphereS8.7R25';
        RPOI = 0.087;
        RBrain = 0.07;
        dRBrain = 0.005;
        title_prefix = 'Measure: Human air,  ';
    case 'mouse_fair0'
        SolutionType = 'Mouse4Sphere~fair0';
        RPOI = 0.0059;
        RBrain = 0.0046;
        dRBrain = 0.0005;
        title_prefix = 'Measure: Mouse Scalp,  ';
    otherwise
        error('false input model')
end

[Options,RSphere] = getSettings(SolutionType,1);
Options{find(strcmpi(Options,'tAir'))+1}=0;
if strcmpi(SolutionType,'Mouse4Sphere~fair0') || strcmpi(SolutionType,'Mouse4Sphere~fair1')
    dRSphere = 0.0005;
else
    dRSphere = 0.005;
end


CSource = Outall(iAmp_sp,idps_sp).Out.CSource;
CSink = Outall(iAmp_sp,idps_sp).Out.CSink;
OSCindices = Outall(iAmp_sp,idps_sp).Out.OSCindices;
POIs = Outall(iAmp_sp,idps_sp).Out.POIs;

if Debug_flag
    Totaldps_plot = min(100,size(CSource,1));
    idx = [1,sort(randperm(size(CSource,1)-1,Totaldps_plot))+1];
    CSource_plot = CSource(idx,:);
    CSink_plot = CSink(idx,:);
    %OSCindices_plot = OSCindices(1:Totaldps_plot,:);
end

Tend = ceil(fus*Tend)/fus; % recalculate Tend such that frequency spectrum contains fus
dt = (resUS*fus)^-1;
Tsim = 0:dt:Tend;
inputSignal = [Tsim;meanI.*Ifun(Tsim)];
mirror = zeros(size(POIs,1),1);
wus = 2*pi*fus;
SF = 1/(wus);

%see declare_vibrationfunc.m
if nos4l_flag

    Dirus_way = Input.Dirus_way;
    Aus_way = Input.Aus_way;
    Aus = Input.Amp(iAmp_sp);
    kAus = Input.kAus;
    k_Aus = Input.kAus(iAmp_sp);
    Phaseus = Input.Phaseus;
    Dirus = Input.Dirus;
    switch lower(Dirus_way)
        case 'dipole_dir'

            Dirus_DOI = Dirus(DOIindice,:);
            Dirus = CSource-CSink;
            Dirus = Dirus./vecnorm(Dirus);
            Dirus(DOIindice,:) = Dirus_DOI;
            
        case 'radial'
            
            Dirus_DOI = Dirus(DOIindice,:);
            Dirus = (CSource+CSink)/2;
            Dirus = Dirus./vecnorm(Dirus);
            Dirus(DOIindice,:) = Dirus_DOI;
        case 'input'
            % no change on input Dirus
            
        otherwise 
            error('false input Dirus_way')

    end
    switch lower(Aus_way)
        case 'kexp'
           
            if isrow(Aus)
                Aus = Aus';
            end
            Aus_DOI = Aus(DOIindice,:);
            kfun = @(x)  Aus_DOI*exp(-k_Aus*x);
            Rdiffdppos = vecnorm((CSource+CSink-CSource(DOIindice,:)-CSink(DOIindice,:))/2,2,2);
            Aus = kfun(Rdiffdppos);
            Aus(DOIindice,:) = Aus_DOI;            
            
        case 'kscale'

            if isrow(Aus)
                Aus = Aus';
            end
            Aus_DOI = Aus(DOIindice,:);
            Aus = Aus_DOI*k_Aus*size(CSource,1);
            Aus(DOIindice,:) = Aus_DOI;
            
        case 'normal'
            % no change on input Aus
        otherwise
            error('false input Aus_way')
    end
    USwave_varlengths = [size(Dirus,1), size(Aus,1), length(wus), size(Phaseus,1)];
    [max_USwave_varlengths] = max(USwave_varlengths);
    idx_mUv_l = max_USwave_varlengths~=USwave_varlengths;
    if  any(USwave_varlengths-USwave_varlengths(1))
        if any(USwave_varlengths(idx_mUv_l)-1)
            error('length of US variables should either be the same for all or one should be bigger while others are one')
        else
            repmat_vals = max_USwave_varlengths-USwave_varlengths+1;
        end
        Dirus = repmat(Dirus,repmat_vals(1),1);
        Aus = repmat(Aus,repmat_vals(2),1);
        wus = repmat(wus,repmat_vals(3),1);
        Phaseus = repmat(Phaseus,repmat_vals(4),1);
        
    end
    USwave_all = @(t,idx) Dirus(idx,:).*Aus(idx,:).*sin(wus(idx).*t+Phaseus(idx,:));
else
vampx = Outall(1,1).Out.vampx;
vampy = Outall(1,1).Out.vampy;
vampz = Outall(1,1).Out.vampz;
vphasex = Outall(1,1).Out.phasex;
vphasey = Outall(1,1).Out.vphasey;
vphasez = Outall(1,1).Out.vphasez;

USwave_all = @(t,idx) SF.*[vampx(idx).*sin(wus.*t+vphasex(idx)),vampy(idx).*sin(wus.*t+vphasey(idx)),vampz(idx).*sin(wus.*t+vphasez(idx))];
end
% determine theta of USwave (see Invest biological noise function line 319)
tvals = linspace(0,1/fus,1e3)'; tvals = tvals(1:end-1);
USwave_sel = USwave_all(tvals,1);
USwave_sel = USwave_sel-mean(USwave_sel);
amp_vals = max(USwave_sel);
USwave_sel = USwave_sel(1:floor(length(tvals)/2),:);
USwave_t0 = USwave_sel(1,:);
dUSwave_t1t0= USwave_sel(2,:)-USwave_sel(1,:);

[~,cross_zero] = min(abs(USwave_sel));
theta_vals(1) = wus(DOIindice)*(1/(2*fus)*double(USwave_t0(1)>0)-tvals(cross_zero(1)));
theta_vals(2) = wus(DOIindice)*(1/(2*fus)*double(USwave_t0(2)>0)-tvals(cross_zero(2)));
theta_vals(3) = wus(DOIindice)*(1/(2*fus)*double(USwave_t0(3)>0)-tvals(cross_zero(3)));
theta_vals(cross_zero==1) = 0; theta_vals(cross_zero==1 & dUSwave_t1t0<0) = pi;
% USwave_t0 = USwave_sel(t0);USwave_t1 = USwave_sel(t1);
% theta_vals = atan2(sin(wus.*t1).*ones(1,3),(USwave_t1./USwave_t0-cos(wus.*t1)));
% % nanmean is for special case when theta is zero ==> sin(theta) is 0 gives
% % devision by zero = nan
% amp_vals = nanmean(vertcat(abs(USwave_t0./sin(theta_vals)),abs(USwave_t1./sin(wus.*t1+theta_vals))));

% circular mean of different directions
avg_vect = nansum(amp_vals.*exp(complex(0,theta_vals)))./nansum(amp_vals);
theta_global = atan2(imag(avg_vect),real(avg_vect));
%% start plotting

%% Plot reconstructed signals



idpss = 9;                %which dipole set needs to be plotted
iAmps = [1,3,5];          % which amplitudes
liAmp = length(iAmps);
if liAmp==3
    eflag = 1;
elseif liAmp>3
    eflag = 1;
    warning ('max three levels')
else
    eflag = 0;
end

%difference between results saved before september 28 2020
% averaged signal is seperate field
NoiseAmp = [10^-2,10^-2*sqrt(2*Window),10^-2*sqrt(1/dt)]; %thermal noise variants
if datenum(date_sel)<=datenum(datetime(2020,9,28,0,0,0)) && datenum(date_sel)>=datenum(datetime(2020,2,28,0,0,0))
    fnamesall = {{'reconsSignalVR','reconsSignalVR','reconsSignalmVR_PsO'},...
        {'reconsSignalVRDOI','reconsSignalVRDOI','reconsSignalmVRDOI_PsO'},...
        {'reconsSignalVRptherm','reconsSignalVRptherm','reconsSignalmVRptherm_PsO'},...
        {'reconsSignalVRptherm','reconsSignalVRptherm','reconsSignalmVRptherm_PsO'},...
        {'reconsSignalVRptherm','reconsSignalVRptherm','reconsSignalmVRptherm_PsO'}};
    ftimes = {'timereconS','timereconS','timereconS_PsO'};
    idx_POIs = [1,3,7;1,3,1]; titles={'POI: 1','3','mean_{PsO}'};
    RMS_fields = {'RMS'};
    
    SNR_fields = {'DOI_noiseall'};
    suffix_extra = horzcat({'Based on V_{allSources}','Based on V_{DOI}'},...
        arrayfun(@(x) sprintf('Based on V_{allSources + N(0,%5.2f)µV}',x),NoiseAmp,'UniformOutput',false));
    
else
    %more RMSE and better SNR definitions
    fnamesall = {{'reconsSignalVR','reconsSignalVR','reconsSignalVR'},...
    {'reconsSignalVRDOIvibr','reconsSignalVRDOIvibr','reconsSignalVRDOIvibr'},...
    {'reconsSignalVRdvpov','reconsSignalVRdvpov','reconsSignalVRdvpov'},...
    {'reconsSignalVRdvpsn','reconsSignalVRdvpsn','reconsSignalVRdvpsn'},...
    {'reconsSignalVRptherm','reconsSignalVRptherm','reconsSignalVRptherm'},...
    {'reconsSignalVRptherm','reconsSignalVRptherm','reconsSignalVRptherm'},...
    {'reconsSignalVRptherm','reconsSignalVRptherm','reconsSignalVRptherm'}};
    ftimes = {'timereconS','timereconS','timereconS'};
    idx_POIs = [1,3,7;1,3,7]; titles={'POI: 1','3','mean_{PsO}'};
    RMS_fields = {'RMS','RMSDOIvr','RMSDOI','RMSDOIvrpOSCvr','RMSDOIvrpstatnoise'}; 
    RMS_ylabels = {'\bf RMSE','\bf RMSE DOI_{vibr}','\bf RMSE DOI',...
        '\bf RMSE DOI+OSC noise','\bf RMSE DOI + statnoise'};
    SNR_fields = {'DOIvr_NoiseAll','DOIvr_NoiseOSCvr','DOIvr_NoiseStatic'};
    SNR_ylabels = {'\bf  SNR','\bf  SNR_{vibr}','\bf  SNR_{static}'};
    suffix_extra = horzcat({'Based on V_{allSources}','Based on V_{DOIvibr}','Based on V_{OSCVibr}','Based on V_{DOI+statnoise}'},...
        arrayfun(@(x) sprintf('Based on V_{allSources + N(0,%5.2f)µV}',x),NoiseAmp,'UniformOutput',false));
end



for ifigs=1:length(fnamesall)
    f_recon = figure(20+ifigs);
    fnames = fnamesall{ifigs};
    
    if PSM_flag
        sf_recon(1) = subplot(liAmp,1+size(idx_POIs,2),1);
        fignr = get(gcf,'number');
        varargin2 = horzcat(Options,{'POI',POIs,'DOI',1,'OSC',1:Totaldps_plot,'scale',0});
        PotentialSphere_Multi(CSource_plot,CSink_plot,10,sigma,RSphere,SphereRes,fignr,varargin2);
        chPSM = get(gca,'children');
    end
    %rmpath(genpath('./Functions'));
    
    sf_recon(2+size(idx_POIs,2)) = subplot(liAmp,1+size(idx_POIs,2),2+size(idx_POIs,2));
    
    if nos4l_flag
        plot([Tsim(1:100)*1e6.*ones(3,1)]',USwave_all(Tsim(1:100)',1)*1e9)
        legend({'x component','y component','z component'})
        ylabel('displacement [nm]')
        xlabel('Time [µs]')
    else
        dAmpx = vampx(1)/wus;
        dAmpy = vampy(1)/wus;
        dAmpz = vampz(1)/wus;
        phasex = vphasex(1);
        phasey = vphasey(1);
        phasez = vphasez(1);
        genvibrplot(CSource(1,:),CSink(1,:),dAmpx,dAmpy,dAmpz,phasex,phasey,phasez,1)
    end
    title('Vibration')
    
    for iAmp = 1:liAmp
        liPOIs = size(idx_POIs,2);
        for iPOIs = 1:liPOIs
            rng(s);
            xvals = Outall(iAmps(iAmp),idpss(1)).Out.(ftimes{iPOIs}){idx_POIs(2,iPOIs)};
            if contains(fnames{iPOIs},'therm')
                % create reconstructed signals with new thermal noise
                % amplitudes
                if contains(fnames{iPOIs},'_PsO')
                    yval =  Outall(iAmps(iAmp),idpss(1)).Out.mVR_PsO;
                else
                    yval =  Outall(iAmps(iAmp),idpss(1)).Out.VR;
                end
                POIsidx = 1:size(yval,1);
                noise = NoiseAmp(3-(length(fnamesall)-ifigs))*randn(size(yval));
                yval = yval+noise;
                
                [reconsSignalVRptherm,powerSignalVRptherm,timerconsSptherm] = calcFourierGetsignal(...
                    yval,Tsim,POIsidx,fus,Window,inputSignal,mirror,theta_global,0);
                Outall(iAmps(iAmp),idpss(1)).Out.(fnames{iPOIs}) = reconsSignalVRptherm;
                xvals = timerconsSptherm{idx_POIs(2,iPOIs)};
            end
            
            % first row
            sf_idx = (iAmp-1)*(1+liPOIs)+1+iPOIs;
            sf_recon(sf_idx) = subplot(liAmp,1+size(idx_POIs,2),sf_idx);
            yyaxis left
            yvals = Outall(iAmps(iAmp),idpss(1)).Out.(fnames{iPOIs}){idx_POIs(2,iPOIs)};
            ratio = min(yvals)/max(yvals);
            if abs(ratio)>1; flip = -1;else flip = 1; end %during reconstruction phase of signal is determined except of pi
            yvals = flip*yvals;
            
            
            ratio = min(yvals)/max(yvals);
            plot(xvals,yvals)
            ylim([min(min(yvals),0),max(yvals)])
            yliml = get(sf_recon(sf_idx),'Ylim');
            xlim([0,Tend])
            if iPOIs==1
                ylabel('Reconstructed Signal [µV]')
            end
            yyaxis right
            plot(xvals,meanI.*Ifun(xvals));
            
            if iPOIs==size(idx_POIs,2)
                ylabel_str = {'Input Current [µA]',['{\color{black}\bf',sprintf('Vibration Amp = %5.2f µm}',Amp(iAmps(iAmp))*1e6)]};
                if nos4l_flag
                    ylabel_str = horzcat(ylabel_str,['{\color{black}\bf',sprintf('spatial distr k_{Aus} = %5.2f cm}',1/kAus(iAmps(iAmp))*100)]);
                end
                  ylabel(ylabel_str)  
            end
            title(titles{iPOIs})
            ylim([min(meanI.*Ifun(xvals)),max(meanI.*Ifun(xvals))])
            % correct axis to allign 0 of left and right y-axis
            ylimr = get(sf_recon(sf_idx),'Ylim');
            if ratio<0
                if ylimr(2)*ratio<ylimr(1)
                    set(sf_recon(sf_idx),'Ylim',[ylimr(2)*ratio ylimr(2)])
                else
                    set(sf_recon(sf_idx),'Ylim',[ylimr(1) ylimr(1)/ratio])
                end
            end
        end
    end
    set(gcf,'Position',figpos)
    legend_sp = get(sf_recon(1),'legend');
    legend_sp.Position = [0.0639,0.9012, 0.1093, 0.0810];
    set(sf_recon(1),'position',[0.0381,0.5368,0.1750,0.3833]);
    legend_sp = get(sf_recon(2+size(idx_POIs,2)),'legend');
    legend_sp.Position = [0.0153,0.0030,0.2728,0.1062];
    set(sf_recon(2+size(idx_POIs,2)),'position',[0.034343406179352,0.202494165316045,0.242220271288621,0.214037439222042]);
    mtit([title_prefix,title_suffix,sprintf(', # dipoles = %5.2e, ',Totaldps(idpss)),suffix_extra{ifigs}],'xoff',+0.17,'yoff',0.04);
end



%% Plot SNR and RMS figures
if Debug_flag
    Totaldps_plot = min(100,size(CSource,1));
    idx = [1,sort(randperm(size(CSource,1)-1,Totaldps_plot))+1];
    CSource_plot = CSource(idx,:);
    CSink_plot = CSink(idx,:);
    %OSCindices_plot = OSCindices(1:Totaldps_plot,:);
end
if PSM_flag
    figure(get(gcf,'number')+1);
    varargin2 = horzcat(Options,{'POI',POIs,'DOI',1,'OSC',1:Totaldps_plot,'scale',0});
    PotentialSphere_Multi(CSource_plot,CSink_plot,10,sigma,RSphere,SphereRes,get(gcf,'number'),varargin2);
    chPSM = get(gca,'children');
end
if ~any(Amp-Amp(1))
    if isfield(Input,'kAus')
        kAus_flag = 1;
        kAus = Input.kAus;
    else
        kAus_flag = 0;
        
    end
else
    kAus_flag = 0;
end
 y_SNR_idx = 1;

if length(RMS_fields) == 1 && length(SNR_fields)==1
    for iAmp = 1:lAmp
        y_RMS = [];
        y_SNR = [];
        for idps = 1:lTdps
            % gather values
            for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                y_RMS(iPOIs,idps) = RMSall(iAmp,idps).RMS(iPOIs).RMS;
                y_SNR(iPOIs,idps) = SNRall(iAmp,idps).SNR.DOI_noiseall(y_SNR_idx,iPOIs);
                if idps==1 && iAmp==1
                    POInames{iPOIs} = RMSall(iAmp,idps).RMS(iPOIs).info;
                    if strcmpi(POInames{iPOIs} ,'meanAllPOIs')
                        POInames{iPOIs} = 'mean_{All}';
                    end
                    if strcmpi(POInames{iPOIs} ,'meanPOIs_sameO')
                        POInames{iPOIs} = 'mean_{sameO}';
                    end
                end
            end
        end
        % plot RMS
        f_RMSSNR = figure(10);
        sf_RMSSNR(iAmp) = subplot(2,lAmp,iAmp);
        for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
            
            plot(Totaldps,y_RMS(iPOIs,:)','*-','DisplayName',POInames{iPOIs})
            if iPOIs==1; hold on; end
            if iPOIs == length(RMSall(iAmp,idps).RMS); hold off; end
        end
        set(gca,'xscale','log')
        %xlabel('number of dipoles')
        if iAmp==1
            if kAus_flag
                title(sprintf('spatial constant: k_{Aus} = %5.2f cm',1/kAus(iAmp)*100))
                ylabel('\bf RMSE wrt Input signal')
            else
                title(sprintf('Vibration Amp = %5.2f µm',Amp(iAmp)*1e6))
                ylabel('\bf RMSE wrt Input signal')
            end
            
        else
            if kAus_flag
                title(sprintf('%5.2f cm',1/kAus(iAmp)*100))
            else
                title(sprintf('%5.2f µm',Amp(iAmp)*1e6))
            end
        end
        if iAmp == lAmp
            lf_RMSSNR = legend('show','location','eastoutside');
        end
        
        % plot SNR
        sf_RMSSNR(iAmp+lAmp) = subplot(2,lAmp,iAmp+lAmp);
        plot(Totaldps,y_SNR','*-')
        set(gca,'xscale','log')
        if iAmp == 1
            ylabel('\bf  SNR')
        end
        xlabel('number of dipoles')
        
        
        
    end
    dl_pl = 0.05;
    figure(10)
    set(gcf,'position',figpos);
    mtit([title_prefix,title_suffix],'xoff',-0.05,'yoff',0.05);
    pause(0.01)
    set(findall(gcf,'-property','Fontsize'),'Fontsize',20)
    set(findobj(gcf,'type','axes'),'Fontsize',15)
    set(findobj(gcf,'type','line'),'Linewidth',2)
    % adjust positions of subplots and limits
    if lAmp>1
        max_x_f1 = sf_RMSSNR(lAmp).Position(1)+sf_RMSSNR(lAmp).Position(3);
        lsp = (max_x_f1-(lAmp-1)*dl_pl-sf_RMSSNR(1).Position(1))/lAmp;
        sf_RMSSNR(1).Position(3) = lsp;
        sf_RMSSNR(1+lAmp).Position(3) = lsp;
        
        YLIMS_RMS(1,:) = get(sf_RMSSNR(1),'ylim');
        YLIMS_SNR(1,:) = get(sf_RMSSNR(1+lAmp),'ylim');
        for iAmp = 2:lAmp
            sf_RMSSNR(iAmp).Position(1) = sf_RMSSNR(iAmp-1).Position(1)+lsp+dl_pl;
            sf_RMSSNR(iAmp+lAmp).Position(1) = sf_RMSSNR(iAmp+lAmp-1).Position(1)+lsp+dl_pl;
            sf_RMSSNR(iAmp).Position(3) = lsp;
            sf_RMSSNR(iAmp+lAmp).Position(3) = lsp;
            
            YLIMS_RMS(iAmp,:) = get(sf_RMSSNR(iAmp),'ylim');
            YLIMS_SNR(iAmp,:) = get(sf_RMSSNR(iAmp+lAmp),'ylim');
        end
        YLIMS_RMS = [min(YLIMS_RMS(:,1)),max(YLIMS_RMS(:,2))];
        YLIMS_SNR = [min(YLIMS_SNR(:,1)),max(YLIMS_SNR(:,2))];
        for iAmp = 1:lAmp
            set(sf_RMSSNR(iAmp),'ylim',YLIMS_RMS,'xtick',Totaldps(1:2:end),...
                'xtickLabel',cellstr(num2str(round(log10(Totaldps(1:2:end)')), '10^%d')))
            set(sf_RMSSNR(iAmp+lAmp),'ylim',YLIMS_SNR,'xtick',Totaldps(1:2:end),...
                'xtickLabel',cellstr(num2str(round(log10(Totaldps(1:2:end)')), '10^%d')))
            
        end
        lf_RMSSNR.Position(2) = 1/2*(sf_RMSSNR(1).Position(2)+sf_RMSSNR(1+lAmp).Position(2)+sf_RMSSNR(1+lAmp).Position(4))-lf_RMSSNR.Position(4)/2;
        lf_RMSSNR.Title.String = 'POI:';
    end
else
    nRMSfigs = ceil(length(RMS_fields)/3);
    nysb = 3*ones(1,nRMSfigs); nysb(end) = mod(length(RMS_fields)-1,3)+1;
    field_idx = 0;
    %plot RMS curves of multiple RMSe fields
    fignr = 9;
    for ifig = 1:nRMSfigs
        fignr = fignr+1;
    for iRMS = 1:nysb(ifig)
        field_idx = field_idx+1;
        for iAmp = 1:lAmp
            y_RMS = [];
            
            for idps = 1:lTdps
                % gather values
                for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                    y_RMS(iPOIs,idps) = RMSall(iAmp,idps).RMS(iPOIs).(RMS_fields{field_idx});
                    
                    if idps==1 && iAmp==1
                        POInames{iPOIs} = RMSall(iAmp,idps).RMS(iPOIs).info;
                        if strcmpi(POInames{iPOIs} ,'meanAllPOIs')
                            POInames{iPOIs} = 'mean_{All}';
                        end
                        if strcmpi(POInames{iPOIs} ,'meanPOIs_sameO')
                            POInames{iPOIs} = 'mean_{sameO}';
                        end
                    end
                end
            end
            % plot RMS
            
            f_RMS = figure(fignr);
            sf_RMS((iRMS-1)*lAmp+iAmp) = subplot(nysb(ifig),lAmp,(iRMS-1)*lAmp+iAmp);
            for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                
                plot(Totaldps,y_RMS(iPOIs,:)','*-','DisplayName',POInames{iPOIs})
                if iPOIs==1; hold on; end
                if iPOIs == length(RMSall(iAmp,idps).RMS); hold off; end
            end
            set(gca,'xscale','log')
            xlabel('number of dipoles')
            if iAmp==1
                ylabel(RMS_ylabels{field_idx})
                if iRMS==1
                    if kAus_flag
                        title(sprintf('spatial constant: k_{Aus} = %5.2f cm',1/kAus(iAmp)*100))
                    else
                        title(sprintf('Vibration Amp = %5.2f µm',Amp(iAmp)*1e6))
                    end
                end
                
            elseif iRMS==1
                if kAus_flag
                    title(sprintf('%5.2f cm',1/kAus(iAmp)*100))
                else
                    title(sprintf('%5.2f µm',Amp(iAmp)*1e6))
                end
            end
            if iAmp == lAmp && iRMS == 1
                lf_RMSSNR = legend('show','location','eastoutside');
            end

        end        
    end
    dl_pl = 0.05;
    figure(fignr)
    set(gcf,'position',figpos);
    mtit([title_prefix,title_suffix],'xoff',-0.05,'yoff',0.05);
    pause(0.01)
    set(findall(gcf,'-property','Fontsize'),'Fontsize',20)
    set(findobj(gcf,'type','axes'),'Fontsize',15)
    set(findobj(gcf,'type','line'),'Linewidth',2)
    % adjust positions of subplots and limits
    for iRMS = 1:nysb(ifig)
        YLIMS_RMS = [];
        if lAmp>1
            max_x_f1 = sf_RMS(lAmp).Position(1)+sf_RMS(lAmp).Position(3);
            lsp = (max_x_f1-(lAmp-1)*dl_pl-sf_RMS(1).Position(1))/lAmp;
            sf_RMS((iRMS-1)*lAmp+1).Position(3) = lsp;
            
            
            YLIMS_RMS(1,:) = get(sf_RMS((iRMS-1)*lAmp+1),'ylim');
            
            for iAmp = 2:lAmp
                sf_RMS((iRMS-1)*lAmp+iAmp).Position(1) = sf_RMS((iRMS-1)*lAmp+iAmp-1).Position(1)+lsp+dl_pl;
                sf_RMS((iRMS-1)*lAmp+iAmp).Position(3) = lsp;
  
                
                YLIMS_RMS(iAmp,:) = get(sf_RMS((iRMS-1)*lAmp+iAmp),'ylim');
            end
            YLIMS_RMS = [min(YLIMS_RMS(:,1)),max(YLIMS_RMS(:,2))];
            for iAmp = 1:lAmp
                set(sf_RMS((iRMS-1)*lAmp+iAmp),'ylim',YLIMS_RMS,'xtick',Totaldps(1:2:end),...
                    'xtickLabel',cellstr(num2str(round(log10(Totaldps(1:2:end)')), '10^%d')))
                
            end
        end
        lf_RMSSNR.Position(2) = 1/2*(sf_RMS(1).Position(2)+sf_RMS(1+lAmp).Position(2)+sf_RMS(1+lAmp).Position(4))-lf_RMSSNR.Position(4)/2;
        lf_RMSSNR.Title.String = 'POI:';
    end
    end
    %plot SNR curves
    nSNRfigs = ceil(length(SNR_fields)/3);
    nysb = 3*ones(1,nSNRfigs); nysb(end) = mod(length(SNR_fields)-1,3)+1;
    field_idx = 0;
    for ifig = 1:nSNRfigs
         fignr = fignr+1;
    for iSNR = 1:nysb(ifig)
        field_idx = field_idx+1;
        for iAmp = 1:lAmp            
            y_SNR = [];
            for idps = 1:lTdps
                % gather values
                for iPOIs = 1:length(RMSall(iAmp,idps).RMS)                    
                    y_SNR(iPOIs,idps) = SNRall(iAmp,idps).SNR.(SNR_fields{field_idx})(y_SNR_idx,iPOIs);
                    if idps==1 && iAmp==1
                        POInames{iPOIs} = RMSall(iAmp,idps).RMS(iPOIs).info;
                        if strcmpi(POInames{iPOIs} ,'meanAllPOIs')
                            POInames{iPOIs} = 'mean_{All}';
                        end
                        if strcmpi(POInames{iPOIs} ,'meanPOIs_sameO')
                            POInames{iPOIs} = 'mean_{sameO}';
                        end
                    end
                end
            end
            % plot RMS
            f_SNR = figure(fignr);
            sf_SNR((iSNR-1)*lAmp+iAmp) = subplot(nysb(ifig),lAmp,(iSNR-1)*lAmp+iAmp);
                   
            plot(Totaldps,y_SNR','*-')
            xlabel('number of dipoles')            
            set(gca,'xscale','log')
            
            if iAmp==1
                ylabel(SNR_ylabels{field_idx})
                if iSNR==1
                    if kAus_flag
                        title(sprintf('spatial constant: k_{Aus} = %5.2f cm',1/kAus(iAmp)*100))
                    else
                        title(sprintf('Vibration Amp = %5.2f µm',Amp(iAmp)*1e6))
                    end
                end
                
            elseif iSNR==1
                if kAus_flag
                    title(sprintf('%5.2f cm',1/kAus(iAmp)*100))
                else
                    title(sprintf('%5.2f µm',Amp(iAmp)*1e6))
                end
            end
            if iAmp == lAmp && iSNR == 1
                lf_RMSSNR = legend(POInames,'location','eastoutside');
            end
        end        
    end
    dl_pl = 0.05;
    figure(fignr)
    set(gcf,'position',figpos);
    mtit([title_prefix,title_suffix],'xoff',-0.05,'yoff',0.05);
    pause(0.01)
    set(findall(gcf,'-property','Fontsize'),'Fontsize',20)
    set(findobj(gcf,'type','axes'),'Fontsize',15)
    set(findobj(gcf,'type','line'),'Linewidth',2)
    for iSNR = 1:nysb(ifig)
        YLIMS_SNR = [];
        if lAmp>1
            max_x_f1 = sf_SNR(lAmp).Position(1)+sf_SNR(lAmp).Position(3);
            lsp = (max_x_f1-(lAmp-1)*dl_pl-sf_SNR(1).Position(1))/lAmp;
            sf_SNR((iSNR-1)*lAmp+1).Position(3) = lsp;            
            
            YLIMS_SNR(1,:) = get(sf_SNR((iSNR-1)*lAmp+1),'ylim');
            
            for iAmp = 2:lAmp
                sf_SNR((iSNR-1)*lAmp+iAmp).Position(1) = sf_SNR((iSNR-1)*lAmp+iAmp-1).Position(1)+lsp+dl_pl;
                sf_SNR((iSNR-1)*lAmp+iAmp).Position(3) = lsp;  
                
                YLIMS_SNR(iAmp,:) = get(sf_SNR((iSNR-1)*lAmp+iAmp),'ylim');
            end
            YLIMS_SNR = [min(YLIMS_SNR(:,1)),max(YLIMS_SNR(:,2))];

            for iAmp = 1:lAmp
                set(sf_SNR((iSNR-1)*lAmp+iAmp),'ylim',YLIMS_SNR,'xtick',Totaldps(1:2:end),...
                    'xtickLabel',cellstr(num2str(round(log10(Totaldps(1:2:end)')), '10^%d')))                
            end
        end
        lf_RMSSNR.Position(2) = 1/2*(sf_SNR(1).Position(2)+sf_SNR(1+lAmp).Position(2)+sf_SNR(1+lAmp).Position(4))-lf_RMSSNR.Position(4)/2;
        lf_RMSSNR.Title.String = 'POI:';
    end
    end
    
    
end

%% plot RMS & SNR vs AMP
if kAus_flag
    Amp = kAus;
end
idps_vals = [1,ceil(lTdps/2),lTdps];
if length(RMS_fields) == 1 && length(SNR_fields)==1
    for idps = idps_vals
        idps_pos = find(idps==idps_vals);
        y_RMS = [];
        y_SNR = [];
        for iAmp = 1:lAmp
            % gather values
            for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                y_RMS(iPOIs,iAmp) = RMSall(iAmp,idps).RMS(iPOIs).RMS;
                y_SNR(iPOIs,iAmp) = SNRall(iAmp,idps).SNR.DOI_noiseall(y_SNR_idx,iPOIs);
                if idps==1 && iAmp==1
                    POInames{iPOIs} = RMSall(iAmp,idps).RMS(iPOIs).info;
                    if strcmpi(POInames{iPOIs} ,'meanAllPOIs')
                        POInames{iPOIs} = 'mean_{All}';
                    end
                    if strcmpi(POInames{iPOIs} ,'meanPOIs_sameO')
                        POInames{iPOIs} = 'mean_{sameO}';
                    end
                end
            end
        end
        % plot RMS
        f_RMSSNR_Amp = figure(11);
        sf_RMSSNR_Amp(idps_pos) = subplot(2,length(idps_vals),idps_pos);
        for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
            
            plot(Amp,y_RMS(iPOIs,:)','*-','DisplayName',POInames{iPOIs})
            if iPOIs==1; hold on; end
            if iPOIs == length(RMSall(iAmp,idps).RMS)
                hold off;
            end
        end
        set(gca,'xscale','log')
        %xlabel('number of dipoles')
        if idps==1
            title(sprintf('Vibration with %5.2e dipoles',Totaldps(idps)))
            ylabel('\bf RMSE wrt Input signal')
        else
            title(sprintf('%5.2e dipoles',Totaldps(idps)))
        end
        if idps == idps_vals(end)
            lf_RMSSNR_Amp = legend('show','location','eastoutside');
        end
        
        % plot SNR
        sf_RMSSNR_Amp(idps_pos+length(idps_vals)) = subplot(2,length(idps_vals),idps_pos+length(idps_vals));
        plot(Amp,y_SNR,'*-')
        set(gca,'xscale','log')
        if idps_pos == 1
            ylabel('\bf  SNR')
        end
        if kAus_flag
            xlabel('spatial constant: k_{Aus} [cm]')
            Xtick = get(gca,'xtick');
                    Xticknew = 1./Xtick*100;
                    Xticklabelnew = arrayfun(@num2str,Xticknew,'UniformOutput',false);
                    set(gca,'xticklabel',Xticklabelnew)
        else
            xlabel('Vibration Amplitude DOI [µm]')
        end
        
        
        
    end
    dl_pl = 0.05;
    figure(11)
    set(gcf,'position',figpos);
    mtit([title_prefix,title_suffix],'xoff',-0.05,'yoff',0.05);
    pause(0.1)
    set(findall(gcf,'-property','Fontsize'),'Fontsize',20)
    set(findobj(gcf,'type','axes'),'Fontsize',15)
    set(findobj(gcf,'type','line'),'Linewidth',2)
    
    if lTdps>1
        max_x_f1 = sf_RMSSNR_Amp(lAmp).Position(1)+sf_RMSSNR_Amp(lAmp).Position(3);
        lsp = (max_x_f1-(lAmp-1)*dl_pl-sf_RMSSNR_Amp(1).Position(1))/lAmp;
        sf_RMSSNR_Amp(1).Position(3) = lsp;
        sf_RMSSNR_Amp(1+lAmp).Position(3) = lsp;
        
        YLIMS_RMS(1,:) = get(sf_RMSSNR_Amp(1),'ylim');
        YLIMS_SNR(1,:) = get(sf_RMSSNR_Amp(1+lAmp),'ylim');
        for iAmp = 2:lAmp
            sf_RMSSNR_Amp(iAmp).Position(1) = sf_RMSSNR_Amp(iAmp-1).Position(1)+lsp+dl_pl;
            sf_RMSSNR_Amp(iAmp+lAmp).Position(1) = sf_RMSSNR_Amp(iAmp+lAmp-1).Position(1)+lsp+dl_pl;
            sf_RMSSNR_Amp(iAmp).Position(3) = lsp;
            sf_RMSSNR_Amp(iAmp+lAmp).Position(3) = lsp;
            
            YLIMS_RMS(iAmp,:) = get(sf_RMSSNR_Amp(iAmp),'ylim');
            YLIMS_SNR(iAmp,:) = get(sf_RMSSNR_Amp(iAmp+lAmp),'ylim');
        end
        YLIMS_RMS = [min(YLIMS_RMS(:,1)),max(YLIMS_RMS(:,2))];
        YLIMS_SNR = [min(YLIMS_SNR(:,1)),max(YLIMS_SNR(:,2))];
        for idps = 1:length(idps_vals)
            if ~kAus_flag
                set(sf_RMSSNR_Amp(idps),'ylim',YLIMS_RMS,'xtick',Amp,...
                    'xtickLabel',arrayfun(@(x) sprintf('%5.2f',x),Amp*1e6,'UniformOutput',false))
                set(sf_RMSSNR_Amp(idps+length(idps_vals)),'ylim',YLIMS_SNR,'xtick',Amp,...
                    'xtickLabel',arrayfun(@(x) sprintf('%5.2f',x),Amp*1e6,'UniformOutput',false))
            else
                set(sf_RMSSNR_Amp(idps),'ylim',YLIMS_RMS,'xtick',Amp,...
                    'xtickLabel',arrayfun(@(x) sprintf('%5f',x),Amp,'UniformOutput',false))                
                    Xtick = get(sf_RMSSNR_Amp(idps),'xtick');
                    Xticknew = 1./Xtick*100;
                    Xticklabelnew = arrayfun(@num2str,Xticknew,'UniformOutput',false);
                    set(sf_RMSSNR_Amp(idps),'xticklabel',Xticklabelnew)
                
                set(sf_RMSSNR_Amp(idps+length(idps_vals)),'ylim',YLIMS_SNR,'xtick',Amp,...
                    'xtickLabel',arrayfun(@(x) sprintf('%5f',x),Amp,'UniformOutput',false))
                    Xtick = get(sf_RMSSNR_Amp(idps+length(idps_vals)),'xtick');
                    Xticknew = 1./Xtick*100;
                    Xticklabelnew = arrayfun(@num2str,Xticknew,'UniformOutput',false);
                    set(sf_RMSSNR_Amp(idps+length(idps_vals)),'xticklabel',Xticklabelnew)
            end
            
            
        end
        lf_RMSSNR_Amp.Position(2) = 1/2*(sf_RMSSNR_Amp(1).Position(2)+sf_RMSSNR_Amp(1+lAmp).Position(2)+sf_RMSSNR_Amp(1+lAmp).Position(4))-lf_RMSSNR_Amp.Position(4)/2;
        lf_RMSSNR_Amp.Title.String = 'POI:';
    end
else
    nRMSfigs = ceil(length(RMS_fields)/3);
    nysb = 3*ones(1,nRMSfigs); nysb(end) = mod(length(RMS_fields)-1,3)+1;
    field_idx = 0;
    for ifig = 1:nRMSfigs
        fignr = fignr+1;
    for iRMS = 1:nysb(ifig)
        field_idx = field_idx+1;
        for idps = idps_vals
            idps_pos = find(idps==idps_vals);
            y_RMS = [];
            for iAmp = 1:lAmp
                % gather values
                for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                    y_RMS(iPOIs,iAmp) = RMSall(iAmp,idps).RMS(iPOIs).(RMS_fields{field_idx});
                    if idps==1 && iAmp==1
                        POInames{iPOIs} = RMSall(iAmp,idps).RMS(iPOIs).info;
                        if strcmpi(POInames{iPOIs} ,'meanAllPOIs')
                            POInames{iPOIs} = 'mean_{All}';
                        end
                        if strcmpi(POInames{iPOIs} ,'meanPOIs_sameO')
                            POInames{iPOIs} = 'mean_{sameO}';
                        end
                    end
                end
            end
            % plot RMS
            f_RMS_Amp = figure(fignr);
            lidps = length(idps_vals);
            sf_RMS_Amp((iRMS-1)*lidps+idps_pos) = subplot(nysb(ifig),length(idps_vals),(iRMS-1)*lidps+idps_pos);
            for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                plot(Amp,y_RMS(iPOIs,:)','*-','DisplayName',POInames{iPOIs})
                if iPOIs==1; hold on; end
                if iPOIs == length(RMSall(iAmp,idps).RMS); hold off; end
            end
            
            set(gca,'xscale','log')
            if kAus_flag
                Xtick = get(gca,'xtick');
                Xticknew = 1./Xtick*100;
                Xticklabelnew = arrayfun(@num2str,Xticknew,'UniformOutput',false);
                set(gca,'xticklabel',Xticklabelnew)
            end
            
            if find(idps==idps_vals)==1
                ylabel(RMS_ylabels{field_idx}) 
            end
                
            if iRMS==1 &&  find(idps==idps_vals)==1           
                title(sprintf('Vibration with %5.2e dipoles',Totaldps(idps)))
            elseif iRMS == 1
                title(sprintf('%5.2e dipoles',Totaldps(idps)))
            end
            if idps == idps_vals(end) && iRMS == 1
                lf_RMS_Amp = legend('show','location','eastoutside');
            end
            if iRMS == nysb(ifig)
                if kAus_flag
                    xlabel('spatial constant: k [cm]')
                else
                    xlabel('Vibration Amplitude DOI [µm]')
                end
            end            
        end
    end
    dl_pl = 0.05;
    figure(fignr)
    set(gcf,'position',figpos);
    mtit([title_prefix,title_suffix],'xoff',-0.05,'yoff',0.05);
    pause(0.1)
    set(findall(gcf,'-property','Fontsize'),'Fontsize',20)
    set(findobj(gcf,'type','axes'),'Fontsize',15)
    set(findobj(gcf,'type','line'),'Linewidth',2)
    for iRMS=1:nysb(ifig)
        YLIMS_RMS = [];
        if lTdps>1
            max_x_f1 = sf_RMS_Amp(lidps).Position(1)+sf_RMS_Amp(lidps).Position(3);
            lsp = (max_x_f1-(lidps-1)*dl_pl-sf_RMS_Amp(1).Position(1))/lidps;
            sf_RMS_Amp((iRMS-1)*lidps+1).Position(3) = lsp;
            
            YLIMS_RMS(1,:) = get(sf_RMS_Amp((iRMS-1)*lidps+1),'ylim');
            for idps_pos = 2:lidps
                sf_RMS_Amp((iRMS-1)*lidps+idps_pos).Position(1) = sf_RMS_Amp((iRMS-1)*lidps+idps_pos-1).Position(1)+lsp+dl_pl; 
                sf_RMS_Amp((iRMS-1)*lidps+idps_pos).Position(3) = lsp;
                
                
                YLIMS_RMS(idps_pos,:) = get(sf_RMS_Amp((iRMS-1)*lidps+idps_pos),'ylim');                
            end
            YLIMS_RMS = [min(YLIMS_RMS(:,1)),max(YLIMS_RMS(:,2))];
            for idps = 1:length(idps_vals)
                if ~kAus_flag
                    set(sf_RMS_Amp((iRMS-1)*lidps+idps),'ylim',YLIMS_RMS,'xtick',Amp,...
                        'xtickLabel',arrayfun(@(x) sprintf('%5.2f',x),Amp*1e6,'UniformOutput',false))                  
                else
                    set(sf_RMS_Amp((iRMS-1)*lidps+idps),'ylim',YLIMS_RMS,'xtick',Amp,...
                        'xtickLabel',arrayfun(@(x) sprintf('%5.0f',x),Amp,'UniformOutput',false))
                    Xtick = get(sf_RMS_Amp((iRMS-1)*lidps+idps),'xtick');
                    Xticknew = 1./Xtick*100;
                    Xticklabelnew = arrayfun(@num2str,Xticknew,'UniformOutput',false);
                    set(sf_RMS_Amp((iRMS-1)*lidps+idps),'xticklabel',Xticklabelnew)

                end
                
                
            end
            lf_RMS_Amp.Position(2) = 1/2*(sf_RMS_Amp(1).Position(2)+sf_RMS_Amp(1+lidps).Position(2)+sf_RMS_Amp(1+lidps).Position(4))-lf_RMS_Amp.Position(4)/2;
            lf_RMS_Amp.Title.String = 'POI:';
        end
    end
    end
    
    
    %plot SNR curves
    nSNRfigs = ceil(length(SNR_fields)/3);
    nysb = 3*ones(1,nSNRfigs); nysb(end) = mod(length(SNR_fields)-1,3)+1;
    field_idx = 0;
    for ifig = 1:nSNRfigs
         fignr = fignr+1;
    for iSNR = 1:nysb(ifig)
        field_idx = field_idx+1;
        for idps = idps_vals
            idps_pos = find(idps==idps_vals);
            y_SNR = [];
            for iAmp = 1:lAmp
                % gather values
                for iPOIs = 1:length(RMSall(iAmp,idps).RMS)
                    y_SNR(iPOIs,iAmp) = SNRall(iAmp,idps).SNR.(SNR_fields{iSNR})(y_SNR_idx,iPOIs);
                    if idps==1 && iAmp==1
                        POInames{iPOIs} = RMSall(iAmp,idps).RMS(iPOIs).info;
                        if strcmpi(POInames{iPOIs} ,'meanAllPOIs')
                            POInames{iPOIs} = 'mean_{All}';
                        end
                        if strcmpi(POInames{iPOIs} ,'meanPOIs_sameO')
                            POInames{iPOIs} = 'mean_{sameO}';
                        end
                    end
                end
            end
            % plot RMS
            f_SNR_Amp = figure(fignr);
            lidps = length(idps_vals);
            sf_SNR_Amp((iSNR-1)*lidps+idps_pos) = subplot(nysb(ifig),length(idps_vals),(iSNR-1)*lidps+idps_pos);
            
            % plot SNR
            plot(Amp,y_SNR,'*-')
            set(gca,'xscale','log')
            
            if find(idps==idps_vals)==1
                ylabel(SNR_ylabels{field_idx})
            end
            if find(idps==idps_vals)==1 && iSNR==1
                title(sprintf('Vibration with %5.2e dipoles',Totaldps(idps)))
            elseif iSNR==1                
                title(sprintf('%5.2e dipoles',Totaldps(idps)))
            end
            if idps == idps_vals(end) && iSNR == 1
                lf_SNR_Amp = legend(POInames,'location','eastoutside');
            end
            if iSNR == nysb(ifig)
                if kAus_flag
                    xlabel('spatial constant: k [cm]')
                else
                    xlabel('Vibration Amplitude DOI [µm]')
                end
            end
        end
    end
    dl_pl = 0.05;
    figure(fignr)
    set(gcf,'position',figpos);
    mtit([title_prefix,title_suffix],'xoff',-0.05,'yoff',0.05);
    pause(0.1)
    set(findall(gcf,'-property','Fontsize'),'Fontsize',20)
    set(findobj(gcf,'type','axes'),'Fontsize',15)
    set(findobj(gcf,'type','line'),'Linewidth',2)
    for iSNR=1:length(SNR_fields)
        YLIMS_SNR = [];
        if lTdps>1
            max_x_f1 = sf_SNR_Amp(lidps).Position(1)+sf_SNR_Amp(lidps).Position(3);
            lsp = (max_x_f1-(lidps-1)*dl_pl-sf_SNR_Amp(1).Position(1))/lidps;
            sf_SNR_Amp((iSNR-1)*lidps+1).Position(3) = lsp;
            
            YLIMS_SNR(1,:) = get(sf_SNR_Amp((iSNR-1)*lidps+1),'ylim');
            
            for idps_pos = 2:lidps
                sf_SNR_Amp((iSNR-1)*lidps+idps_pos).Position(1) = sf_SNR_Amp((iSNR-1)*lidps+idps_pos-1).Position(1)+lsp+dl_pl;
                sf_SNR_Amp((iSNR-1)*lidps+idps_pos).Position(3) = lsp;
                
                
                YLIMS_SNR(idps_pos,:) = get(sf_SNR_Amp((iRMS-1)*lidps+idps_pos),'ylim');
            end
            YLIMS_SNR = [min(YLIMS_SNR(:,1)),max(YLIMS_SNR(:,2))];
            for idps = 1:length(idps_vals)
                if ~kAus_flag
                    set(sf_SNR_Amp((iSNR-1)*lidps+idps),'ylim',YLIMS_SNR,'xtick',Amp,...
                        'xtickLabel',arrayfun(@(x) sprintf('%5.2f',x),Amp*1e6,'UniformOutput',false))

                else
                    set(sf_SNR_Amp((iSNR-1)*lidps+idps),'ylim',YLIMS_SNR,'xtick',Amp,...
                        'xtickLabel',arrayfun(@(x) sprintf('%5.0f',x),Amp,'UniformOutput',false))
                    Xtick = get(sf_SNR_Amp((iSNR-1)*lidps+idps),'xtick');
                    Xticknew = 1./Xtick*100;
                    Xticklabelnew = arrayfun(@num2str,Xticknew,'UniformOutput',false);
                    set(sf_SNR_Amp((iSNR-1)*lidps+idps),'xticklabel',Xticklabelnew)
                end
                
                
            end
            lf_SNR_Amp.Position(2) = 1/2*(sf_SNR_Amp(1).Position(2)+sf_SNR_Amp(1+lidps).Position(2)+sf_SNR_Amp(1+lidps).Position(4))-lf_SNR_Amp.Position(4)/2;
            lf_SNR_Amp.Title.String = 'POI:';
        end
    end
    end
end

end