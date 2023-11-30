% process_DBSTUDY
% create a lot of plots of different metrics. not best and clean code (a
% lot of copy paste but works)
if ~isdeployed
    clear all
    close all
    clc
end
lim_totaldip_flag = 1;    %flag on limitation total dipoles
lim_Amp_flag = 1;         %flag on limitation Amplitudes
PSM_flag = 1;             %flag plot potential sphere multi sources
Debug_flag = 1;           %debug flag: minimize dipoles plot
figsprun = 12;
Old_dir = cd('D:\no backup\EEGUS\DBStudy\Results');
lall = dir('*.mat');
cd(Old_dir);
fPOIstp_flag =1;
fPOIstp = [1:9];          %POIs to plot
Inputnames = {'Input_Mosc_sM_v2_1','Input_Mosc_sM_v2_2','Input_Mosc_sM_v2_3','Input_Mosc_sM_v2_5'};
Colors= [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    0 1 1
    1 0 1
    1 1 1
    ];


line_types = {'-','--','-.',':'};
Amp = logspace(2,4,3).*1e-9;
if lim_Amp_flag
    fprintf('\n !limitation on # Amps!\n')
    Amp(2:end-1) = [];
end

Totaldps = round(logspace(1,5,9));
names = {'Deep single layer','Cortex single layer','Deep 10 layers','Deep Hotspots'};
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
    Totaldps(i) = round(Totaldps(i),orders(i));
end
if lim_totaldip_flag
    fprintf('\n !limitation on # dipoles!')
    Totaldps = logspace(1,4,7);
    orders = -floor(log10(Totaldps));
    for i=1:length(Totaldps)
        Totaldps(i) = round(Totaldps(i),orders(i));
    end
    dps_run = min(100*ones(size(Totaldps)),Totaldps);
end

iInmax = length(lall);
for iIn = 1:iInmax
    RMSREVR = [];
    RMSREVRDOI = [];
    RMSREVRptherm=[];
    RMSREreconSignalVR = [];
    RMSEVR = [];
    RMSEVRDOI = [];
    RMSEVRptherm=[];
    RMSEreconSignalVR = [];
    nRMSEVR = [];
    nRMSEVRDOI = [];
    nRMSEVRptherm=[];
    nRMSEreconSignalVR = [];
    nmmRMSEVR = [];
    nmmRMSEVRDOI = [];
    nmmRMSEVRptherm=[];
    nmmRMSEreconSignalVR = [];
    timevals_n = [];
    timevals_DB = [];
    %[RMSall,SNRall,Outall,TSTOPall,Totaldps] = calcStrengthCurve_mosc_DBfvsnormal(Inputnames{iIn});
    ftl = fullfile(lall(iIn).folder,lall(iIn).name);
    fprintf('\nloading: %s\n',ftl);
    load(ftl);
    
    
    
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
        case 'mouse_fair0'
            SolutionType = 'Mouse4Sphere~fair0';
            RPOI = 0.0059;
            RBrain = 0.0046;
            dRBrain = 0.0005;
            title_prefix = 'Measure: Mouse Scalp,  ';
        otherwise
            error('false input model')
    end
    sigma = 0.33;
    SphereRes = 30;
    [Options,RSphere] = getSettings(SolutionType,1);
    Options{find(strcmpi(Options,'tAir'))+1}=0;
    if strcmpi(SolutionType,'Mouse4Sphere~fair0') || strcmpi(SolutionType,'Mouse4Sphere~fair1')
        dRSphere = 0.0005;
    else
        dRSphere = 0.005;
    end
    
    CSource = Outall(end,end).Out.CSource;
    CSink = Outall(end,end).Out.CSink;
    OSCindices = Outall(end,end).Out.OSCindices;
    POIs = Outall(end,end).Out.POIs;
    if Debug_flag
        Totaldps_plot = 1;
        dps_selec = round(linspace(1,size(CSource,1),Totaldps_plot));
        dps_selec = 1:Totaldps_plot;
        CSource_plot = CSource(dps_selec,:);
        CSink_plot = CSink(dps_selec,:);
        OSCindices_plot = OSCindices(dps_selec,:);
    end
    if PSM_flag
        figure(100+iIn);
        varargin2 = horzcat(Options,{'POI',POIs,'DOI',1,'OSC',1:Totaldps_plot,'scale',0});
        PotentialSphere_Multi(CSource_plot,CSink_plot,10,sigma,RSphere,SphereRes,get(gcf,'number'),varargin2);
        chPSM = get(gca,'children');
    end
    
    
    
    if fPOIstp_flag
        POIstp = fPOIstp;
    else
        POIstp = length(Outall(1,idps).Out.reconsSignalVR);
    end
    
    for idps = 1:length(Totaldps)
        timevals_n(idps) = TSTOPall(1,idps).TSTOP.endsim-TSTOPall(1,idps).TSTOP.Itgen(1);
        timevals_DB(idps) = TSTOPall(1,idps).TSTOP_DB.endsim-TSTOPall(1,idps).TSTOP_DB.Itgen(1);
    end
    %plot time difference normal method and database method slide 18
    %extension 200106
    f1 = figure(1);
    diff_tval_nDB = timevals_n-timevals_DB;
    sf1(1) = subplot(2,2,1);
    hold on
    plot(Totaldps,timevals_n,line_types{1},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, normal',names{iIn}));
    plot(Totaldps,timevals_DB,line_types{2},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, database',names{iIn}));
    hold off
    if iIn==iInmax
        xlabel('number of dipoles')
        ylabel('simulation time (s)')
        set(gca,'xscale','log')
        %legend('show')
        title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
    end
    
    sf1(3) = subplot(2,2,3);
    hold on
    plot(Totaldps,diff_tval_nDB,line_types{1},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, t_n-t_{db}',names{iIn}));
    hold off
    if iIn==iInmax
        xlabel('number of dipoles')
        ylabel('simulation time (s)')
        set(gca,'xscale','log')
        %legend('show','location','eastoutside')
        title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
    end
    
    sf1(2) = subplot(2,2,2);
    for idps = 1:length(Totaldps)
        timevals_n(idps) = TSTOPall(end,idps).TSTOP.endsim-TSTOPall(end,idps).TSTOP.Itgen(1);
        timevals_DB(idps) = TSTOPall(end,idps).TSTOP_DB.endsim-TSTOPall(end,idps).TSTOP_DB.Itgen(1);
    end
    diff_tval_nDB = timevals_n-timevals_DB;
    hold on
    plot(Totaldps,timevals_n,line_types{1},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, normal',names{iIn}));
    plot(Totaldps,timevals_DB,line_types{2},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, database',names{iIn}));
    hold off
    if iIn==iInmax
        xlabel('number of dipoles')
        ylabel('simulation time (s)')
        set(gca,'xscale','log')
        legend('show','location','eastoutside')
        title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
    end
    
    sf1(4) = subplot(2,2,4);
    hold on
    plot(Totaldps,diff_tval_nDB,line_types{1},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, t_n-t_{db}',names{iIn}));
    hold off
    if iIn==iInmax
        xlabel('number of dipoles')
        ylabel('simulation time (s)')
        set(gca,'xscale','log')
        legend('show','location','eastoutside')
        title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
    end
    
    %get different error metrics: RMSE of (VR,VRDOI,VRpthermal noise, reconstructed signal)...
    %normalized errors R values RMSRE nromalized with difference max min or norm with maxed 
    for iAmp=1:length(Amp)
        for idps = 1:length(Totaldps)
            RMSREVR(:,idps) = sqrt(nanmean(((Outall(iAmp,idps).Out.VR-Outall(iAmp,idps).Out_DB.VR)./Outall(iAmp,idps).Out.VR).^2,2));
            nRMSEVR(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VR-Outall(iAmp,idps).Out_DB.VR)).^2,2))./mean(Outall(iAmp,idps).Out.VR,2);
            meanVR = mean(Outall(iAmp,idps).Out.VR,2);
            RsqVR(:,idps) = 1-sum((Outall(iAmp,idps).Out.VR-Outall(iAmp,idps).Out_DB.VR).^2,2)./sum((Outall(iAmp,idps).Out.VR-meanVR).^2,2);
            nmmRMSEVR(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VR-Outall(iAmp,idps).Out_DB.VR)).^2,2))./(max(Outall(iAmp,idps).Out.VR,[],2)-min(Outall(iAmp,idps).Out.VR,[],2));
            RMSEVR(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VR-Outall(iAmp,idps).Out_DB.VR)).^2,2));
            nmRMSEVR(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VR-Outall(iAmp,idps).Out_DB.VR)).^2,2))./(max(Outall(iAmp,idps).Out.VR,[],2));
            
            RMSREVRDOI(:,idps) = sqrt(nanmean(((Outall(iAmp,idps).Out.VRDOI-Outall(iAmp,idps).Out_DB.VRDOI)./Outall(iAmp,idps).Out.VRDOI).^2,2));
            nRMSEVRDOI(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRDOI-Outall(iAmp,idps).Out_DB.VRDOI)).^2,2))./mean(Outall(iAmp,idps).Out.VRDOI,2);
            meanVRDOI = mean(Outall(iAmp,idps).Out.VRDOI,2);
            RsqVRDOI(:,idps) = 1-sum((Outall(iAmp,idps).Out.VRDOI-Outall(iAmp,idps).Out_DB.VRDOI).^2,2)./sum((Outall(iAmp,idps).Out.VRDOI-meanVRDOI).^2,2);
            nmmRMSEVRDOI(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRDOI-Outall(iAmp,idps).Out_DB.VRDOI)).^2,2))./(max(Outall(iAmp,idps).Out.VRDOI,[],2)-min(Outall(iAmp,idps).Out.VRDOI,[],2));
            RMSEVRDOI(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRDOI-Outall(iAmp,idps).Out_DB.VRDOI)).^2,2));
            nmRMSEVRDOI(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRDOI-Outall(iAmp,idps).Out_DB.VRDOI)).^2,2))./(max(Outall(iAmp,idps).Out.VRDOI,[],2));
            
            
            RMSREVRptherm(:,idps) = sqrt(nanmean(((Outall(iAmp,idps).Out.VRptherm-Outall(iAmp,idps).Out_DB.VRptherm)./Outall(iAmp,idps).Out.VRptherm).^2,2));
            nRMSEVRptherm(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRptherm-Outall(iAmp,idps).Out_DB.VRptherm)).^2,2))./mean(Outall(iAmp,idps).Out.VRptherm,2);            
            meanVRptherm = mean(Outall(iAmp,idps).Out.VRptherm,2);
            RsqVRptherm(:,idps) = 1-sum((Outall(iAmp,idps).Out.VRptherm-Outall(iAmp,idps).Out_DB.VRptherm).^2,2)./sum((Outall(iAmp,idps).Out.VRptherm-meanVRptherm).^2,2);
            nmmRMSEVRptherm(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRptherm-Outall(iAmp,idps).Out_DB.VRptherm)).^2,2))./(max(Outall(iAmp,idps).Out.VRptherm,[],2)-min(Outall(iAmp,idps).Out.VRptherm,[],2));
            RMSEVRptherm(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRptherm-Outall(iAmp,idps).Out_DB.VRptherm)).^2,2));
            nmRMSEVRptherm(:,idps) = sqrt(mean(((Outall(iAmp,idps).Out.VRptherm-Outall(iAmp,idps).Out_DB.VRptherm)).^2,2))./(max(Outall(iAmp,idps).Out.VRptherm,[],2));
            
            
            for iPOI = POIstp
                RMSREreconSignalVR(iPOI,idps) = sqrt(nanmean(((Outall(iAmp,idps).Out.reconsSignalVR{iPOI}-Outall(iAmp,idps).Out_DB.reconsSignalVR{iPOI})./Outall(iAmp,idps).Out.reconsSignalVR{iPOI}).^2,2));
                nRMSEreconSignalVR(iPOI,idps) = sqrt(mean(((Outall(iAmp,idps).Out.reconsSignalVR{iPOI}-Outall(iAmp,idps).Out_DB.reconsSignalVR{iPOI})).^2,2))./mean(Outall(iAmp,idps).Out.reconsSignalVR{iPOI},2);
                meanVRreconS(iPOI,idps) = mean(Outall(iAmp,idps).Out.reconsSignalVR{iPOI},2);
                RsqVRreconS(iPOI,idps) = 1-sum((Outall(iAmp,idps).Out.reconsSignalVR{iPOI}-Outall(iAmp,idps).Out_DB.reconsSignalVR{iPOI}).^2,2)./sum((Outall(iAmp,idps).Out.reconsSignalVR{iPOI}-meanVRreconS(iPOI,idps)).^2,2);
                nmmRMSEreconS(iPOI,idps) = sqrt(mean(((Outall(iAmp,idps).Out.reconsSignalVR{iPOI}-Outall(iAmp,idps).Out_DB.reconsSignalVR{iPOI})).^2,2))./(max(Outall(iAmp,idps).Out.reconsSignalVR{iPOI},[],2)-min(Outall(iAmp,idps).Out.reconsSignalVR{iPOI},[],2));
                RMSEreconS(iPOI,idps) = sqrt(mean(((Outall(iAmp,idps).Out.reconsSignalVR{iPOI}-Outall(iAmp,idps).Out_DB.reconsSignalVR{iPOI})).^2,2));
                nmRMSEreconS(iPOI,idps) = sqrt(mean(((Outall(iAmp,idps).Out.reconsSignalVR{iPOI}-Outall(iAmp,idps).Out_DB.reconsSignalVR{iPOI})).^2,2))./(max(Outall(iAmp,idps).Out.reconsSignalVR{iPOI},[],2));
            end
        end
        
        
        % figure 2 using RMSRE: plot see slide 19 raportations extension
        % 200106
        f2(iAmp).f = figure(2+(iAmp-1)*figsprun);
        if iIn<=2
            
            f2(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn == iInmax
                xlabel('number of dipoles')
                ylabel('RMSRE Potential')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f2(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSRE Potential DOI only')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f2(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                xlabel('number of dipoles')
                ylabel('RMSRE Potential')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f2(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSRE Potential DOI only')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        
        
        % figure 3 using RMSRE ptherm en recon
        f3(iAmp).f = figure(3+(iAmp-1)*figsprun);
        if iIn<=2
            
            f3(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSRE Potential + Noise_{therm}')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f3(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSRE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f3(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('RMSRE Potential + Noise_{therm}')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f3(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSREreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSRE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        
        % figure 4 using nRMSE vr en doi
        f4(iAmp).f = figure(4+(iAmp-1)*figsprun);
        if iIn<=2
            
            f4(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('nRMSE Potential')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f4(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('nRMSE Potential DOI only')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f4(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                xlabel('number of dipoles')
                ylabel('nRMSE Potential')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f4(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('nRMSE Potential DOI only')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        
        
        % figure 5 nRMSE  recon ptherm
        f5(iAmp).f = figure(5+(iAmp-1)*figsprun);
        if iIn<=2
            
            f5(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('nRMSE Potential + Noise_{therm}')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f5(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('nRMSE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f5(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('nRMSE Potential + Noise_{therm}')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f5(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nRMSEreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('nRMSE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        pause(1)
        % figure 6 R^2 doi VR
        f6(iAmp).f = figure(6+(iAmp-1)*figsprun);
        if iIn<=2
            
            f6(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('R^2 Potential')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f6(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('R^2 VR DOI')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f6(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('R^2 Potential')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f6(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('R^2 Potenetial DOI')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        % figure7 R2 ptherm and recon
        f7(iAmp).f = figure(7+(iAmp-1)*figsprun);
        if iIn<=2
            
            f7(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('R^2 Potential + Noise_{therm}')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f7(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVRreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('R^2 Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f7(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('R^2 Potential + Noise_{therm}')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f7(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RsqVRreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('R^2 Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        
        % figure 8 nmmRMSE doi VR
        f8(iAmp).f = figure(8+(iAmp-1)*figsprun);
        if iIn<=2
            
            f8(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE   Potential')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f8(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE VR DOI')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f8(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE Potential')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f8(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE Potenetial DOI')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        % figure9 nmmRMSE ptherm and recon
        f9(iAmp).f = figure(9+(iAmp-1)*figsprun);
        if iIn<=2
            
            f9(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE Potential + Noise_{therm}')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f9(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f9(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE Potential + Noise_{therm}')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f9(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmmRMSEreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max-min}RMSE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        
        % figure 10 RMSE doi VR
        f10(iAmp).f = figure(10+(iAmp-1)*figsprun);
        if iIn<=2
            
            f10(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSE Potential [µV]')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f10(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSE VR DOI [µV]')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f10(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('RMSE Potential [µV]')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f10(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSE Potenetial DOI [µV]')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        % figure11 RMSE ptherm and recon
        f11(iAmp).f = figure(11+(iAmp-1)*figsprun);
        if iIn<=2
            
            f11(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSE Potential + Noise_{therm} [µV]')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f11(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSE Reconstructed Signal [µV]')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f11(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('RMSE Potential + Noise_{therm} [µV]')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f11(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,RMSEreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('RMSE Reconstructed Signal [µV]')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        
     % figure 12 RMSE doi VR
        f12(iAmp).f = figure(12+(iAmp-1)*figsprun);
        if iIn<=2
            
            f12(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE Potential')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f12(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE VR DOI ')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f12(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE Potential ')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f12(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE Potenetial DOI')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end
        % figure13 nmRMSE ptherm and recon
        f13(iAmp).f = figure(13+(iAmp-1)*figsprun);
        if iIn<=2
            
            f13(iAmp).sf(1) = subplot(2,2,1);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE Potential + Noise_{therm} [µV]')
                set(gca,'xscale','log')
                %legend('show')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
            f13(iAmp).sf(3) = subplot(2,2,3);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            if iIn==2 || iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
            
        else
            
            f13(iAmp).sf(2) = subplot(2,2,2);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                hold on
                for iPOI = 1:length(POIstp)
                    ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(POIstp(iPOI),:),...
                        'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(POIstp(iPOI)).info));
                end
                for iIn2=1:iInmax
                    ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                        'DisplayName',sprintf('%s',names{iIn2}));
                end
                hold off
                
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE Potential + Noise_{therm} ')
                set(gca,'xscale','log')
                legend(ph(1:end),'location','eastoutside')
                title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
                ph = [];
            end
            
            f13(iAmp).sf(4) = subplot(2,2,4);
            hold on
            for iPOI = POIstp
                plot(Totaldps,nmRMSEreconS(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
            end
            hold off
            
            
            if iIn==iInmax
                xlabel('number of dipoles')
                ylabel('n_{max}RMSE Reconstructed Signal')
                set(gca,'xscale','log')
                %legend('show')
                %title(sprintf('vibration amplitude: %3.2e µm',Amp(iAmp)*1e6))
            end
        end    
        
        
    end
    
end
subp = {sf1};
figs = {f1};
fignames = {'time'};
fignprefix = {'RMSREVR-VRDOI', 'RMSREVRphterm-recon','nRMSEVR-VRDOI',...
    'nRMSEVRptherm-recon','RsqVR_VRDOI','RsqVRptherm-recon','nmmRMSEVR_VRDOI','nmmRMSEVRptherm-recon',...
    'RMSEVR_VRDOI','RMSEVRptherm-recon','nmRMSEVR_VRDOI','nmRMSEVRptherm-recon'};
for iAmp = 1:length(Amp)
subp = horzcat(subp,{f2(iAmp).sf,f3(iAmp).sf,f4(iAmp).sf,f5(iAmp).sf,f6(iAmp).sf,f7(iAmp).sf,...
    f8(iAmp).sf,f9(iAmp).sf,f10(iAmp).sf,f11(iAmp).sf,f12(iAmp).sf,f13(iAmp).sf});
figs = horzcat(figs,{f2(iAmp).f,f3(iAmp).f,f4(iAmp).f,f5(iAmp).f,f6(iAmp).f,f7(iAmp).f,...
    f8(iAmp).f,f9(iAmp).f,f10(iAmp).f,f11(iAmp).f,f12(iAmp).f,f13(iAmp).f});
fignames = horzcat(fignames,cellfun(@(x) [x,num2str(Amp(iAmp)*1e6),'µm'],fignprefix,'UniformOutput',false));   
end
%change figure settings 
dl_pl = 0.05;
fig_pos = [-1701,121,1465,748];
for isubp = 1:length(subp)
    set(figs{isubp},'position',fig_pos)
    set(figs{isubp},'position',fig_pos);
    pause(0.01)
    set(findall(figs{isubp},'-property','Fontsize'),'Fontsize',20)
    set(findobj(figs{isubp},'type','axes'),'Fontsize',15)
    set(findobj(figs{isubp},'type','line'),'Linewidth',1.5)
    
    max_x_f1 = subp{isubp}(2).Position(1)+subp{isubp}(2).Position(3);
    lsp = (max_x_f1-dl_pl-subp{isubp}(1).Position(1))/2;
    pos_sp2 = subp{isubp}(1).Position(1)+lsp+dl_pl;
    %adjust width sp1
    subp{isubp}(1).Position(3) = lsp;
    %adjust width sp3
    subp{isubp}(3).Position(3) = lsp;
    
    YLIMS_time(1,:) = get(subp{isubp}(1),'ylim');
    YLIMS_difftime(1,:) = get(subp{isubp}(3),'ylim');
    nrc_sp = floor(length(subp{isubp})/2);
    for isp=2:floor(length(subp{isubp})/2)
    %adjust pos and width sp2
    subp{isubp}(isp).Position(1) = subp{isubp}(isp-1).Position(1)+lsp+dl_pl;
    subp{isubp}(isp).Position(3) = lsp;
    
    %adjust pos and width sp2
    subp{isubp}(isp+nrc_sp).Position(1) = subp{isubp}(isp+nrc_sp-1).Position(1)+lsp+dl_pl;
    subp{isubp}(isp+nrc_sp).Position(3) = lsp;
    
    YLIMS_time(isp,:) = get(subp{isubp}(isp),'ylim');
    YLIMS_difftime(isp,:) = get(subp{isubp}(isp+nrc_sp),'ylim');
    end
    YLIMS_time = [min(YLIMS_time(:,1)),max(YLIMS_time(:,2))];
    YLIMS_difftime = [min(YLIMS_difftime(:,1)),max(YLIMS_difftime(:,2))];
    for isp = 1:nrc_sp
        set(subp{isubp}(isp),'ylim',YLIMS_time,'xtick',Totaldps(1:2:end),...
            'xtickLabel',cellstr(num2str(round(log10(Totaldps(1:2:end)')), '10^%d')))
        set(subp{isubp}(isp+nrc_sp),'ylim',YLIMS_difftime,'xtick',Totaldps(1:2:end),...
            'xtickLabel',cellstr(num2str(round(log10(Totaldps(1:2:end)')), '10^%d')))
        
    end
    
end

%save figures
suffix = datestr(now,'yymmddHHMM');
fignames = cellfun(@(str) ['./figDBstudy/',str,'_',suffix],fignames,'uniformOutput',false);
for ifig = 1:length(figs)
    savefig(figs{ifig},[fignames{ifig},'.fig']);
end