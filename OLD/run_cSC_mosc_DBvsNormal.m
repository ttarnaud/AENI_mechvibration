%run_calcStrengthCurve_mosc_DBvsnormal
if ~isdeployed
clear all
close all
clc
end
Inputnames = {'Input_Mosc_sM_v2_1','Input_Mosc_sM_v2_2','Input_Mosc_sM_v2_3','Input_Mosc_sM_v2_5'};
Colors= [        0    0.4470    0.7410
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
Totaldps = round(logspace(1,5,9));
names = {'Deep single layer','Cortex single layer','Deep 10 layers','Deep Hotspots'};
orders = -floor(log10(Totaldps));
for i=1:length(Totaldps)
    Totaldps(i) = round(Totaldps(i),orders(i));
end
for iIn = 1:length(Inputnames)
    RMSREVR = [];
    RMSREVRDOI = [];
    RMSREVRptherm=[];
    RMSREreconSignalVR = [];
    nRMSEVR = [];
    nRMSEVRDOI = [];
    nRMSEVRptherm=[];
    nRMSEreconSignalVR = [];
    timevals_n = [];
    timevals_DB = [];
    [RMSall,SNRall,Outall,TSTOPall,Totaldps] = calcStrengthCurve_mosc_DBvsnormal(Inputnames{iIn});
    
    for idps = 1:length(Totaldps)
    timevals_n(idps) = TSTOPall(1,idps).TSTOP.endsim-TSTOPall(1,idps).TSTOP.Itgen(1);
    timevals_DB(idps) = TSTOPall(1,idps).TSTOP_DB.endsim-TSTOPall(1,idps).TSTOP_DB.Itgen(1);
    end
    f1 = figure(1);
    diff_tval_nDB = timevals_n-timevals_DB;
    sf1(1) = subplot(2,2,1);
    hold on
    plot(Totaldps,timevals_n,line_types{1},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, normal',names{iIn}));
    plot(Totaldps,timevals_DB,line_types{2},'Color',Colors(iIn,:),'DisplayName',sprintf('%s, database',names{iIn}));
    hold off
    if iIn==length(Inputnames)
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
    if iIn==length(Inputnames)
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
    if iIn==length(Inputnames)
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
    if iIn==length(Inputnames)
        xlabel('number of dipoles')
        ylabel('simulation time (s)')
        set(gca,'xscale','log')
        legend('show','location','eastoutside')
        title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
    end
    
    % calculate data fig 2,3,4&5
    for idps = 1:length(Totaldps)
        RMSREVR(:,idps) = sqrt(nanmean(((Outall(1,idps).Out.VRDOI-Outall(1,idps).Out_DB.VRDOI)./Outall(1,idps).Out.VRDOI).^2,2));
        RMSREVRDOI(:,idps) = sqrt(nanmean(((Outall(1,idps).Out.VRDOI-Outall(1,idps).Out_DB.VRDOI)./Outall(1,idps).Out.VRDOI).^2,2));
        RMSREVRptherm(:,idps) = sqrt(nanmean(((Outall(1,idps).Out.VRptherm-Outall(1,idps).Out_DB.VRptherm)./Outall(1,idps).Out.VRptherm).^2,2));
        nRMSEVR(:,idps) = sqrt(mean(((Outall(1,idps).Out.VRDOI-Outall(1,idps).Out_DB.VRDOI)).^2,2))./mean(Outall(1,idps).Out.VRDOI,2);
        nRMSEVRDOI(:,idps) = sqrt(mean(((Outall(1,idps).Out.VRDOI-Outall(1,idps).Out_DB.VRDOI)).^2,2))./mean(Outall(1,idps).Out.VRDOI,2);
        nRMSEVRptherm(:,idps) = sqrt(mean(((Outall(1,idps).Out.VRptherm-Outall(1,idps).Out_DB.VRptherm)).^2,2))./mean(Outall(1,idps).Out.VRptherm,2);
        for iPOI = 1:length(Outall(1,idps).Out.reconsSignalVR)
        RMSREreconSignalVR(iPOI,idps) = sqrt(nanmean(((Outall(1,idps).Out.reconsSignalVR{iPOI}-Outall(1,idps).Out_DB.reconsSignalVR{iPOI})./Outall(1,idps).Out.reconsSignalVR{iPOI}).^2,2));
        nRMSEreconSignalVR(iPOI,idps) = sqrt(mean(((Outall(1,idps).Out.reconsSignalVR{iPOI}-Outall(1,idps).Out_DB.reconsSignalVR{iPOI})).^2,2))./mean(Outall(1,idps).Out.reconsSignalVR{iPOI},2);
        end
    end
    
    
    % figure 2 using RMSRE
    f2 = figure(2);
    if iIn<=2
        
        sf2(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMSRE Potential')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
        sf2(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMSRE Potential DOI only')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf2(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(RMSREVR,1)
                ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                    'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                    'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off
            xlabel('number of dipoles')
            ylabel('RMSRE Potential')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
            ph = [];
        end
        
        sf2(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('RMSRE Potential DOI only')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
    
    
    % figure 3 using RMSRE
    f3 = figure(3);
    if iIn<=2
        
        sf3(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMSRE Potential + Noise_{therm}')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
        sf3(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMSRE Reconstructed Signal')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf3(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(RMSREVR,1)
            ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off

            xlabel('number of dipoles')
            ylabel('RMSRE Potential + Noise_{therm}')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
            ph = [];
        end
        
        sf3(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('RMSRE Reconstructed Signal')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
    
    % figure 4 using nRMSE
    f4 = figure(4);
    if iIn<=2
        
        sf4(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Potential')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
        sf4(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Potential DOI only')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf4(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(nRMSEVR,1)
                ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                    'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                    'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off
            xlabel('number of dipoles')
            ylabel('nRMSE Potential')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
            ph = [];
        end
        
        sf4(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('nRMSE Potential DOI only')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
    
    
    % figure 3
    f5 = figure(5);
    if iIn<=2
        
        sf5(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Potential + Noise_{therm}')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
        sf5(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Reconstructed Signal')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf5(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(nRMSEVR,1)
            ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off

            xlabel('number of dipoles')
            ylabel('nRMSE Potential + Noise_{therm}')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
            ph = [];
        end
        
        sf5(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:))%,...
                %'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('nRMSE Reconstructed Signal')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
    pause(1)
    
    % calculate data fig 6,7,8&9
    for idps = 1:length(Totaldps)
        RMSREVR(:,idps) = sqrt(nanmean(((Outall(end,idps).Out.VRDOI-Outall(end,idps).Out_DB.VRDOI)./Outall(end,idps).Out.VRDOI).^2,2));
        RMSREVRDOI(:,idps) = sqrt(nanmean(((Outall(end,idps).Out.VRDOI-Outall(end,idps).Out_DB.VRDOI)./Outall(end,idps).Out.VRDOI).^2,2));
        RMSREVRptherm(:,idps) = sqrt(nanmean(((Outall(end,idps).Out.VRptherm-Outall(end,idps).Out_DB.VRptherm)./Outall(end,idps).Out.VRptherm).^2,2));
        nRMSEVR(:,idps) = sqrt(mean(((Outall(end,idps).Out.VRDOI-Outall(end,idps).Out_DB.VRDOI)).^2,2))./mean(Outall(end,idps).Out.VRDOI,2);
        nRMSEVRDOI(:,idps) = sqrt(mean(((Outall(end,idps).Out.VRDOI-Outall(end,idps).Out_DB.VRDOI)).^2,2))./mean(Outall(end,idps).Out.VRDOI,2);
        nRMSEVRptherm(:,idps) = sqrt(mean(((Outall(end,idps).Out.VRptherm-Outall(end,idps).Out_DB.VRptherm)).^2,2))./mean(Outall(end,idps).Out.VRptherm,2);
        for iPOI = 1:length(Outall(end,idps).Out.reconsSignalVR)
        RMSREreconSignalVR(iPOI,idps) = sqrt(nanmean(((Outall(end,idps).Out.reconsSignalVR{iPOI}-Outall(end,idps).Out_DB.reconsSignalVR{iPOI})./Outall(end,idps).Out.reconsSignalVR{iPOI}).^2,2));
        nRMSEreconSignalVR(iPOI,idps) = sqrt(mean(((Outall(end,idps).Out.reconsSignalVR{iPOI}-Outall(end,idps).Out_DB.reconsSignalVR{iPOI})).^2,2))./mean(Outall(end,idps).Out.reconsSignalVR{iPOI},2);
        end
    end
    
    
    % figure 6
    f6 = figure(6);
    if iIn<=2
        
        sf6(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMS Potential (µV)')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
        end
        
        sf6(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMS Potential DOI only (µV)')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf6(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(RMSREVR,1)
                ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                    'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                    'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off
            xlabel('number of dipoles')
            ylabel('RMS Potential (µV)')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
            ph = [];
        end
        
        sf6(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('RMS Potential DOI only (µV)')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
    
    
    % figure 7
    f7 = figure(7);
    if iIn<=2
        
        sf7(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMS Potential + Noise_{therm} (µV)')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
        end
        
        sf7(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('RMS Reconstructed Signal (pV)')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf7(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(RMSREVR,1)
                ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                    'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                    'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off
            xlabel('number of dipoles')
            ylabel('RMS Potential + Noise_{therm} (µV)')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
            ph = [];
        end
        
        sf7(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(RMSREVR,1)
            plot(Totaldps,RMSREreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('RMS Reconstructed Signal (pV)')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
    
    % figure 8 nRMSE
    f8 = figure(8);
    if iIn<=2
        
        sf8(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Potential')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
        end
        
        sf8(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVRDOI(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Potential DOI only')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf8(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(RMSREVR,1)
                ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                    'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                    'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off
            xlabel('number of dipoles')
            ylabel('nRMSE Potential')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
            ph = [];
        end
        
        sf8(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('nRMSE Potential DOI only')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
    
    
    % figure 9
    f9 = figure(9);
    if iIn<=2
        
        sf9(1) = subplot(2,2,1);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
                
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Potential + Noise_{therm}')
            set(gca,'xscale','log')
            %legend('show')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
        end
        
        sf9(3) = subplot(2,2,3);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        if iIn==2
            xlabel('number of dipoles')
            ylabel('nRMSE Reconstructed Signal')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
        
    else
        
        sf9(2) = subplot(2,2,2);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEVRptherm(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            hold on
            for iPOI = 1:size(nRMSEVR,1)
                ph(iPOI) = plot(Totaldps,nan*Totaldps,'Color',Colors(iPOI,:),...
                    'DisplayName',sprintf('POI: %s',RMSall(1,1).RMS(iPOI).info));
            end
            for iIn2=1:length(Inputnames)
                ph(end+1) = plot(Totaldps,nan*Totaldps,['k',line_types{iIn2}],...
                    'DisplayName',sprintf('%s',names{iIn2}));
            end
            hold off
            xlabel('number of dipoles')
            ylabel('nRMSE Potential + Noise_{therm}')
            set(gca,'xscale','log')
            legend(ph(1:end),'location','southeastoutside')
            title(sprintf('vibration amplitude: %3.2e µm',Amp(end)*1e6))
            ph = [];
        end
        
        sf9(4) = subplot(2,2,4);
        hold on
        for iPOI = 1:size(nRMSEVR,1)
            plot(Totaldps,nRMSEreconSignalVR(iPOI,:),line_types{iIn},'Color',Colors(iPOI,:),...
                'DisplayName',sprintf('%s, POI: %s',names{iIn},RMSall(1,1).RMS(iPOI).info));
        end
        hold off
        
        
        if iIn==4
            xlabel('number of dipoles')
            ylabel('nRMSE Reconstructed Signal')
            set(gca,'xscale','log')
            %legend('show')
            %title(sprintf('vibration amplitude: %3.2e µm',Amp(1)*1e6))
        end
    end
end

subp = {sf1,sf2,sf3,sf4,sf5,sf6,sf7,sf8,sf9};
for isubp = 1:length(subp)
    max_x_f1 = subp{isubp}(2).Position(1)+subp{isubp}(2).Position(3);
    lsp = (max_x_f1-0.1-subp{isubp}(1).Position(1))/2;
    pos_sp2 = subp{isubp}(1).Position(1)+lsp+0.1;
    
    %adjust width sp1
    subp{isubp}(1).Position(3) = lsp;
    %adjust pos and width sp2
    subp{isubp}(2).Position(1) = pos_sp2;
    subp{isubp}(2).Position(3) = lsp;
    %adjust width sp3
    subp{isubp}(3).Position(3) = lsp;
    %adjust pos and width sp4
    subp{isubp}(4).Position(1) = pos_sp2;
    subp{isubp}(4).Position(3) = lsp;

end

figs = {f1,f2,f3,f4,f5,f6,f7,f8,f9};
fignames = {'time', 'RMSREVR-VRDOI_minA', 'RMSREVRphterm-recon_minA','nRMSEVR-VRDOI_minA',...
    'nRMSEVRphterm-recon_minA', 'RMSREVR-VRDOI_maxA', 'RMSREVRphterm-recon_maxA', 'nRMSEVR-VRDOI_maxA', 'nRSMEVRphterm-recon_maxA'};
suffix = datestr(now,'yymmddHHMM');
fignames = cellfun(@(str) [str,'_',suffix],fignames,'uniformOutput',false); 
for ifig = 1:length(figs)
    savefig(figs{ifig},[fignames{ifig},'.fig']);
end