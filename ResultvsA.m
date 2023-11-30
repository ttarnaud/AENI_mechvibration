% With this scrip we test the linear dependece wrt to the vibration
close all
clear all
clc
addpath(genpath('./Functions'));
addpath(genpath('./Inputs'));
Colors= [        0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
% solution types
% SolutionTypes = {'3SphereS8.2R25','3SphereS8.2R~f','4SphereS8.7R25',...
%     '4SphereS8.7R~f','Mouse4Sphere~fair0','Mouse4Sphere~fair1'};
SolutionTypes2p = {'3SphereS8.2R25_X','3SphereS8.2R25_Y','3SphereS8.2R25_Z','3SphereS8.2R25_XYZ',};
% SolutionTypes2p = {'3SphereS8.2R25_X','3SphereS8.2R25_Y','3SphereS8.2R25_Z','3SphereS8.2R25_XYZ','4SphereS8.7R25_X',...
%     '4SphereS8.7R25_Y','4SphereS8.7R25_Z','4SphereS8.7R25_XYZ','Mouse4Sphere~fair1_X',...
%     'Mouse4Sphere~fair1_Y','Mouse4Sphere~fair1_Z','Mouse4Sphere~fair1_XYZ'};
Allinfo = struct();      %structure that contains info of first two subscripts and is saved
SCATTERPLOT = 1;         %create scatterplot in second subscript
fignr2 = 0;              %
PLOT = 1;                % PLOT_flag
calclin_flag = 1;        %calculate linear relationship and plot
Tend = 10e-3;            %end of each simulation: has to be large enough to contain enough datapoints in reconstruction
DEBUG_flag = 1;
if DEBUG_flag
    warning('DEBUG_FLAG: ON!!!')
end

% Declaration POIs
POIs = [];
dtheta = 4;
for itheta = 0:pi/dtheta:pi
    dphi = pi/4;
    if itheta==0 || itheta==pi
        dphi = 2*pi;
    end
    for iphi=[0:dphi:2*pi-dphi]
        POIs = [POIs;cos(iphi)*sin(itheta),sin(iphi)*sin(itheta),cos(itheta)];
    end
end

%% this subscript: we check for two dipole positions ([0,0,],[0,0,Rmax]) if output is linear wrt vibration amplitude
% for a fixed vibration direction. Vibration direction is given in the
% solutionTypes2p name. we loop over Aus2p = logspace(0,Aus_max,n) values
% the dipole is radially oriented
if DEBUG_flag
    isol_end = 1;
else
    isol_end = length(SolutionTypes2p);
end

for isol=1:isol_end
    % first check for two positions at which Aus not linear anymore
    
    if contains(SolutionTypes2p{isol},'Mouse')
        Rmax = 4.1;
        if DEBUG_flag
            dR = 1;
        else
            dR=0.1;
        end
        Aus_max = 5;
    else
        Rmax = 65;
        if DEBUG_flag
            dR = 15;
        else 
            dR = 1;
        end
        Aus_max = 6;
    end
    dp_pos2p = [0,0,0;0,0,Rmax]; dp_pos2p = dp_pos2p/1000;
    orien_dp2p = [0,0,1;0,0,1];
    if DEBUG_flag
        Aus2p = logspace(0,Aus_max,3)*1e-9;
    else
        Aus2p = logspace(0,Aus_max,20)*1e-9;
    end
    %Aus2p=linspace(1e-3,1e-2,20);
    
    % Declare other settings
    idx_Dir_flag = strfind(SolutionTypes2p{isol},'_');
    Dir_flag = SolutionTypes2p{isol}(idx_Dir_flag+1:end);
    switch Dir_flag
        case 'X'; Dirus = [1,0,0];
        case 'Y'; Dirus = [0,1,0];
        case 'Z'; Dirus = [0,0,1];
        case 'XYZ'; Dirus = [1,1,1]./(sqrt(3));
        otherwise; error('wrong input');
    end
    SolutionType = SolutionTypes2p{isol}(1:idx_Dir_flag-1);
    Ifun = @(t) 1+0*t;
    fus = 1e6;
    
    
    
    disp(datetime('now','Format','HH:mm:ss'))
    Pflag = 0;
    diffrecon2p = zeros(size(dp_pos2p,1),length(Aus2p)-1,size(POIs,1)+2); %difference reconstructed signal of two positions
    reldiffrecon2p = zeros(size(dp_pos2p,1),length(Aus2p)-1,size(POIs,1)+2); %relative difference
    diffAmp2p = zeros(size(dp_pos2p,1),length(Aus2p)-1,size(POIs,1));  %difference of amplitude at 1MHz
    reldiffAmp2p = zeros(size(dp_pos2p,1),length(Aus2p)-1,size(POIs,1)); %relative difference
    dA02p = zeros(size(dp_pos2p,1),length(Aus2p)-1,size(POIs,1));   %real outcome at certain intensity (not linear scaled)
    dr02p = zeros(size(dp_pos2p,1),length(Aus2p)-1,size(POIs,1)+2);
    for idp=1:size(dp_pos2p,1)
        % first calculate for lowest Vibration amplitude
        [~,~,~,Out1] = investBiologicalNoise('posDp',dp_pos2p(idp,:),'POIs',POIs,'OrienDipole',orien_dp2p(idp,:),...
            'Plot',0,'Display',0,'totaldps',1,'dps_run',1,'showSphere','','Solutiontype',SolutionType,...
            'Ifun',Ifun,'fus',fus,'Aus',Aus2p(1),'Dirus',Dirus,'Tend',Tend);
        AmpVR1 = (max(Out1.VR,[],2)-min(Out1.VR,[],2))/2;   %because Ifun is constant amplitude can be defined like this
        
        dA02p(idp,1,:) = AmpVR1;
        dr02p(idp,1,:) = cellfun(@mean,Out1.reconsSignalVRDOI);
        
        parfor iaus=2:length(Aus2p)
            %calculate for other amplitudes
            [~,~,~,Out] = investBiologicalNoise('posDp',dp_pos2p(idp,:),'POIs',POIs,'OrienDipole',orien_dp2p(idp,:),...
                'Plot',0,'Display',0,'totaldps',1,'dps_run',1,'showSphere','','Solutiontype',SolutionType,...
                'Ifun',Ifun,'fus',fus,'Aus',Aus2p(iaus),'Dirus',Dirus,'Tend',Tend);
            AmpVR = (max(Out.VR,[],2)-min(Out.VR,[],2))/2;
            diffAmp2p(idp,iaus-1,:) = (Aus2p(iaus)/Aus2p(1))*AmpVR1-AmpVR;  %linear scale result when lowest amplitude selected and take differenc3e
            reldiffAmp2p(idp,iaus-1,:) = squeeze(diffAmp2p(idp,iaus-1,:))./AmpVR;
            dA02p(idp,iaus,:) = AmpVR;
            
            diffreconarray = zeros(1,size(POIs,1)+2);
            reldiffreconarray = zeros(1,size(POIs,1)+2);
            dr0array = zeros(1,size(POIs,1)+2);
            for ipoi=1:size(POIs,1)+2
                diffreconarray(ipoi) = mean((Aus2p(iaus)/Aus2p(1))*(Out1.reconsSignalVRDOI{ipoi})-Out.reconsSignalVRDOI{ipoi});
                reldiffreconarray(ipoi) = diffreconarray(ipoi)/mean(Out.reconsSignalVRDOI{ipoi});
                dr0array(ipoi) = mean(Out.reconsSignalVRDOI{ipoi});
                
            end
            diffrecon2p(idp,iaus-1,:) = diffreconarray;
            reldiffrecon2p(idp,iaus-1,:) = reldiffreconarray;
            dr02p(idp,iaus,:) = dr0array;
        end
        Pval = round(idp/size(dp_pos2p,1)*100,0);
        if Pval >= Pflag
            disp(datetime('now','Format','HH:mm:ss'))
            disp(['Progress: ',num2str(Pval),'%'])
            Pflag = Pflag+5;
        end
        
    end
    % stor info to be saved
    Allinfo.twop(isol).diffAmp2p = diffAmp2p;
    Allinfo.twop(isol).reldiffAmp2p = reldiffAmp2p;
    Allinfo.twop(isol).diffrecon2p = diffrecon2p;
    Allinfo.twop(isol).reldiffrecon2p = reldiffrecon2p;
    Allinfo.twop(isol).dp_pos2p = dp_pos2p;
    Allinfo.twop(isol).orien_dp2p = orien_dp2p;
    Allinfo.twop(isol).POIs=POIs;
    Allinfo.twop(isol).SolutionType = SolutionTypes2p{isol};
    Allinfo.twop(isol).dA02p = dA02p;
    Allinfo.twop(isol).dr02p = dr02p;
    
    if PLOT
        for idp=1:size(dp_pos2p,1)
            fignr1 = fignr2+1;  %figure 1 is for amplitude
            fignr2 = fignr1+1;  %figure 2 is for reconstructed signal
            
            for iPOI = 1:size(POIs,1)
                % plot 4 lines in each subplot except last one (if not
                % devidable by 4)
                if iPOI>floor(size(POIs,1)/4)*4
                    spnr = floor(size(POIs,1)/4);
                else
                    spnr = ceil(iPOI/4);
                end
                figure(fignr1)
                subplot(floor(size(POIs,1)/4)/2,2,spnr)
                hold on
                plot(Aus2p,squeeze(dA02p(idp,:,iPOI)),'Color',Colors(mod(iPOI-1,4)+1,:),'DisplayName',...
                    sprintf('POI: \\phi = %3.0f°, \\theta = %3.0f°',rad2deg(atan2(POIs(iPOI,2),POIs(iPOI,1))),rad2deg(acos(POIs(iPOI,3)))));
                pause(0.1)
                legend('show','location','eastoutside')
                xlabel('Amplitude Vibration')
                ylabel('Signal strength [µV]')
                              
                
                if calclin_flag
                    if length(Aus2p)>15
                        idxend = 10;
                    else
                        idxend = 0;
                    end
                    b=Aus2p(1:end-idxend)'\squeeze(dA02p(idp,1:end-idxend,iPOI))';
                    yfun =@(x) b.*x;
                    p1=plot(Aus2p,yfun(Aus2p),'--','Color',Colors(mod(iPOI-1,4)+1,:));
                    pause(0.1)
                    set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');                    
                end
                pause(0.1)
                hold off
                set(gca,'xscale','log','yscale','log') 
                
                %plot results of reconstructed signal
                figure(fignr2)
                subplot(floor(size(POIs,1)/4)/2,2,spnr)
                hold on
                plot(Aus2p,squeeze(abs(dr02p(idp,:,iPOI))),'Color',Colors(mod(iPOI-1,4)+1,:),'DisplayName',...
                    sprintf('POI: \\phi = %3.0f°, \\theta = %3.0f°',rad2deg(atan2(POIs(iPOI,2),POIs(iPOI,1))),rad2deg(acos(POIs(iPOI,3)))));
                pause(0.1)
                legend('show','location','eastoutside')
                xlabel('Amplitude Vibration')
                ylabel('mean reconSignal')
                
                if calclin_flag
                    if length(Aus2p)>15
                        idxend = 10;
                    else
                        idxend = 0;
                    end
                    b=Aus2p(1:end-idxend)'\squeeze(abs(dr02p(idp,1:end-idxend,iPOI)))';
                    yfun =@(x) b.*x;
                    p1=plot(Aus2p,yfun(Aus2p),'--','Color',Colors(mod(iPOI-1,4)+1,:));
                    pause(0.1)
                    set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');                    
                end
                pause(0.1)
                hold off
                set(gca,'xscale','log','yscale','log')
            end
            %set(gcf,'position',[-1902,50,1889,941])
            mtit(sprintf('recon, Model: %s, dp_{pos} = [%3.1f,%3.1f,%3.1f]mm',SolutionTypes2p{isol},dp_pos2p(idp,:)*1000),'xoff',0,'yoff',0.05)
            figure(fignr1)
            %set(gcf,'position',[-1902,50,1889,941])
            mtit(sprintf('A0, Model: %s, dp_{pos} = [%3.1f,%3.1f,%3.1f]mm',SolutionTypes2p{isol},dp_pos2p(idp,:)*1000),'xoff',0,'yoff',0.05)
            pause(0.1)
        end
        
        
    end
    
end




%% check for high amount of positions if linear wrt vibration amplitude
if DEBUG_flag
    SolutionTypes = {'3SphereS8.2R25'};
else
    SolutionTypes = {'3SphereS8.2R25','4SphereS8.7R25','Mouse4Sphere~fair1'};
end


for isol=1:length(SolutionTypes)
    if contains(SolutionTypes{isol},'Mouse')
        Rmax = 4.1;
        if DEBUG_flag
            dR = 1;
        else
            dR = 0.1;            
        end
        Aus_max = 5;
    else
        Rmax = 65;
        if DEBUG_flag
            dR = 15;
        else
            dR=1;
        end
        Aus_max = 6;
    end
    
    % Declare dipole positions
    dp_pos = [];
    orien_dp = [];
    if DEBUG_flag
        dtheta = 4;
        Aus = [100e-9,1e-9,10e-6];
    else
        dtheta = 20;
        Aus = [100e-9,1e-9,10e-9,10e-6,100e-6];
    end
    for itheta =0:pi/dtheta:pi/2
        for ir=0:dR:Rmax
            dp_pos = [dp_pos;ir*cos(0)*sin(itheta), ir*sin(0)*sin(itheta),ir*cos(itheta)];
            orien_dp = [orien_dp;cos(0)*sin(itheta), cos(0)*sin(itheta),cos(itheta)];
        end
    end
    dp_pos = dp_pos/1000; % from mm to m
    
    
    % Declare other settingsTes
    SolutionType = SolutionTypes{isol};
    Ifun = @(t) 1+0*t;
    fus = 1e6;
    Dirus = [0,0,1];
    

    disp(datetime('now','Format','HH:mm:ss'))
    Pflag = 0;
    diffrecon = zeros(size(dp_pos,1),length(Aus)-1,size(POIs,1)+2);
    reldiffrecon = zeros(size(dp_pos,1),length(Aus)-1,size(POIs,1)+2);
    diffAmp = zeros(size(dp_pos,1),length(Aus)-1,size(POIs,1));
    reldiffAmp = zeros(size(dp_pos,1),length(Aus)-1,size(POIs,1));
    dA0 = zeros(size(dp_pos,1),length(Aus)-1,size(POIs,1));
    dr0 = zeros(size(dp_pos,1),length(Aus)-1,size(POIs,1)+2);
    for idp=1:size(dp_pos,1)
        %run for first Aus
        [~,~,~,Out1] = investBiologicalNoise('posDp',dp_pos(idp,:),'POIs',POIs,'OrienDipole',orien_dp(idp,:),...
            'Plot',0,'Display',0,'totaldps',1,'dps_run',1,'showSphere','','Solutiontype',SolutionType,...
            'Ifun',Ifun,'fus',fus,'Aus',Aus(1),'Dirus',Dirus);
        AmpVR1 = (max(Out1.VR,[],2)-min(Out1.VR,[],2))/2;
        % turn of warning
        warning('off','dps_run:notidealnumber')
        parfor iaus=2:length(Aus)
            %run other aus and compare linear scaled differences
            [~,~,~,Out] = investBiologicalNoise('posDp',dp_pos(idp,:),'POIs',POIs,'OrienDipole',orien_dp(idp,:),...
                'Plot',0,'Display',0,'totaldps',1,'dps_run',1,'showSphere','','Solutiontype',SolutionType,...
                'Ifun',Ifun,'fus',fus,'Aus',Aus(iaus),'Dirus',Dirus);
            AmpVR = (max(Out.VR,[],2)-min(Out.VR,[],2))/2;
            diffAmp(idp,iaus-1,:) = (Aus(iaus)/Aus(1))*AmpVR1-AmpVR;
            reldiffAmp(idp,iaus-1,:) = squeeze(diffAmp(idp,iaus-1,:))./AmpVR;
            dA0(idp,iaus-1,:) = AmpVR;
            
            diffreconarray = zeros(1,size(POIs,1)+2);
            reldiffreconarray = zeros(1,size(POIs,1)+2);
            dr0array = zeros(1,size(POIs,1)+2);
            for ipoi=1:size(POIs,1)+2 %reconstructed cell contains also signal of averaged and averaged over same Order of maginitude POIs
                diffreconarray(ipoi) = mean((Aus(iaus)/Aus(1))*(Out1.reconsSignalVRDOI{ipoi})-Out.reconsSignalVRDOI{ipoi});
                reldiffreconarray(ipoi) = diffreconarray(ipoi)/mean(Out.reconsSignalVRDOI{ipoi});
                dr0array(ipoi) = mean(Out.reconsSignalVRDOI{ipoi});
                
            end
            diffrecon(idp,iaus-1,:) = diffreconarray;
            reldiffrecon(idp,iaus-1,:) = reldiffreconarray;
            dr0(idp,iaus-1,:) = dr0array;
        end
        Pval = round(idp/size(dp_pos,1)*100,0);
        if Pval >= Pflag
            disp(datetime('now','Format','HH:mm:ss'))
            disp(['Progress: ',num2str(Pval),'%'])
            Pflag = Pflag+5;
        end
        
    end
    Allinfo.multip(isol).diffAmp = diffAmp;
    Allinfo.multip(isol).reldiffAmp = reldiffAmp;
    Allinfo.multip(isol).diffrecon = diffrecon;
    Allinfo.multip(isol).reldiffrecon = reldiffrecon;
    Allinfo.multip(isol).dp_pos = dp_pos;
    Allinfo.multip(isol).orien_dp = orien_dp;
    Allinfo.multip(isol).POIs=POIs;
    Allinfo.multip(isol).SolutionType = SolutionTypes{isol};
    Allinfo.multip(isol).dA0 = dA0;
    Allinfo.multip(isol).dr0 = dr0;
    
    %% plot results
    if SCATTERPLOT
        for ifig=1:size(POIs,1)
            fignr = (isol-1)*size(POIs,1)+ifig+fignr2
            figure(fignr)
            for iplot=1:length(Aus)-1
                subplot(3,2,(iplot-1)*2+1)
                scatter3(dp_pos(:,1),dp_pos(:,2),dp_pos(:,3),[],reldiffrecon(:,iplot,ifig))
                hold on
                scatter3(0.082*POIs(ifig,1),0.082*POIs(ifig,2),0.082*POIs(ifig,3),300,'rX')
                hold off
                colormap jet
                colorbar
                
                subplot(3,2,iplot*2)
                scatter3(dp_pos(:,1),dp_pos(:,2),dp_pos(:,3),[],reldiffAmp(:,iplot,ifig))
                hold on
                scatter3(0.082*POIs(ifig,1),0.082*POIs(ifig,2),0.082*POIs(ifig,3),300,'rX')
                hold off
                colormap jet
                colorbar
                
            end
            pause(0.1)
        end
    end
end
save(['OutputLinearitystudy_',datestr(now,'yyyyMMdd-hhmm'),'.mat'],'Allinfo');
%% test if linear in direction results of 2point study
if isol_end == 4
d3cvs1_dA0 = 1/sqrt(3).*(Allinfo.twop(1).dA02p+Allinfo.twop(2).dA02p+Allinfo.twop(3).dA02p)-Allinfo.twop(4).dA02p;
d3cvs1000_dA0 = squeeze(d3cvs1_dA0(1,:,:));
dA02p000_dir111 = squeeze(Allinfo.twop(4).dA02p(1,:,:));
dA02p000_dir111_comb = squeeze(1/sqrt(3).*(Allinfo.twop(1).dA02p(1,:,:)+Allinfo.twop(2).dA02p(1,:,:)+Allinfo.twop(3).dA02p(1,:,:)));
dA02p000_dir100 = squeeze(1/sqrt(3).*(Allinfo.twop(1).dA02p(1,:,:)));
dA02p000_dir010 = squeeze(1/sqrt(3).*(Allinfo.twop(2).dA02p(1,:,:)));
dA02p000_dir001 = squeeze(1/sqrt(3).*(Allinfo.twop(3).dA02p(1,:,:)));
rd3cvs1000 = d3cvs1000_dA0./dA02p000_dir111;
rd3cvs1_dA0 = d3cvs1_dA0./Allinfo.twop(4).dA02p;
perc = 5/100;
fprintf('\ndifferences not equal to zero: %i\n',sum(d3cvs1_dA0(:)~=0))
fprintf('\ndifferences not smaller than %3.2f%%: %i\n',perc*100,sum(rd3cvs1_dA0(:)>perc))
end

%% numerical investigation if signal strength linear dependent of vibration amplitude: this for multiple dipole positions
%results given in table
info_invest = Allinfo.multip(1);
testpercent = 1;
abs_diff = 1e-9;
for i=1:length(Aus)-1
    vals_rdA = info_invest.reldiffAmp(:,i,:);
    vals_rdA_noinf = vals_rdA; vals_rdA_noinf(isinf(vals_rdA_noinf)) = 0; 
    pos_rdA_h(i) = sum(vals_rdA(:)>testpercent/100);
    pos_rdA_noinf_h(i) = sum(vals_rdA_noinf(:)>testpercent/100);
    fprintf('\nAus = %2.2e, number of positions with reldiffAmp > %3.2f%% = %i \n',Aus(i),testpercent,pos_rdA_h(i))
    fprintf('\nAus = %2.2e, number of positions with reldiffAmp_noinf > %3.2f%% = %i \n',Aus(i),testpercent,pos_rdA_noinf_h(i))
    
    vals_dA = info_invest.diffAmp(:,i,:);
    pos_dA_h(i) = sum(vals_dA(:)<=abs_diff);
    %fprintf('\nAus = %2.2e, number of positions with diffAmp > %3.2f%% = %i \n',Aus(i),testpercent,pos_dA_h(i))
    vals_rdA_noinf_0diff = vals_rdA_noinf; vals_rdA_noinf_0diff(vals_dA<=abs_diff) = 0;
    pos_rdA_noinf_0diff_h(i) = sum(vals_rdA_noinf_0diff(:)>testpercent/100); 
    
    vals_rdr = info_invest.reldiffrecon(:,i,:);
    vals_rdr_noinf = vals_rdr; vals_rdr_noinf(isinf(vals_rdr_noinf)) = 0;
    pos_rdr_h(i) = sum(vals_rdr(:)>testpercent/100);
    pos_rdr_noinf_h(i) = sum(vals_rdr_noinf(:)>testpercent/100);
    fprintf('\nAus = %2.2e, number of positions with reldiffrecon > %3.2f%% = %i \n',Aus(i),testpercent,pos_rdr_h(i))
    fprintf('\nAus = %2.2e, number of positions with reldiffrecon_noinf > %3.2f%% = %i \n',Aus(i),testpercent,pos_rdr_noinf_h(i))
    
    vals_dr = info_invest.diffrecon(:,i,:);
    pos_dr_h(i) = sum(vals_dr(:)<abs_diff);
    %fprintf('\nAus = %2.2e, number of positions with diffrecon > %3.2f%% = %i \n',Aus(i),testpercent,pos_dr_h(i))
    vals_rdr_noinf_0diff = vals_rdr_noinf; vals_rdr_noinf_0diff(vals_dr<=abs_diff) = 0;
    pos_rdr_noinf_0diff_h(i) = sum(vals_rdr_noinf_0diff(:)>testpercent/100); 
    
end

T=table(pos_rdA_h',pos_rdA_noinf_h',pos_dA_h',pos_rdA_noinf_0diff_h',pos_rdr_h',pos_rdr_noinf_h',pos_dr_h',pos_rdr_noinf_0diff_h',...
    'Variablenames',{'rdA','rdA_einf','dA','rdA_noinf_0diff','rdr','rdr_woinf','dr','rdr_noinf_0diff'},'Rownames',arrayfun(@num2str,Aus(2:end),'UniformOutput',false))


%% dipole in center and and near surface. create plots of strength visible at sphere boundary like in FlucDistributionCorticalColumn

% POIs
if exist('fignr2','var')
Fignr = fignr2+1;
else
    Fignr = 20;
end
if DEBUG_flag
    Theta = [0:pi/20:pi];
    dPhi = pi/20;
    isol_end = 4;
else
    Theta = [0:pi/100:pi];
    dPhi = pi/50;
    isol_end = length(SolutionTypes2p)
end
[phi,theta]=meshgrid(0:dPhi:2*pi,Theta);
phisc = phi(:); thetasc = theta(:);       
 POIs = 0.082*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];
Pflag = 0;
disp(datetime('now','Format','HH:mm:ss'))
for isol=1:isol_end
    %% first check for two positions at which Aus not linear anymore
    
    if contains(SolutionTypes2p{isol},'Mouse')
        Rmax = 4.1;
        if DEBUG_flag
            dR=1;
        else
            dR=0.1;
        end
        Aus_max = 5;
    else
        Rmax = 65;
        if DEBUG_flag
            dR=15;
        else
            dR=1;
        end
        
        Aus_max = 6;
    end
    dp_pos2p = [0,0,0;0,0,Rmax]; dp_pos2p = dp_pos2p/1000;
    orien_dp2p = [0,0,1;0,0,1];
    Aus = 1e-6;
    
    % Declare other settings
    idx_Dir_flag = strfind(SolutionTypes2p{isol},'_');
    Dir_flag = SolutionTypes2p{isol}(idx_Dir_flag+1:end);
    switch Dir_flag
        case 'X'; Dirus = [1,0,0];
        case 'Y'; Dirus = [0,1,0];
        case 'Z'; Dirus = [0,0,1];
        case 'XYZ'; Dirus = [1,1,1]./(sqrt(3));
        otherwise; error('wrong input');
    end
    SolutionType = SolutionTypes2p{isol}(1:idx_Dir_flag-1);
    Ifun = @(t) 1+0*t;
    fus = 1e6;
    
    
    
    
    for idp=1:size(dp_pos2p,1)
        [~,~,~,Out1] = investBiologicalNoise('posDp',dp_pos2p(idp,:),'POIs',POIs,'OrienDipole',orien_dp2p(idp,:),...
            'Plot',0,'Display',0,'totaldps',1,'dps_run',1,'showSphere','','Solutiontype',SolutionType,...
            'Ifun',Ifun,'fus',fus,'Aus',Aus,'Dirus',Dirus,'Tend',2*fus^-1,'noRecon',1);
        AmpVR1 = (max(Out1.VR,[],2)-min(Out1.VR,[],2))/2*10^6;
        VRflucSurf = reshape(AmpVR1,size(phi));
        dA0mPOI(idp,isol,:) = AmpVR1;
    end
    Pval = round(isol/length(SolutionTypes2p)*100,0);
    if Pval >= Pflag
        disp(datetime('now','Format','HH:mm:ss'))
        disp(['Progress: ',num2str(Pval),'%'])
        Pflag = Pflag+5;
    end
end
%% creating surface plots of strength measurable on the outer sphere (similar to plots created in Flucdistribution Coritcial column
Diruss=[1,0,0;0,1,0;0,0,1;1,1,1;1,1,1];
Diruss = Diruss./vecnorm(Diruss,2,2);
for idp=1:size(dp_pos2p,1) 
    figure(Fignr+idp)
    for isol=1:length(SolutionTypes2p)
        VRflucSurf = reshape(dA0mPOI(idp,isol,:),size(phi));
        sp{idp,isol} = subplot(2,3,isol);
        surf(phi,theta,VRflucSurf,'EdgeColor','none')
        xlabel('Azimuth')
        ylabel('Co-elevation')
        zlabel('Fluctuation amplitude [pV]')
        set(gca,'XTick',0:pi/2:2*pi) 
        set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
        set(gca,'YTick',0:pi/4:pi) 
        set(gca,'YTickLabel',{'0','pi/4','pi/2','3*pi/4','pi'})
        colormap('jet')
        c{idp,isol} = colorbar;
        c{idp,isol}.Label.String = 'Fluctuation Amplitude [pV]';
        xlim([0,2*pi]);
        ylim([0,pi]);
        title(['Vibration [',strjoin(arrayfun(@(x) num2str(x),round(Diruss(isol,:),2),'UniformOutput',false),','),']'])
        view(0,90)
        pause(0.1)
        if isol==length(SolutionTypes2p)
            VRflucSurf = reshape(Diruss(4,:)*squeeze(dA0mPOI(idp,1:3,:)),size(phi));
            
            sp{idp,5} = subplot(2,3,5);
            surf(phi,theta,VRflucSurf,'EdgeColor','none')
            xlabel('Azimuth')
            ylabel('Co-elevation')
            zlabel('Fluctuation amplitude [pV]')
            set(gca,'XTick',0:pi/2:2*pi)
            set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
            set(gca,'YTick',0:pi/4:pi)
            set(gca,'YTickLabel',{'0','pi/4','pi/2','3*pi/4','pi'})
            colormap('jet')
            c{idp,5} = colorbar;
            c{idp,5}.Label.String = 'Fluctuation Amplitude [pV]';
            xlim([0,2*pi]);
            ylim([0,pi]);
            title(['combined: Vibration [',strjoin(arrayfun(@(x) num2str(x),round(Diruss(5,:),2),'UniformOutput',false),','),']'])
            view(0,90)
            %set(gcf,'position',[-1902,50,1889,941])
        end
    end

end
pause(0.1)
for idp=1:size(dp_pos2p,1)
for isol = 1:size(sp,2)
    LIMS(isol,:) = get(c{idp,isol},'Limits');
end
LIMS = [min(LIMS(:,1)),max(LIMS(:,2))];
for isol = 1:size(sp,2)
    set(sp{idp,isol},'CLIM',LIMS);
    set(sp{idp,isol}.Colorbar,'Limits',LIMS);
    if isol~=3 && isol~=size(sp,2)
    set(sp{idp,isol}.Colorbar,'Visible', 'off')
    end
end
possp5 = get(sp{idp,5},'position');
possp4 = get(sp{idp,4},'position');
dpsppos = (possp5-possp4)/2;
set(sp{idp,4},'position',possp4+dpsppos);
pause(1)
set(sp{idp,5},'position',possp5+dpsppos);
pause(0.1)
end

Fignr = get(gcf,'number');
Fignr = Fignr+1;
%% test if taking 10 from set of 1000 has same result as new distri
%her we do a small test is a small subset of uniform distribution is same
%distributed if small subset uniformly distributed themselfs. this could be
%usefull to upgrade a database. If not the case, a whole new set needs to
%be determined
N = 100;
indices = [0:N-1]+0.5;
theta = acos(1-2*indices'/N);
phi = pi*(1+5^(1/2))*indices';
Vibdir = [cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta)];
scatter3(Vibdir(:,1),Vibdir(:,2),Vibdir(:,3))
% allcomb = combnk(1:N,3);
% for i=1:size(allcomb,1)
%     selcomb = allcomb(i,:);
%     surf1(i) = 1/2*norm(cross(Vibdir(selcomb(2),:)-Vibdir(selcomb(1),:),Vibdir(selcomb(3),:)-Vibdir(selcomb(1),:)));
% end
% minsurf1 = sqrt(N)*min(surf1)

%calculate the minimal distance between points generated
all2comb = combnk(1:N,2);
for i=1:size(all2comb,1)
    selcomb = all2comb(i,:);
    dist(i) = norm(Vibdir(selcomb(2),:)-Vibdir(selcomb(1),:));
end
mindist1 = min(dist)

selVibdir = Vibdir(1:10:end,:);
allcomb = combnk(1:10,3);
for i=1:size(allcomb,1)
    selcomb = allcomb(i,:);
    surf2(i) = 1/2*norm(cross(selVibdir(selcomb(2),:)-selVibdir(selcomb(1),:),selVibdir(selcomb(3),:)-selVibdir(selcomb(1),:)));
end
minsurf2 = sqrt(10)*min(surf2)
figure
scatter3(selVibdir(:,1),selVibdir(:,2),selVibdir(:,3))
hold on
N = 10;
indices = [0:N-1]+0.5;
theta = acos(1-2*indices'/N);
phi = pi*(1+5^(1/2))*indices';
Vibdir = [cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta)];
scatter3(Vibdir(:,1),Vibdir(:,2),Vibdir(:,3))
hold off
axis square
allcomb = combnk(1:N,3);
for i=1:size(allcomb,1)
    selcomb = allcomb(i,:);
    surf3(i) = 1/2*norm(cross(Vibdir(selcomb(2),:)-Vibdir(selcomb(1),:),Vibdir(selcomb(3),:)-Vibdir(selcomb(1),:)));
end
minsurf3 = sqrt(10)*min(surf3) %sqrt(10)?

%% Create 5 lookup tables consisting of N vibration directions and dipole at two positions
%signal measured at several POIs. few sets are uniformly idstributed oe is
%random Tho vibration directions of interest are createdd with the for loop
% compare last with interpolation results of lookup tables
Vibdir = {};
theta = {};
phi = {};

N = 100;
indices = [0:N-1]+0.5;
theta{1} = acos(1-2*indices'/N);
phi{1} = pi*(1+5^(1/2))*indices';
Vibdir{1} = [cos(phi{1}).*sin(theta{1}),sin(phi{1}).*sin(theta{1}),cos(theta{1})];

N = 1000;
indices = [0:N-1]+0.5;
theta{2} = acos(1-2*indices'/N);
phi{2} = pi*(1+5^(1/2))*indices';
Vibdir{2} = [cos(phi{2}).*sin(theta{2}),sin(phi{2}).*sin(theta{2}),cos(theta{2})];

N = 2000;
indices = [0:N-1]+0.5;
theta{3} = acos(1-2*indices'/N);
phi{3} = pi*(1+5^(1/2))*indices';
Vibdir{3} = [cos(phi{3}).*sin(theta{3}),sin(phi{3}).*sin(theta{3}),cos(theta{3})];

N = 3000;
indices = [0:N-1]+0.5;
theta{4} = acos(1-2*indices'/N);
phi{4} = pi*(1+5^(1/2))*indices';
Vibdir{4} = [cos(phi{4}).*sin(theta{4}),sin(phi{4}).*sin(theta{4}),cos(theta{4})];

N = 100;
u = rand(N,1);
v = rand(N,1);
nr_rand = 5;
phi{nr_rand} = 2*pi*u;
theta{nr_rand} = acos(2*v-1);
Vibdir{nr_rand} = [cos(phi{nr_rand}).*sin(theta{nr_rand}),sin(phi{nr_rand}).*sin(theta{nr_rand}),cos(theta{nr_rand})];

nr_ord=6;
Vibdir_inter = [];
phi_inter = [];
theta_inter = [];
for itheta = linspace(0,pi,10)
    dphi = pi/10;
    if itheta==0 || itheta==pi
        dphi = 2*pi;
    end
    for iphi=[0:dphi:2*pi-dphi]
        phi_inter = [phi_inter;iphi];
        theta_inter = [theta_inter;itheta];
        Vibdir_inter = [Vibdir_inter ;cos(iphi)*sin(itheta),sin(iphi)*sin(itheta),cos(itheta)];
    end
end
Vibdir{nr_ord} = Vibdir_inter;
phi{nr_ord} = phi_inter;
theta{nr_ord} = theta_inter;
SolutionType = '3SphereS8.2R25';
dA0all = {};
% Declaration POIs
POIs = [];
for itheta = linspace(0,pi,409)
%     dphi = pi/4;
%     if itheta==0 || itheta==pi
%         dphi = 2*pi;
%     end
    for iphi=0%[0:dphi:2*pi-dphi]
        POIs = [POIs;cos(iphi)*sin(itheta),sin(iphi)*sin(itheta),cos(itheta)];
    end
end
if contains(SolutionType,'Mouse')
    Rmax = 4.1;
    dR=0.1;
    Aus_max = 5;
else
    Rmax = 65;
    dR=1;
    Aus_max = 6;
end
dp_pos2p = [0,0,0;0,0,Rmax]; dp_pos2p = dp_pos2p/1000;
orien_dp2p = [0,0,1;0,0,1];
Aus = 1e-6;
%Aus2p=linspace(1e-3,1e-2,20);
Ifun = @(t) 1+0*t;
fus = 1e6;
%%
for icoll=1:nr_ord
    disp([sprintf('set %i: ',icoll),datestr(now,'HH:mm:ss')])
    Pflag = 0;
    for idp=1:size(dp_pos2p,1)
        parfor ivib=1:length(theta{icoll})
            Dirus = Vibdir{icoll}(ivib,:);            
            
            [~,~,~,Out1] = investBiologicalNoise('posDp',dp_pos2p(idp,:),'POIs',POIs,'OrienDipole',orien_dp2p(idp,:),...
                'Plot',0,'Display',0,'totaldps',1,'dps_run',1,'showSphere','','Solutiontype',SolutionType,...
                'Ifun',Ifun,'fus',fus,'Aus',Aus,'Dirus',Dirus,'Tend',2*fus^-1,'noRecon',1);
            AmpVR1 = (max(Out1.VR,[],2)-min(Out1.VR,[],2))/2*10^6;
            dA0mPOI(ivib,:) = AmpVR1;
        end
        dA0all{idp,icoll} = dA0mPOI;
        dA0mPOI = [];
        Pval = round(((icoll-1)*size(dp_pos2p,1)+idp)/(size(dp_pos2p,1)*nr_ord)*100,0);
        if Pval >= Pflag
            disp(datetime('now','Format','HH:mm:ss'))
            disp(['Progress: ',num2str(Pval),'%'])
            Pflag = Pflag+5;
        end
    end
end
save(['OutputLinearityvib_',datestr(now,'yyyyMMdd-hhmm'),'.mat'],'dA0all','theta','phi','Vibdir','POIs');
%% compare difference between interpolation
load('Linearity/OutputLinearityvib_20203807-2201.mat')
addpath(genpath('D:\Users\rschoeters\Documents\MATLAB\Interestingfunctions'))
perc = 5/100; %error lim
abs_tol = [1e-3,1]; %abs tol in pV (indicative for small signal strengths and perhaps numerical errors)

for idp=1:size(dp_pos2p,1)
    for icoll = 1:nr_rand-1
        [rdAinterp_rand{idp,icoll},used_points_randall{idp,icoll}] = slerp2_sph(phi{icoll},theta{icoll},dA0all{idp,icoll},phi{nr_rand},theta{nr_rand});
        adA_rand{idp,icoll} = abs(dA0all{idp,nr_rand} - rdAinterp_rand{idp,icoll});
        rdA_rand{idp,icoll} = adA_rand{idp,icoll}./ abs(dA0all{idp,nr_rand});
        tr_rand(idp,icoll) = sum(sum(rdA_rand{idp,icoll}>perc));
        for iabstol = 1:length(abs_tol)
            tq0_rand(idp,icoll) = sum(sum(adA_rand{idp,icoll}<abs_tol(iabstol)));
            
            trnq0_rand(idp,icoll,iabstol) = sum(sum(rdA_rand{idp,icoll}(dA0all{idp,nr_rand}>abs_tol(iabstol) & ...
                rdAinterp_rand{idp,icoll}>abs_tol(iabstol))>perc));
        end
        
        [rdAinterp_ord{idp,icoll},used_points_ordall{idp,icoll}] = slerp2_sph(phi{icoll},theta{icoll},dA0all{idp,icoll},phi{nr_ord},theta{nr_ord});
        adA_ord{idp,icoll} = abs(dA0all{idp,nr_ord} - rdAinterp_ord{idp,icoll});
        rdA_ord{idp,icoll} = adA_ord{idp,icoll}./ abs(dA0all{idp,nr_ord});
        tr_ord(idp,icoll) = sum(sum(rdA_ord{idp,icoll}>perc));
        for iabstol = 1:length(abs_tol)
            tq0_ord(idp,icoll) = sum(sum(adA_ord{idp,icoll}<abs_tol(iabstol)));
            
            trnq0_ord(idp,icoll,iabstol) = sum(sum(rdA_ord{idp,icoll}(dA0all{idp,nr_ord}>abs_tol(iabstol) & ...
                rdAinterp_ord{idp,icoll}>abs_tol(iabstol))>perc));
        end
        
    end
end
%% compare difference between interpolation other interpolation method

perc = 5/100;
abs_tol = 1e-9;
for idp=1:size(dp_pos2p,1)
    for icoll = 1:nr_rand-1        
        rdAinterp_ordnnb{idp,icoll} = slerp2_sph(phi{icoll},theta{icoll},dA0all{idp,icoll},phi{nr_ord},theta{nr_ord},'method','nnb');
        adA_ordnnb{idp,icoll} = abs(dA0all{idp,nr_ord} - rdAinterp_ordnnb{idp,icoll});
        rdA_ordnnb{idp,icoll} = adA_ordnnb{idp,icoll}./ abs(dA0all{idp,nr_ord});
        tr_ordnnb(idp,icoll) = sum(sum(rdA_ordnnb{idp,icoll}>perc));
        tq0_ordnnb(idp,icoll) = sum(sum(adA_ordnnb{idp,icoll}<abs_tol));
        trnq0_ordnnb(idp,icoll) = sum(sum(rdA_ordnnb{idp,icoll}(dA0all{idp,nr_ord}>abs_tol & ...
            rdAinterp_ordnnb{idp,icoll}>abs_tol)>perc));
    end
    
end
%% compare difference between interpolation
%othe rinterpolation method
perc = 5/100;
abs_tol = 1e-9;
for idp=1:size(dp_pos2p,1)
    for icoll = 1:nr_rand-1        
        rdAinterp_ordinterp2{idp,icoll} = slerp2_sph(phi{icoll},theta{icoll},dA0all{idp,icoll},phi{nr_ord},theta{nr_ord},'method','slerp2');
        adA_ordinterp2{idp,icoll} = abs(dA0all{idp,nr_ord} - rdAinterp_ordinterp2{idp,icoll});
        rdA_ordinterp2{idp,icoll} = adA_ordinterp2{idp,icoll}./ abs(dA0all{idp,nr_ord});
        tr_ordinterp2(idp,icoll) = sum(sum(rdA_ordinterp2{idp,icoll}>perc));
        tq0_ordinterp2(idp,icoll) = sum(sum(adA_ordinterp2{idp,icoll}<abs_tol));
        trnq0_ordinterp2(idp,icoll) = sum(sum(rdA_ordinterp2{idp,icoll}(dA0all{idp,nr_ord}>abs_tol & ...
            rdAinterp_ordinterp2{idp,icoll}>abs_tol)>perc));
    end
    
end
%%
Colors= [        0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

figure
nr = 2;
scatter3(Vibdir{nr}(:,1),Vibdir{nr}(:,2),Vibdir{nr}(:,3))
title('DataBase: Uniform distributed (N=1000)')
axis square

figure
nr = 5;
scatter3(Vibdir{nr}(:,1),Vibdir{nr}(:,2),Vibdir{nr}(:,3))
title('Test: Random (N=100)')
axis square

figure
nr = 6;
scatter3(Vibdir{nr}(:,1),Vibdir{nr}(:,2),Vibdir{nr}(:,3))
title('Test: Ordered (N=162)')
axis square

Nvals = [100,1000,2000,3000];
Names = {'[0,0,0]mm','[0,0,65]mm'};
types = {'o--','*-.'}
abstol= arrayfun(@num2str,abs_tol,'UniformOutput',false);
figure
for iabstol = 1:2
for i=1:size(trnq0_ord,1)
        subplot(1,2,1)
        hold on
        plot(Nvals,trnq0_rand(i,:,iabstol)/numel(rdA_rand{1})*100,types{iabstol},'Color',Colors(i,:),'Displayname',[Names{i},',  ',abstol{iabstol},'pV'])
        hold off
        ylabel('# errors (%)')
        xlabel('Database size')
        title('test: rand')
%         leg = legend('show','location','southoutside');
%         title(leg,'dipole positon & abs tol')
%         
        subplot(1,2,2)
        hold on
        plot(Nvals,trnq0_ord(i,:,iabstol)/numel(rdA_ord{1})*100,types{iabstol},'Color',Colors(i,:),'Displayname',[Names{i},',  ',abstol{iabstol},'pV'])
        hold off
        ylabel('# errors (%)')
        xlabel('Database size')
        title('test: ordered')
        leg = legend('show','location','eastoutside');
        title(leg,'dipole positon & abs tol')
end
end
%% Check if correct interpolation
Vibdiridx = 45;
POIidx = 311;
InterpolPoints = [used_points_ordall{2,4}{1,1}(Vibdiridx,:);...
    used_points_ordall{2,4}{1,2}(Vibdiridx,:);used_points_ordall{2,4}{1,3}(Vibdiridx,:);...
    used_points_ordall{2,4}{1,4}(Vibdiridx,:)];
figure
nr = 6;
scatter3(Vibdir{nr}(:,1),Vibdir{nr}(:,2),Vibdir{nr}(:,3))
hold on
scatter3(Vibdir{nr}(Vibdiridx,1),Vibdir{nr}(Vibdiridx,2),Vibdir{nr}(Vibdiridx,3),'rx')
scatter3(POIs(POIidx,1),POIs(POIidx,2),POIs(POIidx,3),'ks')
scatter3(InterpolPoints(:,1),InterpolPoints(:,2),InterpolPoints(:,3),'g*')
hold off
title('Test: Ordered (N=162)')
axis square


