% Compare different models:
% different model type, different layer thickness
% All results are for constant moving dipole with variable radial
% orientation and radial position. Dipole is thus always radially oriented
if ~isdeployed
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Inputs'))
end
if exist('Fignr','var')
    clearvars -except 'Fignr'
    Fignr = Fignr+1;
else
    clear all
    Fignr = 10;
end
CM = 'jet';     %Colormap
equalColor = 1  %equal colorbar along plots
PLOT_DC = 0;    %plot DC map (not oscillation effect)
% Electrical properties
Npercolumn = 1e4;
I = 1*1e-3*Npercolumn; %µA
sigma = 0.33;    %S/m
d = 500*10^-6; %distance between current source and current sink
dI = I*d;

% US options
fus = 1e6; %Hz
wus = 2*pi*fus; %rad/s
Aus = 10e-9; %m
Phaseus = 0; %rad
Dirus = [0,0,1]; Dirus = Dirus/norm(Dirus);
USwave = @(t,idx) Dirus.*Aus.*sin(wus.*t+Phaseus);

% POIs
Theta = [0:pi/50:pi];
dPhi = pi/10;
[phi,theta]=meshgrid(0:dPhi:pi,Theta);
phisc = phi(:); thetasc = theta(:);

% Simulation settings
resUS = 20;
dt = (resUS*fus)^-1;
Tend = 1*fus^-1;
Tsim = 0:dt:Tend;
resolution = 10;


%model sets
%MODELS = {'ClosedBounded2','3SphereS7R1','4SphereasClosedBounded2'}
%MODELS = {'4SphereasClosedBounded2'}

%MODELS = {'3SphereS8.2R1','3SphereS8.2R25','3SphereS8.2R~f'};
%MODELS = {'4SphereS8.7R~f'};
MODELS = {'3SphereS8.2R25','3SphereS8.2R~f','4SphereS8.7R25','4SphereS8.7R~f','Mouse4Sphere~fair0','Mouse4Sphere~fair1'};
%MODELS  ={'3SphereS8.7R25_{VAL}','4SphereS8.7R25_{VAL}','4SphereS8.7R~f_{VAL}'};



tic
%% generate plots for dipole with variable distance to sphere + variable orientation to vibration
for imodel = 1:length(MODELS)
    [Options,RSphere,RPOI,RTICKS,Rvals,Angles] = getSettings(MODELS{imodel},resolution);
    VR = [];
    VRfluc = [];
    maxVRfluc = [];
    Pflag = 0;
    POIs = RPOI.*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];
    
    [X,Y] = meshgrid(Angles,Rvals);
    X_coll{imodel} = X; Y_coll{imodel} = Y;  %collection strings for plots
    Rvals_coll{imodel} = Rvals;
    Angles_coll{imodel} = Angles;
    RPOI_coll{imodel} = RPOI;
    RTICKS_coll{imodel} = RTICKS;

    for iR = 1:length(Rvals)
        for iorien = 1:length(Angles)
            % Select Sources and sinks + select vibrator
            OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
            LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
            CSource = LocDipole+d/2*OrienDipole;
            CSink = LocDipole-d/2*OrienDipole; 
            Settings = horzcat(Options,{'POI',POIs,'resUS',resUS,'Display',double(iR==1&&iorien==1),'scale',0});
            
            % Run simulation
            VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,1,Settings);
            
            if contains(lower(MODELS{imodel}),'mouse')
                VRfluc = (max(VR,[],2)-min(VR,[],2))/2*10^3;
            else
                VRfluc = (max(VR,[],2)-min(VR,[],2))/2*10^6;
            end            
            VRmean = (max(VR,[],2)+min(VR,[],2))/2;
            maxVRmean(iR,iorien) = max(VRmean);          
            maxVRfluc(iR,iorien) = max(VRfluc);
            maxVRmean_coll{imodel} = maxVRmean;
            maxVRfluc_coll{imodel} = maxVRfluc;
        
            
        end
        PComp = round(iR/length(Rvals)*100);
        if PComp>=Pflag
            disp(['progress:',num2str(PComp),'%'])
            Pflag=Pflag+1;
        end
    end

end
toc
%% Plot
Fignr = get(gcf,'number');
for imodel = 1:length(MODELS)
 X = X_coll{imodel}; Y = Y_coll{imodel};
 maxVRfluc = maxVRfluc_coll{imodel};  maxVRmean = maxVRmean_coll{imodel};
 Rvals = Rvals_coll{imodel}; RPOI = RPOI_coll{imodel}; RTICKS = RTICKS_coll{imodel};
    figure(Fignr+1)
    if imodel == 1
        %set(gcf,'position',[-1919,41,1920,963])
    end
    %save subplots: easier to alter properties
    sp_fluc{imodel} = subplot(2+PLOT_DC+(1+PLOT_DC)*floor((length(MODELS)-0.1)/3),min(length(MODELS),3),imodel);
    surf(X,Y,maxVRfluc,'EdgeColor','none','FaceColor','interp')
    xlabel('Clockwise rotation in yz plane')
    ylabel('Distance dipole from outer sphere [cm]')
    zlabel('Fluctuation amplitude [pV]')
    set(gca,'XTick',0:pi/4:pi)
    set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
    set(gca,'YTick',RTICKS)
    set(gca,'YTickLabel',arrayfun(@num2str,(RPOI-get(gca,'YTick'))*100,'UniformOutput',false))
    colormap('jet')
    xlim([0,pi]);
    ylim([0,RPOI]);
    c = colorbar;
    if any(contains(lower(MODELS),'mouse'))
        if contains(lower(MODELS{imodel}),'mouse')
            c.Label.String = 'Fluctuation Amplitude [nV]';
        else
            c.Label.String = 'Fluctuation Amplitude [pV]';
        end
    else
        if imodel == length(MODELS)
            c.Label.String = 'Fluctuation Amplitude [pV]';
        end
    end
    if ~equalColor
    set(gca,'colorscale','log')
    c.Ticks = [min(maxVRfluc(:)),c.Ticks,max(maxVRfluc(:))];
    c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
    end
    title(['Maximal Amplitude: ',MODELS{imodel}])
    view(0,90)
    LIMS_fluc(imodel,:) =  [min(maxVRfluc(:)),max(maxVRfluc(:))];
    
    if PLOT_DC
        sp_DC{imodel} = subplot(2+PLOT_DC+(1+PLOT_DC)*floor((length(MODELS)-0.1)/3),min(length(MODELS),3),(1+floor((length(MODELS)-0.1)/3))*min(length(MODELS),3)+imodel);
        surf(X,Y,maxVRmean,'EdgeColor','none','FaceColor','interp')
        xlabel('Clockwise rotation in yz plane')
        ylabel('Distance dipole from outer sphere [cm]')
        zlabel('Fluctuation amplitude [µV]')
        set(gca,'XTick',0:pi/4:pi)
        set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
        set(gca,'YTick',RTICKS)
        set(gca,'YTickLabel',arrayfun(@num2str,(RPOI-get(gca,'YTick'))*100,'UniformOutput',false))
        colormap('jet')
        xlim([0,pi]);
        ylim([0,RPOI]);
        c = colorbar;
        if imodel == length(MODELS)
            c.Label.String = 'Fluctuation Amplitude [µV]';
        end
        if ~equalColor
        set(gca,'colorscale','log')
        c.Ticks = [min(maxVRmean(:)),c.Ticks,max(maxVRmean(:))];
        c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
        end
        title(['Maximal Amplitude: ',MODELS{imodel}])
        view(0,90)
        LIMS_DC(imodel,:) =  [min(maxVRmean(:)),max(maxVRmean(:))];
    end
    
    subplot(2+PLOT_DC+(1+PLOT_DC)*floor((length(MODELS)-0.1)/3),min(length(MODELS),3),[(1+PLOT_DC)*(1+floor((length(MODELS)-0.1)/3))*min(length(MODELS),3)+1:(1+(1+PLOT_DC)*(1+floor((length(MODELS)-0.1)/3)))*min(length(MODELS),3)]);
    hold on
    %plot(RPOI-Rvals,fliplr(maxVRfluc(:,1)),'Displayname',MODELS{imodel})
    plot(Rvals,maxVRfluc(:,1),'Displayname',MODELS{imodel})
    hold off
    pause(0.1)
end
legend('show','location','eastoutside')
%xlabel('Distance dipole from outer sphere [cm]')
xlabel('dipole position from centre [cm]')
ylabel({'Maximal measurable amplitude [pV]','dipole parallel with oscillation'})
set(gca,'yscale','log')%,'xtick',RTICKS)
set(findobj('type','axes'),'fontsize',12)

if equalColor
    idx_mouse = contains(lower(MODELS),'mouse');
    idx_human = ~idx_mouse;
    if any(idx_human)
        LIMS_fluc_h = [min(LIMS_fluc(idx_human,1)),max(LIMS_fluc(idx_human,2))];
        for isp = find(idx_human)
            c = sp_fluc{isp}.Colorbar;
            set(sp_fluc{isp},'CLIM',LIMS_fluc_h)
            set(sp_fluc{isp}.Colorbar,'Limits',LIMS_fluc_h);
            set(sp_fluc{isp},'colorscale','log')
            ticks = sort([LIMS_fluc(isp,1),c.Ticks,LIMS_fluc(isp,2)]); ticks(ticks>LIMS_fluc(isp,2)) = [];
            %c.Ticks = unique([round(LIMS_fluc_h(1),2,'significant'),ticks,LIMS_fluc_h(2)]);
            c.Ticks = ticks;
            c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
        end
        if PLOT_DC
            LIMS_DC_h = [min(LIMS_DC(idx_human,1)),max(LIMS_DC(idx_human,2))];
            for isp = find(idx_human)
                c = sp_DC{isp}.Colorbar;
                set(sp_DC{isp},'CLIM',LIMS_DC_h)
                set(sp_DC{isp}.Colorbar,'Limits',LIMS_DC_h);
                set(sp_DC{isp},'colorscale','log')
                ticks = sort([LIMS_DC(isp,1),c.Ticks,LIMS_DC(isp,2)]); ticks(ticks>LIMS_DC(isp,2)) = [];
                %c.Ticks = unique([round(LIMS_fluc_h(1),2,'significant'),ticks,LIMS_fluc_h(2)]);
                c.Ticks = ticks;
                c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
            end
        end
    end
    if any(idx_mouse)
         LIMS_fluc_m = [min(LIMS_fluc(idx_mouse,1)),max(LIMS_fluc(idx_mouse,2))];
        for isp = find(idx_mouse)
            c = sp_fluc{isp}.Colorbar;
            set(sp_fluc{isp},'CLIM',LIMS_fluc_m)
            set(sp_fluc{isp}.Colorbar,'Limits',LIMS_fluc_m);
            set(sp_fluc{isp},'colorscale','log')
            ticks = sort([LIMS_fluc(isp,1),c.Ticks,LIMS_fluc(isp,2)]); ticks(ticks>LIMS_fluc(isp,2)) = [];
            %c.Ticks = unique([round(LIMS_fluc_h(1),2,'significant'),ticks,LIMS_fluc_h(2)]);
            c.Ticks = ticks;
            c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
        end
        if PLOT_DC
            LIMS_DC_m = [min(LIMS_DC(idx_mouse,1)),max(LIMS_DC(idx_mouse,2))];
            for isp = find(idx_mouse)
                c = sp_DC{isp}.Colorbar;
                set(sp_DC{isp},'CLIM',LIMS_DC_m)
                set(sp_DC{isp}.Colorbar,'Limits',LIMS_DC_m);
                set(sp_DC{isp},'colorscale','log')
                ticks = sort([LIMS_DCisp,1,c.Ticks,LIMS_DC(isp,2)]); ticks(ticks>LIMS_DC(isp,2)) = [];
                %c.Ticks = unique([round(LIMS_fluc_h(1),2,'significant'),ticks,LIMS_fluc_h(2)]);
                c.Ticks = ticks;
                c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
            end
        end
    end
end

