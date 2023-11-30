% determine fluctution effect based on elevation and azimuth
%MINI COLUMN
% Run simulations for moving dipole
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
if exist('Fignr','var')
    clearvars -except 'Fignr'
    Fignr = Fignr+1;
else
    clear all
    Fignr = 10;
end
CM = 'jet';
GENGIF = false;


% Electrical properties
Npercolumn = 1e2;
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

% Simulation settings
resUS = 10;
dt = (resUS*fus)^-1;
Tend = 2*fus^-1;
Tsim = 0:dt:Tend;

% Geometric  properties
SolutionType = '4Sphere~f';
if strcmpi(SolutionType,'3Sphere')
    RSphere = 0.058; %m
    RPOI = RSphere + 0.012;
    Rvals = linspace(0,RPOI-0.002,21);
    RTICKS = 0:0.01:RPOI; RTICKS(end:end+1) = [Rvals(end),RTICKS(end)];
elseif strcmpi(SolutionType,'4Sphere~f')
    RSphere = 0.07;
    RPOI = RSphere+0.017;
    SolutionType = '4Sphere';
    RatioSkT = 1/25;
    Rvals = linspace(0,RSphere-0.002,resolution+1);
    RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
    fDependence = 1;
    scale_flag = 0;
else
    RSphere = 0.07; %m
    RPOI = 0.07;
    Rvals = linspace(0,RPOI-0.002,21);
    RTICKS = 0:0.01:RPOI; RTICKS(end:end+1) = [Rvals(end),RTICKS(end)];
    
end
[xs,ys,zs] = sphere(100);
xsphere = xs*RSphere; ysphere = ys*RSphere; zsphere = zs*RSphere;
Inputtype = '';

% dipoles on a sphere
Rd = 0.03; % concentric sphere on which dipoles are situated
dTheta = pi/4;
dPhi = pi/4;
xd = []; yd = []; zd = [];
idx = 0;
LocDipoleRd=[];
for i=[0:dTheta:pi]
    if i==0 || i==pi
        Angles = 0;
    else
        Angles = [dPhi:dPhi:2*pi];
    end
    for j= 1:length(Angles)
        idx=idx+1;
        LocDipoleRd(idx,:) = Rd*[cos(Angles(j))*sin(i),sin(Angles(j))*sin(i),cos(i)];
    end
end

% POIs

Theta = [0:pi/100:pi];
dPhi = pi/4;
[phi,theta]=meshgrid(0:dPhi:pi,Theta);
phisc = phi(:); thetasc = theta(:);       
 POIs = RSphere*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];
    


%% dipole in center and yz-direction
tic
OrienDipole = [0,1/sqrt(2),1/sqrt(2)];
idxOscDip = 1;

% Select Sources and sinks + select vibrator
CSource = (0.01+d/2)*OrienDipole;
CSink = (0.01-d/2)*OrienDipole;

% Run simulation
VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,'SolutionType',SolutionType,'resUS',resUS,'scale',scale_flag);
VRfluc = (max(VR,[],2)-min(VR,[],2))*10^6;

% PLOT RESULTS
VRflucSurf = reshape(VRfluc,size(phi));
figure(Fignr+2)
surf(phi,theta,VRflucSurf/2,'EdgeColor','none')
xlabel('Azimuth')
ylabel('Co-elevation')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/2:2*pi) 
 set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
 set(gca,'YTick',0:pi/4:pi) 
 set(gca,'YTickLabel',{'0','pi/4','pi/2','3*pi/4','pi'})
 colormap('jet')
 xlim([0,2*pi]);
  ylim([0,pi]);
 Fignr = get(gcf,'number');
Fignr = Fignr+1;
toc
%% Multiple orientation with and or without scaling

OrienDipoles = [0,0,1;0,1,0;0,1/sqrt(2),1/sqrt(2);0,1/sqrt(2),-1/sqrt(2)];
idxOscDip = 1;
figure(Fignr)
sp = {};
c = {};
iscale_end = 1;
for iscale = 0:iscale_end
for ifig = 1:size(OrienDipoles,1)
sp{ifig} = subplot(1+iscale_end,size(OrienDipoles,1),iscale*size(OrienDipoles,1)+ifig);
OrienDipole = OrienDipoles(ifig,:);
% Select Sources and sinks + select vibrator
CSource = d/2*OrienDipole;
CSink = -d/2*OrienDipole;
% Run simulation
Tsim = 0:dt:Tend;
VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'resUS',resUS,'POI',POIs,'SolutionType',SolutionType,'scale',iscale);
VRfluc = (max(VR,[],2)-min(VR,[],2))*10^6;

%PLOT 
VRflucSurf = reshape(VRfluc,size(phi));
surf(phi,theta,VRflucSurf/2,'EdgeColor','none')
xlabel('Azimuth')
ylabel('Co-elevation')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/2:2*pi) 
 set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'})
 set(gca,'YTick',0:pi/4:pi) 
 set(gca,'YTickLabel',{'0','pi/4','pi/2','3*pi/4','pi'})
 colormap('jet')
 c{ifig} = colorbar;
 c{ifig}.Label.String = 'Fluctuation Amplitude [pV]';
 xlim([0,2*pi]);
 ylim([0,pi]);
title(['Dipole [',num2str(OrienDipole),']'])
view(0,90)
end
for ifig = 1:size(OrienDipoles,1)
    LIMS(ifig,:) = get(c{ifig},'Limits');
end
LIMS = [min(LIMS(:,1)),max(LIMS(:,2))];
for ifig = 1:size(OrienDipoles,1)
    set(sp{ifig},'CLIM',LIMS);
    set(sp{ifig}.Colorbar,'Limits',LIMS);
    if ifig<size(OrienDipoles,1)
    set(sp{ifig}.Colorbar,'Visible', 'off')
    end
end
end
Fignr = get(gcf,'number');
Fignr = Fignr+1;

%% generate plots for dipole with variable distance to sphere + variable orientation to vibration
%four plots one in direction of vibration, one in direction y-axis one in
%direction oriention dipole and one in direction inbetween
%vibration still in z-direction

VR = [];
VRfluc = [];
idxOscDip = 1;
dAngle = pi/(10*resolution);
Angles = 0:dAngle:pi;
Pflag = 0;
for iR = 1:length(Rvals)
    parfor iorien = 1:length(Angles)
        % Select Sources and sinks + select vibrator
        OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
        LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
        CSource = LocDipole+d/2*OrienDipole;
        CSink = LocDipole-d/2*OrienDipole;
        
        % point of interests
        MeanDirection = ([0,0,1]+OrienDipole)/2;
        if Angles(iorien)<pi/2
            OptimalPOI = MeanDirection;
        else
            OptimalPOI = ([0,0,-1]+OrienDipole)/2;
        end
        POIs = RSphere.*[0,0,1;0,1,0;OrienDipole;OptimalPOI./norm(OptimalPOI)];
        
        % Run simulation
        Tsim = 0:dt:Tend;
        
        VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,'Display',0,...
            'resUS',resUS,'RSphere',RSphere,'SolutionType',SolutionType);
        
        VRfluc(iR,iorien,:) = (max(VR,[],2)-min(VR,[],2))/2*10^6;
    end
    PComp = round(iR/length(Rvals)*100);
    if PComp>=Pflag
        disp([num2str(PComp),'%'])
        Pflag=Pflag+1;
    end
end
[X,Y] = meshgrid(Angles,Rvals);
figure(Fignr+1)
sp={};
c={};
titles={'[0 0 1]','[0 1 0]','Orientation p','mean z-axis and orientation dipole'};
for ifig=1:4
    figure(Fignr+1)
sp{ifig}=subplot(ceil(4/2),2,ifig);
surf(X,Y,VRfluc(:,:,ifig),'EdgeColor','none')
xlabel('Clockwise rotation in yz plane')
ylabel('Distance dipole form centre [cm]')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',RTICKS) 
set(gca,'YTickLabel',arrayfun(@num2str,(RPOI-get(gca,'YTick'))*100,'UniformOutput',false))
colormap('jet')
xlim([0,pi]);
ylim([0,RPOI]);
set(gca,'colorscale','log')
c{ifig} = colorbar;
c{ifig}.Label.String = 'Fluctuation Amplitude [pV]';
title(['POI : ',titles{ifig}])
view(0,90)
figure(Fignr+2)
hold on
plot(X(floor(end/2),:),VRfluc(floor(end/2),:,ifig),'DisplayName',titles{ifig})
xlabel('Clockwise rotation in yz plane')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','pi/4','pi/2','3*pi/4','pi'})
colormap('jet')
xlim([0,pi]);
hold off
legend('show')
end
LIMS=[];
for ifig = 1:4
    LIMS(ifig,:) = get(c{ifig},'Limits');
end
LIMS = [min(LIMS(:,1)),max(LIMS(:,2))];
for ifig = 1:4
    set(sp{ifig},'CLIM',LIMS);
    set(sp{ifig},'ZLIM',LIMS);
    set(sp{ifig}.Colorbar,'Limits',LIMS);
    c{ifig}.Ticks = [min(c{ifig}.Limits),c{ifig}.Ticks,max(c{ifig}.Limits)];
    c{ifig}.TickLabels = arrayfun(@num2str,c{ifig}.Ticks,'UniformOutput',false);
    if mod(ifig,2)
    set(sp{ifig}.Colorbar,'Visible', 'off')
    end
end
drawnow
pause(0.1)
Fignr = get(gcf,'number');
Fignr = Fignr+1;
%% generate plots for dipole with variable distance to sphere + variable orientation to vibration
%3 plots one with maximal amplitude, one with optimal azimuth and one with
% optimal co-elevation
%vibration still in z-direction
tic
VR = [];
VRfluc = [];
maxVRfluc = [];
idxOscDip = 1;
PHI=[];
THETA=[];
Pflag = 0;
% POIs

Theta = [0:pi/100:pi];
dPhi = pi/4;
[phi,theta]=meshgrid(0:dPhi:pi,Theta);
phisc = phi(:); thetasc = theta(:);
POIs=[];
dAngle = pi/(10*resolution);
Angles = 0:dAngle:pi;
iscale_end=2;
Show_sphere={'',''};
titles = {'POIs on Brain','POIs on Scalp','POIs 5mm from scalp'};
for iscale=0:iscale_end

for iR = 1:length(Rvals)
    parfor iorien = 1:length(Angles)
        % Select Sources and sinks + select vibrator
        OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
        LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
        CSource = LocDipole+d/2*OrienDipole;
        CSink = LocDipole-d/2*OrienDipole;
        
        POIs = RSphere.*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];
        if iscale == 1
            POIs = POIs/RSphere*(RSphere+0.012);
        end
        
        % Run simulation
        Tsim = 0:dt:Tend;
        VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,...
            'resUS',resUS,'Display',double(iorien==1&&iR==1),'SolutionType',SolutionType,'scale',floor((2-iscale)/2),'ShowSphere',Show_sphere{double(iR==1&&iorien==1)+1});
        VRfluc = (max(VR,[],2)-min(VR,[],2))/2*10^6;
        
        [maxVRfluc(iR,iorien),idxmaxVRfluc]=max(VRfluc);
        if POIs(idxmaxVRfluc,1) <= 1e-3 && (POIs(idxmaxVRfluc,2)~=0||POIs(idxmaxVRfluc,3)~=0)
            PHI(iR,iorien) = pi/2;
        elseif POIs(idxmaxVRfluc,1) <= 1e-3 && POIs(idxmaxVRfluc,2) <= 1e-3
            PHI(iR,iorien) = pi/2;
        else
            PHI(iR,iorien) = atan(POIs(idxmaxVRfluc,2)/POIs(idxmaxVRfluc,1))+...
                double(POIs(idxmaxVRfluc,1)<0)*pi+double(POIs(idxmaxVRfluc,1)>0&&POIs(idxmaxVRfluc,2)<0)*2*pi;
        end
        if Angles(iorien) == 0 || Angles(iorien) == pi
            THETA(iR,iorien) = Angles(iorien);
        else
            THETA(iR,iorien) = acos(POIs(idxmaxVRfluc,3)/RSphere);
        end

    end
    PComp = round(iR/length(Rvals)*100);
    if PComp>=Pflag
        disp([num2str(PComp),'%'])
        Pflag=Pflag+1;
    end
end

[X,Y] = meshgrid(Angles,Rvals);
figure(Fignr+1)
sp={};
c={};
figure(Fignr+1)
subplot(1+iscale_end,1+double(iscale_end==0)*2,1+iscale);
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
c.Label.String = 'Fluctuation Amplitude [pV]';
set(gca,'colorscale','log')
c.Ticks = [min(maxVRfluc(:)),c.Ticks,max(maxVRfluc(:))];
c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
title(['Maximal Amplitude: ', titles{3-iscale}])
view(0,90)
if iscale_end==0
subplot(1,3,2)
surf(X,Y,PHI,'EdgeColor','none','FaceColor','interp')
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
c=colorbar;
title('Azimuth optimal POI')
c.Ticks = [0:pi/2:pi];
c.TickLabels = {'0','\pi/2','\pi'};
view(0,90)
subplot(1,3,3)
surf(X,Y,THETA,'EdgeColor','none')
xlabel('Clockwise rotation in yz plane')
ylabel('Distance dipole from oute sphere [cm]')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',RTICKS) 
set(gca,'YTickLabel',arrayfun(@num2str,(RPOI-get(gca,'YTick'))*100,'UniformOutput',false))
colormap('jet')
xlim([0,pi]);
ylim([0,RPOI]);
c = colorbar;
c.Ticks = [0:pi/4:pi];
c.TickLabels = {'0','\pi/4','\pi/2','3\pi/4','\pi'};
title('Co-elevation optimal POI')
view(0,90)
end
end
Fignr = get(gcf,'number');
Fignr = Fignr+1;
toc
