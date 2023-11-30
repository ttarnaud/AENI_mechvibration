% determine fluctution effect based on elevation and azimuth
%MINI COLUMN
% Run simulations for moving dipole
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
dt = (1e2*fus)^-1;
Tend = 2*1e-6;

% Geometric  properties
RSphere = 0.07;  %m
[xs,ys,zs] = sphere(1000);
xsphere = xs*RSphere; ysphere = ys*RSphere; zsphere = zs*RSphere;
Inputtype = '';
SolutionType = 'ClosedBounded';

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
dPhi = pi/100;
[phi,theta]=meshgrid(0:dPhi:2*pi,Theta);
phisc = phi(:); thetasc = theta(:);
POIs=[];

    
        
 POIs = RSphere*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];
    


%% dipole in center and yz-direction
OrienDipole = [0,1/sqrt(2),1/sqrt(2)];
idxOscDip = 1;
figure(Fignr)
ax0 = subplot(4,2,[1,3,5,7]);
%PotentialSphere_Multi(d/2*LocDipoleZ,-d/2*LocDipoleZ,I,sigma,RSphere,1000,Fignr);

% Select Sources and sinks + select vibrator
CSource = d/2*OrienDipole;
CSink = -d/2*OrienDipole;
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);
% Run simulation
Tsim = 0:dt:Tend;

for i=1:length(Tsim)
    t = Tsim(i);
    CSource(idxOscDip,:) = CSourceus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    CSink(idxOscDip,:) = CSinkus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    [VRMATSphere_tot,VRMATPOI_tot]=PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,1,[],'POI',POIs);
    VR(:,i) = VRMATPOI_tot;
    disp([num2str(round(i/length(Tsim)*100,2)),'%'])
end

%plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,10,16],VR)

VRfluc = (max(VR,[],2)-min(VR,[],2))*10^6;
% figure(Fignr+1)
% for i=1:size(VR,1)
%      plot3(i*ones(size(Tsim)),Tsim*1e6,(VR(i,:)-mean(VR(i,:)))*10^6)
%      hold on
%      VRfluc(i)=(max(VR(i,:))-min(VR(i,:)))*10^6;
% end
%  hold off
%  xlabel('POI')
%  ylabel('Time [µs]')
%  zlabel('Voltage [pV]')
%  figure(Fignr+2)
%  plot(1:length(VRfluc),VRfluc)
%  xlabel('POI')
%  ylabel('Potential variation [pV]')
VRflucSurf = reshape(VRfluc,size(phi));
figure(Fignr+1)
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
%% dipole in center and yz-direction

OrienDipoles = [0,0,1;0,1,0;0,1/sqrt(2),1/sqrt(2);0,1/sqrt(2),-1/sqrt(2)];
idxOscDip = 1;
figure(Fignr)
sp = {};
c = {};
for ifig = 1:size(OrienDipoles,1)
    sp{ifig} = subplot(1,size(OrienDipoles,1),ifig);
    OrienDipole = OrienDipoles(ifig,:);
% Select Sources and sinks + select vibrator
CSource = d/2*OrienDipole;
CSink = -d/2*OrienDipole;
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);
% Run simulation
Tsim = 0:dt:Tend;

for i=1:length(Tsim)
    t = Tsim(i);
    CSource(idxOscDip,:) = CSourceus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    CSink(idxOscDip,:) = CSinkus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    [VRMATSphere_tot,VRMATPOI_tot]=PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,1,[],'POI',POIs);
    VR(:,i) = VRMATPOI_tot;
    disp([num2str(round(i/length(Tsim)*100,2)),'%'])
end



VRfluc = (max(VR,[],2)-min(VR,[],2))*10^6;
% figure(Fignr+1)
% for i=1:size(VR,1)
%      plot3(i*ones(size(Tsim)),Tsim*1e6,(VR(i,:)-mean(VR(i,:)))*10^6)
%      hold on
%      VRfluc(i)=(max(VR(i,:))-min(VR(i,:)))*10^6;
% end
%  hold off
%  xlabel('POI')
%  ylabel('Time [µs]')
%  zlabel('Voltage [pV]')
%  figure(Fignr+2)
%  plot(1:length(VRfluc),VRfluc)
%  xlabel('POI')
%  ylabel('Potential variation [pV]')

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

Fignr = get(gcf,'number');
Fignr = Fignr+1;

%% generate plots for dipole with variable distance to sphere + variable orientation to vibration
%four plots one in direction of vibration, one in direction y-axis one in
%direction oriention dipole and one in direction inbetween
%vibration still in z-direction
resolution = 100;
VR = [];
VRfluc = [];
Pflag = 0;
idxOscDip = 1;
dR = RSphere/resolution;
dAngle = pi/(10*resolution);
Rvals = 0:dR:RSphere-0.005;
Angles = 0:dAngle:pi;
for iR = 1:length(Rvals)
for iorien = 1:length(Angles)
% Select Sources and sinks + select vibrator
OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
CSource = LocDipole+d/2*OrienDipole;
CSink = LocDipole-d/2*OrienDipole;
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);
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

for i=1:length(Tsim)
    t = Tsim(i);
    CSource(idxOscDip,:) = CSourceus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    CSink(idxOscDip,:) = CSinkus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    [VRMATSphere_tot,VRMATPOI_tot]=PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,1,[],'POI',POIs);
    VR(:,i) = VRMATPOI_tot;
    
end
VRfluc(iR,iorien,:) = (max(VR,[],2)-min(VR,[],2))/2*10^6;
PComp = round((iorien+(iR-1)*length(Angles))/(length(Rvals)*length(Angles))*100);
if PComp>=Pflag
disp([num2str(Pflag),'%'])
Pflag=Pflag+1;
end
end
end
%plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,10,16],VR)


% figure(Fignr+1)
% for i=1:size(VR,1)
%      plot3(i*ones(size(Tsim)),Tsim*1e6,(VR(i,:)-mean(VR(i,:)))*10^6)
%      hold on
%      VRfluc(i)=(max(VR(i,:))-min(VR(i,:)))*10^6;
% end
%  hold off
%  xlabel('POI')
%  ylabel('Time [µs]')
%  zlabel('Voltage [pV]')
%  figure(Fignr+2)
%  plot(1:length(VRfluc),VRfluc)
%  xlabel('POI')
%  ylabel('Potential variation [pV]')
[X,Y] = meshgrid(Angles,Rvals);
figure(Fignr+1)
sp={};
c={};
titles={'[0 0 1]','[0 1 0]','Orientation p','mean z-axis and orientation dipole'};
for ifig=1:size(POIs,1)
    figure(Fignr+1)
sp{ifig}=subplot(ceil(size(POIs,1)/2),2,ifig);
surf(X,Y,VRfluc(:,:,ifig),'EdgeColor','none')
xlabel('Clockwise rotation in yz plane')
ylabel('Distance dipole form centre [cm]')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',0:RSphere/10:RSphere) 
set(gca,'YTickLabel',arrayfun(@num2str,get(gca,'YTick')*100,'UniformOutput',false))
colormap('jet')
xlim([0,pi]);
ylim([0,RSphere]);
set(gca,'colorscale','log')
c{ifig} = colorbar;
c{ifig}.Label.String = 'Fluctuation Amplitude [pV]';
title(['POI : ',titles{ifig}])
view(0,90)
figure(Fignr+2)
hold on
plot(X(5,:),VRfluc(5,:,ifig),'DisplayName',titles{ifig})
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
for ifig = 1:size(POIs,1)
    LIMS(ifig,:) = get(c{ifig},'Limits');
end
LIMS = [min(LIMS(:,1)),max(LIMS(:,2))];
for ifig = 1:size(POIs,1)
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
pause(1)
Fignr = get(gcf,'number');
Fignr = Fignr+1;
%% generate plots for dipole with variable distance to sphere + variable orientation to vibration
% 3 plots one with maximal amplitude, one with optimal azimuth and one with
% optimal co-elevation
%vibration still in z-direction
VR = [];
VRfluc = [];
maxVRfluc = [];
idxOscDip = 1;
PHI=[];
THETA=[];
Pflag = 0;
% POIs

Theta = [0:pi/100:pi];
dPhi = pi/100;
[phi,theta]=meshgrid(0:dPhi:pi,Theta);
phisc = phi(:); thetasc = theta(:);
POIs=[];
dR = RSphere/resolution;
dAngle = pi/(10*resolution);
Rvals = 0:dR:RSphere-0.005;
Angles = 0:dAngle:pi;
for iR = 1:length(Rvals)
for iorien = 1:length(Angles)
% Select Sources and sinks + select vibrator
OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
CSource = LocDipole+d/2*OrienDipole;
CSink = LocDipole-d/2*OrienDipole;
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);
% point of interests
MeanDirection = ([0,0,1]+OrienDipole)/2;
if Angles(iorien)<pi/2
    OptimalPOI = MeanDirection;
else
    OptimalPOI = ([0,0,-1]+OrienDipole)/2;
end
POIs = RSphere.*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];

% Run simulation
Tsim = 0:dt:Tend;

for i=1:length(Tsim)
    t = Tsim(i);
    CSource(idxOscDip,:) = CSourceus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    CSink(idxOscDip,:) = CSinkus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    [VRMATSphere_tot,VRMATPOI_tot]=PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,1,[],'POI',POIs);
    VR(:,i) = VRMATPOI_tot;
    
end
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
PComp = round((iorien+(iR-1)*length(Angles))/(length(Rvals)*length(Angles))*100);
if PComp>=Pflag
disp([num2str(Pflag),'%'])
Pflag=Pflag+1;
end
end
end
%plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,10,16],VR)


% figure(Fignr+1)
% for i=1:size(VR,1)
%      plot3(i*ones(size(Tsim)),Tsim*1e6,(VR(i,:)-mean(VR(i,:)))*10^6)
%      hold on
%      VRfluc(i)=(max(VR(i,:))-min(VR(i,:)))*10^6;
% end
%  hold off
%  xlabel('POI')
%  ylabel('Time [µs]')
%  zlabel('Voltage [pV]')
%  figure(Fignr+2)
%  plot(1:length(VRfluc),VRfluc)
%  xlabel('POI')
%  ylabel('Potential variation [pV]')
[X,Y] = meshgrid(Angles,Rvals);
figure(Fignr+1)
sp={};
c={};
figure(Fignr+1)
subplot(1,3,1);
surf(X,Y,maxVRfluc,'EdgeColor','none','FaceColor','interp')
xlabel('Clockwise rotation in yz plane')
ylabel('Distance dipole form centre [cm]')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',0:RSphere/10:RSphere) 
set(gca,'YTickLabel',arrayfun(@num2str,get(gca,'YTick')*100,'UniformOutput',false))
colormap('jet')
xlim([0,pi]);
ylim([0,RSphere]);
c = colorbar;
c.Label.String = 'Fluctuation Amplitude [pV]';
set(gca,'colorscale','log')
c.Ticks = [min(c.Limits),c.Ticks,max(c.Limits)];
c.TickLabels = arrayfun(@num2str,c.Ticks,'UniformOutput',false);
title('Maximal Amplitude')
view(0,90)
subplot(1,3,2)
surf(X,Y,PHI,'EdgeColor','none','FaceColor','interp')
xlabel('Clockwise rotation in yz plane')
ylabel('Distance dipole form centre [cm]')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',0:RSphere/10:RSphere) 
set(gca,'YTickLabel',arrayfun(@num2str,get(gca,'YTick')*100,'UniformOutput',false))
colormap('jet')
xlim([0,pi]);
ylim([0,RSphere]);
c=colorbar;
title('Azimuth optimal POI')
c.Ticks = [0:pi/2:pi];
c.TickLabels = {'0','\pi/2','\pi'};
view(0,90)
subplot(1,3,3)
surf(X,Y,THETA,'EdgeColor','none')
xlabel('Clockwise rotation in yz plane')
ylabel('Distance dipole form centre [cm]')
zlabel('Fluctuation amplitude [pV]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',0:RSphere/10:RSphere) 
set(gca,'YTickLabel',arrayfun(@num2str,get(gca,'YTick')*100,'UniformOutput',false))
colormap('jet')
xlim([0,pi]);
ylim([0,RSphere]);
c = colorbar;
c.Ticks = [0:pi/4:pi];
c.TickLabels = {'0','\pi/4','\pi/2','3\pi/4','\pi'};
title('Co-elevation optimal POI')
view(0,90)

Fignr = get(gcf,'number');
Fignr = Fignr+1;
%%
function plotinfo(Tsim,Dirus,Aus,wus,Phaseus,POIs,VR)


subplot(4,2,2);
plot([Tsim*1e6.*ones(3,1)]',[Dirus'.*Aus.*sin(wus.*Tsim+Phaseus)]'*1e9)
legend({'x component','y component','z component'})
ylabel('displacement [nm]')
xlabel('Time [µs]')



ax1 = subplot(4,2,4);
nr1 = POIs(1);
plot(Tsim*1e6,(VR(nr1,:)-mean(VR(nr1,:)))*1e6)
xlabel('Time [µs]')
ylabel('Potential fluctuation [pV]')
title(['POI = ',num2str(nr1)])


ax2 = subplot(4,2,6);
nr2 = POIs(2);
plot(Tsim*1e6,(VR(nr2,:)-mean(VR(nr2,:)))*1e6)
xlabel('Time [µs]')
ylabel('Potential fluctuation [pV]')
title(['POI = ',num2str(nr2)])


ax3 = subplot(4,2,8);
nr3 = POIs(3);
plot(Tsim*1e6,(VR(nr3,:)-mean(VR(nr3,:)))*1e6)
xlabel('Time [µs]')
ylabel('Potential fluctuation [pV]')
title(['POI = ',num2str(nr3)])
LIMS = vertcat(get(ax1,'YLIM'),get(ax2,'YLIM'),get(ax3,'YLIM'));
YLIM = [min(LIMS(:,1)),max(LIMS(:,2))];
set(ax1,'YLIM',YLIM);set(ax2,'YLIM',YLIM);set(ax3,'YLIM',YLIM)

axes(ax1)
YTICKS = get(gca,'ytick');
yyaxis right
set(gca,'ylim',YLIM,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)+mean(VR(nr1,:)),12))))
ylabel('Absolute V [µV]')
axes(ax2)
YTICKS = get(gca,'ytick');
yyaxis right
set(gca,'ylim',YLIM,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)+mean(VR(nr2,:)),12))))
ylabel('Absolute V [µV]')
axes(ax3)
YTICKS = get(gca,'ytick');
yyaxis right
set(gca,'ylim',YLIM,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)+mean(VR(nr3,:)),12))))
ylabel('Absolute V [µV]')
end