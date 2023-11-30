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

% Simulation settings
resUS = 100;
dt = (resUS*fus)^-1;
Tend = 2*fus^-1;
Tsim = 0:dt:Tend;
resolution = 3;

% Geometric  properties
SolutionType = '3SphereS8.2R~f';
PlotType = '3Sphere';
[Options,RSphere] = getSettings(SolutionType,resolution);
SphereRes = 100;


% dipoles on a sphere
Rd = 0.03; % concentric sphere on which dipoles are situated
dTheta = pi/3;
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
Theta = [0:pi/4:pi/2];
dPhi = pi/4;
idx = 0;
POIs=[];
for i=Theta
    if i==0 || i==pi
        Angles = 0;
    else
        Angles = [dPhi:dPhi:2*pi];
    end
    for j= 1:length(Angles)
        idx=idx+1;
        POIs(idx,:) = RSphere*[cos(Angles(j))*sin(i),sin(Angles(j))*sin(i),cos(i)];
    end
end

%% dipole in center and z-direction
LocDipoleZ = [0,0,1];
idxOscDip = 1;
CSource = [0,1,1]; CSource = CSource/norm(CSource)*0.065*(1+d/2);
CSink = [0,1,1]; CSink = CSink/norm(CSink)*0.065*(1-d/2);

figure(Fignr)
ax0 = subplot(4,2,[1,3,5,7]);
VRMATSphere_tot = PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,SphereRes,Fignr,'DOI',idxOscDip,'POI',POIs,'SolutionType',PlotType);
title(['Bounded Medium: max = ',num2str(max(VRMATSphere_tot(:)))])

% Run simulation
Settings = horzcat(Options,{'POI',POIs,'resUS',resUS});
VR =  SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,Settings);
if idxOscDip>size(CSource,1)
    plotinfo(Tsim,Dirus,0,wus,Phaseus,[1,8,16],VR,0)
else
    plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,8,16],VR,0)
end
mtit(SolutionType)
figure(Fignr+1)
for i=1:size(VR,1)
     plot3(i*ones(size(Tsim)),Tsim*1e6,(VR(i,:)-mean(VR(i,:)))*10^6)
     hold on
     VRfluc(i)=(max(VR(i,:))-min(VR(i,:)))/2*10^6;
end
 hold off
 xlabel('POI')
 ylabel('Time [µs]')
 zlabel('Voltage [pV]')

 figure(Fignr+2)
 plot(1:length(VRfluc),VRfluc)
 xlabel('POI')
 ylabel('Potential variation [pV]')
 Fignr = get(gcf,'number');
Fignr = Fignr+1;



Fignr = get(gcf,'number');
Fignr = Fignr+1;
%% dipole in center and y-direction
LocDipoleZ = [0,1,0];
idxOscDip = 1;
figure(Fignr)
ax0 = subplot(4,2,[1,3,5,7]);
PotentialSphere_Multi(d/2*LocDipoleZ,-d/2*LocDipoleZ,I,sigma,RSphere,SphereRes,Fignr,'DOI',idxOscDip,'POI',POIs,'SolutionType',PlotType);

% Select Sources and sinks + select vibrator
CSource = d/2*LocDipoleZ;
CSink = -d/2*LocDipoleZ;
% Run simulation
Settings = horzcat(Options,{'POI',POIs,'resUS',resUS});
VR =  SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,Settings);

plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,8,16],VR,0)

Fignr = get(gcf,'number');
Fignr = Fignr+1;
%% multiple dipoles + vibrator in center
% Select Sources and sinks + select vibrator
CSource = vertcat((Rd+d/2)/Rd*LocDipoleRd,[0,0,d/2]);
CSink = vertcat((Rd-d/2)/Rd*LocDipoleRd,[0,0,-d/2]);
idxOscDip = size(CSource,1);

figure(Fignr)
ax0 = subplot(4,2,[1,3,5,7]);
PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,SphereRes,Fignr,'POI',POIs,'DOI',idxOscDip,'SolutionType',PlotType);
Settings = horzcat(Options,{'POI',POIs,'resUS',resUS});
VR =  SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,Settings);


plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,8,16],VR,0)
Fignr = get(gcf,'number');
Fignr = Fignr+1;
%% multiple dipoles + vibrator in center

% Select Sources and sinks + select vibrator
CSource = vertcat((Rd+d/2)/Rd*LocDipoleRd,[0,d/2,0]);
CSink = vertcat((Rd-d/2)/Rd*LocDipoleRd,[0,-d/2,0]);
idxOscDip = size(CSource,1);
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);

figure(Fignr)
ax0 = subplot(4,2,[1,3,5,7]);
PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,SphereRes,Fignr,'POI',POIs,'DOI',idxOscDip,'SolutionType',PlotType);
Settings = horzcat(Options,{'POI',POIs,'resUS',resUS});
VR =  SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,Settings);
plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,8,16],VR,0)
Fignr = get(gcf,'number');
Fignr = Fignr+1;
%% multiple dipoles + vibrator in at the side

% Select Sources and sinks + select vibrator
CSource = vertcat((Rd+d/2)/Rd*LocDipoleRd,[0,0,d/2]);
CSink = vertcat((Rd-d/2)/Rd*LocDipoleRd,[0,0,-d/2]);
idxOscDip = size(CSource,1)-5;
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);

figure(Fignr)
ax0 = subplot(4,2,[1,3,5,7]);
PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,SphereRes,Fignr,'POI',POIs,'DOI',idxOscDip,'SolutionType',PlotType);
Settings = horzcat(Options,{'POI',POIs,'resUS',resUS});
VR =  SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,Settings);
plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,8,16],VR,0)
Fignr = get(gcf,'number');
Fignr = Fignr+1;

%% 
function plotinfo(Tsim,Dirus,Aus,wus,Phaseus,POIs,VR,YLIM_flag)


subplot(4,2,2);
plot([Tsim*1e6.*ones(3,1)]',[Dirus'.*Aus.*sin(wus.*Tsim+Phaseus)]'*1e9)
legend({'x component','y component','z component'})
ylabel('displacement [nm]')
xlabel('Time [µs]')



ax1 = subplot(4,2,4);
nr1 = POIs(1);
plot(Tsim*1e6,round((VR(nr1,:)-mean(VR(nr1,:)))*1e6,7))
xlabel('Time [µs]')
ylabel('Potential fluctuation [pV]')
title(['POI = ',num2str(nr1)])


ax2 = subplot(4,2,6);
nr2 = POIs(2);
plot(Tsim*1e6,round((VR(nr2,:)-mean(VR(nr2,:)))*1e6,7))
xlabel('Time [µs]')
ylabel('Potential fluctuation [pV]')
title(['POI = ',num2str(nr2)])


ax3 = subplot(4,2,8);
nr3 = POIs(3);
plot(Tsim*1e6,round((VR(nr3,:)-mean(VR(nr3,:)))*1e6,7))
xlabel('Time [µs]')
ylabel('Potential fluctuation [pV]')
title(['POI = ',num2str(nr3)])
LIMS = vertcat(get(ax1,'YLIM'),get(ax2,'YLIM'),get(ax3,'YLIM'));
YLIM = [min(LIMS(:,1)),max(LIMS(:,2))];
if YLIM_flag
set(ax1,'YLIM',YLIM);set(ax2,'YLIM',YLIM);set(ax3,'YLIM',YLIM)
axes(ax1)
YTICKS = get(gca,'ytick');
yyaxis right
set(gca,'ylim',YLIM,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)*1e-6+mean(VR(nr1,:)),12))))
ylabel('Absolute V [µV]')
axes(ax2)
YTICKS = get(gca,'ytick');
yyaxis right
set(gca,'ylim',YLIM,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)*1e-6+mean(VR(nr2,:)),12))))
ylabel('Absolute V [µV]')
axes(ax3)
YTICKS = get(gca,'ytick');
yyaxis right
set(gca,'ylim',YLIM,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)*1e-6+mean(VR(nr3,:)),12))))
ylabel('Absolute V [µV]')
else
axes(ax1)
YTICKS = get(gca,'ytick');
YLIMS = get(gca,'ylim');
yyaxis right
set(gca,'ylim',YLIMS,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)*1e-6+mean(VR(nr1,:)),12))))
ylabel('Absolute V [µV]')
axes(ax2)
YTICKS = get(gca,'ytick');
YLIMS = get(gca,'ylim');
yyaxis right
set(gca,'ylim',YLIMS,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)*1e-6+mean(VR(nr2,:)),12))))
ylabel('Absolute V [µV]')
axes(ax3)
YTICKS = get(gca,'ytick');
YLIMS = get(gca,'ylim');
yyaxis right
set(gca,'ylim',YLIMS,'ytick',YTICKS,'yticklabel',cellstr(num2str(round(YTICKS(:)*1e-6+mean(VR(nr3,:)),12))))
ylabel('Absolute V [µV]')
end
end