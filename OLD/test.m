% test results of SimDipoleOSC.m

% determine fluctution effect based on elevation and azimuth
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

% Simulation settings
resUS = 100;
dt = (resUS*fus)^-1;
Tend = 2*1e-6;

% Geometric  properties
RSphere = 0.07;  %m
[xs,ys,zs] = sphere(1000);
xsphere = xs*RSphere; ysphere = ys*RSphere; zsphere = zs*RSphere;
Inputtype = '';
SolutionType = '3Sphere';

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
tic
OrienDipole = [0,1/sqrt(2),1/sqrt(2)];

% figure(Fignr)
% ax0 = subplot(4,2,[1,3,5,7]);
%PotentialSphere_Multi(d/2*LocDipoleZ,-d/2*LocDipoleZ,I,sigma,RSphere,1000,Fignr);

%Select Sources and sinks + select vibrator
CSource = vertcat((Rd+d/2)/Rd*LocDipoleRd,(0.01+d/2)*OrienDipole);
CSink = vertcat((Rd-d/2)/Rd*LocDipoleRd,(0.01-d/2)*OrienDipole);
idxOscDip = size(CSource,1);


% with new function
tic
USwave = @(t) Dirus.*Aus.*sin(wus.*t+Phaseus);
VRPOIs=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,'SolutionType',SolutionType,...
    'resUS',resUS,'display',1,'ShowSphere','noPOI');
VRfluc2 = (max(VRPOIs,[],2)-min(VRPOIs,[],2))*10^6;

VRflucSurf = reshape(VRfluc2,size(phi));
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
toc
%% with new function & PC  
% Usefull for when using 3Sphere model
tic
USwave = @(t) Dirus.*Aus.*sin(wus.*t+Phaseus);
VRPOIsPC=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,'SolutionType',SolutionType,'resUS',resUS,...
    'ParallelCompute',1,'display',1,'ShowSphere','noPOI');
VRfluc2PC = (max(VRPOIsPC,[],2)-min(VRPOIsPC,[],2))*10^6;

VRflucSurf = reshape(VRfluc2PC,size(phi));
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
toc
%% OLD METHOD
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);
% Run simulation
Tsim = 0:dt:Tend;

for i=1:length(Tsim)
    t = Tsim(i);
    CSource(idxOscDip,:) = CSourceus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    CSink(idxOscDip,:) = CSinkus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    [VRMATSphere_tot,VRMATPOI_tot]=PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,1,[],'POI',POIs,'SolutionType',SolutionType);
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
toc