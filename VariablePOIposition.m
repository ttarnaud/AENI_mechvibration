% Run simulations for moving dipole
% Varying Radius position of Point of interest => see signal decay 
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
resUS = 10;
dt = (1e2*fus)^-1;
Tend = 2*fus^-1;
Tsim = 0:dt:Tend;

% Geometric  properties
tAir = 0.005;
RSphere = 0.07;
RPOI = RSphere+0.012+tAir;
SolutionType = '4Sphere';
RatioSkT = 1/25;
Rdp = [0.05,0.060];
fDependence = 1;
SphereRes = 100;
figure
for iRdp = 1:length(Rdp)
OrienDipole = [0,0,1];
LocDipole = Rdp(iRdp).*[0,0,1];
CSource = LocDipole+d/2*OrienDipole;
CSink = LocDipole-d/2*OrienDipole;

POIs = [];
RPOIs = sort([linspace(0.065,RPOI,101),RSphere,RSphere+0.005,RSphere+0.012]);
%RPOIs(1) = RPOIs(1)+1.1*d;
POIs = [0,0,1].*RPOIs';

VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,1,'RatioSkT',RatioSkT,'POI',POIs,...
    'resUS',resUS,'RSphere',RSphere,'SolutionType',SolutionType,'fDependence',fDependence,'tAir',tAir,'scale',0);

VRfluc(:,iRdp) = (max(VR,[],2)-min(VR,[],2))/2*10^6;
VRmean(:,iRdp) = (max(VR,[],2)+min(VR,[],2))/2;
subplot(3,2,1)
hold on
plot(RPOIs,VRmean(:,iRdp),'DisplayName',['dipole @ ',num2str(Rdp(iRdp))])
xlabel('positions POI')
ylabel('V_{DC} [µV]')
hold off
subplot(3,2,2)
hold on
plot(RPOIs,VRmean(:,iRdp),'DisplayName',['dipole @ ',num2str(Rdp(iRdp))])
xlabel('positions POI')
ylabel('V_{DC} [µV]')
legend('show')
xlim([0.08,RPOIs(end)])
hold off
subplot(3,2,3)
hold on
plot(RPOIs,VRfluc(:,iRdp),'DisplayName',['dipole @ ',num2str(Rdp(iRdp))])
xlabel('positions POI')
ylabel('V_{1MHz} [pV]')
hold off
subplot(3,2,4)
hold on
plot(RPOIs,VRfluc(:,iRdp),'DisplayName',['dipole @ ',num2str(Rdp(iRdp))])
xlabel('positions POI')
ylabel('V_{1MHz} [pV]')
legend('show')
xlim([0.08,RPOIs(end)])
hold off



end
for iRdp = 1:2
subplot(3,2,5)
hold on
plot(RPOIs',VRmean(:,iRdp)./(max(VRmean(:))),'DisplayName',['DC - dipole @ ',num2str(Rdp(iRdp))])
plot(RPOIs',VRfluc(:,iRdp)./(max(VRfluc(:))),'DisplayName',['1MHz - dipole @ ',num2str(Rdp(iRdp))])
xlabel('positions POI')
ylabel('normalized V_{DC} [-]')
hold off
subplot(3,2,6)
hold on
plot(RPOIs',VRmean(:,iRdp)./(max(VRmean(:))),'DisplayName',['DC dipole @ ',num2str(Rdp(iRdp))])
plot(RPOIs',VRfluc(:,iRdp)./(max(VRfluc(:))),'DisplayName',['1MHz - dipole @ ',num2str(Rdp(iRdp))])
xlabel('positions POI')
ylabel('normalized V_{DC} [-]')
legend('show')
xlim([0.08,RPOIs(end)])
hold off
end