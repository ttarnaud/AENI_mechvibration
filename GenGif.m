% Script for generation of Gif needed for presentation
close all
clear all
clc
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Inputs'))
if exist('Fignr','var')
    clearvars -except 'Fignr'
    Fignr = Fignr+1;
else
    clear all
    Fignr = 10;
end
CM = 'jet';
PLOT_DC = 0;
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
USwave = @(t) Dirus.*Aus.*sin(wus.*t+Phaseus);

% POIs
[xPOI,yPOI,zPOI] = sphere(100);


% Simulation settings
resUS = 8;
dt = (resUS*fus)^-1;
Tend = 2*fus^-1;
Tsim = 0:dt:Tend;
resolution = 3;
%%
%Rds = [0,0.01,0.06];
Rds = [0.06];
for iRd = 1:length(Rds)
Fignr = iRd;
VRMAT = [];
%Dipole location
LocDipole = [0,0,1];
Rdp = Rds(iRd);
CSource = (Rdp+d/2)*LocDipole;
CSink = (Rdp-d/2)*LocDipole;
% get settings
[Options,RSphere,RPOI,~,~,Angles] = getSettings('4SphereS8.7R~f',1);
% Create POIs
POIsScalp = (RSphere+0.012)*[xPOI(:),yPOI(:),zPOI(:)];
POIsAir = RPOI*[xPOI(:),yPOI(:),zPOI(:)];
% Variable vibration orientation
filename = ['StrengthfOrientation',num2str(Rdp),'.gif'];
for iorien=1:length(Angles)
    % Change wave
    Dirus = [0,sin(Angles(iorien)),cos(Angles(iorien))]; Dirus = Dirus/norm(Dirus);
    USwave = @(t) Dirus.*Aus.*sin(wus.*t+Phaseus);
    % Run simulations
    SettingsScalp = horzcat(Options,{'POI',POIsScalp,'resUS',resUS,'Display',double(iorien==0),'scale',0});
    SettingsAir = horzcat(Options,{'POI',POIsAir,'resUS',resUS,'Display',double(iorien==0),'scale',0});
    VRScalp = SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,1,SettingsScalp);
    VRAir = SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,1,SettingsAir);
    VRflucScalp = (max(VRScalp,[],2)-min(VRScalp,[],2))/2*10^6;
    VRflucAir = (max(VRAir,[],2)-min(VRAir,[],2))/2*10^6;
    if iorien == 1
        LIMS = [min(VRflucScalp(:)),max(VRflucScalp(:))];
    end
    figure(Fignr)
    subplot(1,2,1)
    PLOTSphere(POIsScalp,VRflucScalp,size(xPOI),CSource,CSink,Dirus,LIMS)
    ax = gca;
    set(ax.Colorbar,'visible','off')
    title('Signal at scalp')
    drawnow
    subplot(1,2,2)
    PLOTSphere(POIsAir,VRflucAir,size(xPOI),CSource,CSink,Dirus,LIMS)
    title('Signal 5 mm from scalp in air')
    drawnow
    %annotation('textbox', [0.9, 0.9, 0.1, 0.1], 'String', ['\theta = ',num2str(Angles(iorien)/pi),'\pi'])
    h = gcf;
    set(h,'position',[-1690,310,1344,649],'color','w')
    mtit( ['\theta = ',num2str(Angles(iorien)/pi),'\pi'],'fontsize',20)
    %LIMS=[-1,2];
    %ax = gca;
    %set(ax,'CLIM',LIMS);
    %set(ax.Colorbar,'Limits',LIMS);
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if iorien == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
    close(h)
    disp([num2str(iorien/length(Angles)*100),'%']);
end
end
%%
function PLOTSphere(POIspos,vals,reshapeSize,CSource,CSink,vibOrien,LIMS)
Scale = 5;
Isize = 6;
% Plot arrow dipole
d = CSource-CSink;
cdipole = 1/2*(CSource+CSink);
Arrowpoints = [(cdipole-d/norm(d)*0.07/Scale)',(cdipole+d/norm(d)*0.07/Scale)'];

ARROW = arrow3d(Arrowpoints(1,:),Arrowpoints(2,:),Arrowpoints(3,:),0.75,2*Isize*0.07/400,2*Isize*3*0.07/400);
set(ARROW,'facecolor',[0,0,0])
hold on
%plot arrow vibration
d = vibOrien;
cdipole = 1/2*(CSource+CSink);
Arrowpoints = [(cdipole-d/norm(d)*0.07/Scale)',(cdipole+d/norm(d)*0.07/Scale)'];
ARROW = arrow3d(Arrowpoints(1,:),Arrowpoints(2,:),Arrowpoints(3,:),0.75,2*Isize*0.07/400,2*Isize*3*0.07/400);
set(ARROW,'facecolor','m')


x = reshape(POIspos(:,1),reshapeSize);
y = reshape(POIspos(:,2),reshapeSize);
z = reshape(POIspos(:,3),reshapeSize);
vals = reshape(vals,reshapeSize);
surf(x,y,z,vals,'FaceAlpha',0.5,'EdgeColor','none','SpecularStrength',0);
view(120,30)
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
axis equal
hAxis=gca;
hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
c=colorbar;
c.Label.String = 'Potential at sphere boundary [pV]';
c.Label.FontSize = 20;
colormap('jet')
%LIMS = [floor(min(VRMATSphere_tot(:))),ceil(max(VRMATSphere_tot(:)))];
%LIMS = [-0.4,0.4];
%LIMS = [min(vals(:)),max(vals(:))];
if any(LIMS)
    set(hAxis,'CLIM',LIMS);
    set(c,'Limits',LIMS);
end
hold off
set(gca,'fontsize',20)
end