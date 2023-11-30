% Figure 1 results paper
% code copied form Flucdistribution and RunSimulation US
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
clear all
close all
clc
Fignr = 1;
CM = inferno;
GENGIF = false;
resolution = 30;

% Electrical properties
Npercolumn = 1e4;  %Npercolumn = 1e2; for minicolumn
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
USwave = @(t,idp) Dirus.*Aus.*sin(wus.*t+Phaseus);

% Simulation settings
resUS = 20;
dt = (resUS*fus)^-1;
Tend = 2*fus^-1;
Tsim = 0:dt:Tend;

RPOIs= [];
% Geometric  properties
RSphere = 0.07;
RPOI = RSphere+0.017;
SolutionType = '4Sphere';
RatioSkT = 1/25;
Rvals = linspace(0,RSphere-0.005,resolution+1);
RTICKS = unique([0:0.01:Rvals(end),Rvals(end),fliplr(RPOI:-0.01:Rvals(end))],'stable');
fDependence = 0;
scale_flag = 0;
RPOIs = [0,0.012,0.017]+RSphere;
    
% POIs    
OrienDipole = [0,1,2]; OrienDipole = OrienDipole/norm(OrienDipole);
POIs = [1,0,0;0,0,1;0,1,0;OrienDipole;1,1,1];
POIs = POIs./vecnorm(POIs,2,2);
POIs = RPOI.*POIs;
[xSphere,ySphere,zSphere] = sphere(40);
POIs_intm = [POIs;RPOI*[xSphere(:),ySphere(:),zSphere(:)]];
tic

idxOscDip = 1;
%% collect data first second subplots
% Select Sources and sinks + select vibrator
RDOI = 0.02;
CSource = (RDOI+d/2)*OrienDipole;
CSink = (RDOI-d/2)*OrienDipole;
% Run simulation
Tsim = 0:dt:Tend;

VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs_intm,'SolutionType',SolutionType,'resUS',resUS,'scale',scale_flag,'fDependence', fDependence);
VRPOI = VR(1:length(POIs),:);
VRSphere = VR(length(POIs)+1:end,:);
VRfluc = (max(VRSphere,[],2)-min(VRSphere,[],2))*10^6/2; %uV -> pV
%% create first subplot
% create sphere with color the oscillation amplitude
figure(Fignr)
CM = inferno;
% draw dipole
dp = CSource-CSink;
cdipole = 1/2*(CSource+CSink);
Arrowpoints = [(cdipole-dp/norm(dp)*RSphere/15)',(cdipole+dp/norm(dp)*RSphere/15)'];
Iarrow = I*5;
ARROW = arrow3d(Arrowpoints(1,:),Arrowpoints(2,:),Arrowpoints(3,:),0.75,2*(Iarrow/10+1)*RSphere/400,(Iarrow/10+1)*6*RSphere/400);
set(ARROW,'facecolor',[1,0,0])
hold on
surf(RPOI*xSphere,RPOI*ySphere,RPOI*zSphere,reshape(VRfluc,size(xSphere)),'FaceAlpha',0.2,'EdgeColor','none','SpecularStrength',0);
hold off
view(110,30)
xlabel('xaxis')
ylabel('yaxis')
zlabel('zaxis')
axis equal
hAxis=gca;
hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
c=colorbar;
c.Label.String = 'Signal Strength @ 1 MHz [pV]';
c.Label.FontSize = 12;
cmap = colormap(CM);
%LIMS = [floor(min(VRMATSphere_tot(:))),ceil(max(VRMATSphere_tot(:)))];
%LIMS = [-0.4,0.4];
LIMS = [min(VRfluc(:)),max(VRfluc(:))];
set(hAxis,'CLIM',LIMS);
set(c,{'Limits','fontsize'},{LIMS,10});
c.Label.FontSize = 12;

DrawElectrodes(POIs,'radius',0.01)
set(gca,{'color','box'},{'None','off'})
grid off
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})

set(gca,'position',[-0.0083    0.0659    0.7472    0.9066])
set(c,'position',[0.6987    0.0789    0.0489    0.8542])
xl = get(gca,'xlabel'); set(xl,'position', [0.0297   -0.0488   -0.1056]);
yl = get(gca,'ylabel'); set(yl,'position', [-0.0108    0.0549   -0.1017]);
zl = get(gca,'zlabel'); set(zl,'position',  [0.0241   -0.0946    0.0174]);
numbers = findall(gca,'type','text');
set(numbers(1),'position',[0.0563 0.0485 0.0563])
set(numbers(2),'position',[-0.0086 0.0566 0.0906])
set(numbers(3),'position',[-0.0023 0.0942 0.0042])
set(numbers(4),'position',[-0.0088 -0.0067 0.1133])
set(numbers(5),'position',[0.0933 -0.0024 0.0042])

%% second subplot
figure()
ax1 = subplot(3,1,1);
yvals = [Dirus'.*Aus.*sin(wus.*Tsim+Phaseus)]'*1e9;
xvals = [Tsim*1e6.*ones(3,1)]';
ls = {'-','-','-'};
lnames = {'x', 'y', 'z'};
for i=1:3
plot(xvals(:,i),yvals(:,i),'color',((i-1)*0.3+0.1)*ones(1,3),'linestyle',ls{i}, 'linewidth', 1+double(i==1),'HandleVisibility','off')
hold on
plot(nan,nan,'color',((i-1)*0.2+0.1)*ones(1,3), 'linewidth', 1, 'DisplayName', lnames{i})

end
hold off
l1 = legend('show','box','off','location','eastoutside');
ylabel('\Delta{\bf r} [nm]')
%xlabel('time [\mus]')
set(gca,'box','off')

cmap = turbo(11);
cmap = cmap(7:end,:);
cmap = flare;
ax2 = subplot(3,1,[2,3]);
for nr1 = 1:5
yvals = round((VRPOI(nr1,:)-mean(VRPOI(nr1,:)))*1e6,7);
amp = max(yvals)-min(yvals)/2;
ci = round(min(max((length(cmap)-1)/(diff(LIMS))*(amp-LIMS(1))+1,1),length(cmap)));
hold on
plot(Tsim*1e6,yvals,'color',cmap(ci,:),'linewidth',1,'DisplayName', sprintf('POI_{%s}',num2str(nr1)))
end
hold off
xlabel('time [\mus]')
ylabel(' [pV]')
l2 = legend('show','box','off','NumColumns',1,'location','eastoutside');
toc
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
set(ax1,'position',[0.1639    0.8206    0.5603    0.1044])
set(ax2,'position',[0.1639    0.1821    0.5603    0.4433])
set(l2,'position',[0.7406    0.1620    0.2492    0.4261])
set(l1,'position',[0.7406    0.7642    0.1807    0.1972])
set(get(ax1,'ylabel'),'position',[-0.3346   -3.0769   -1.0000])
set(get(ax2,'ylabel'),'position',[-0.3346   -0.0018   -1.0000])
set(get(ax2,'xlabel'),'position',[1.0126   -0.1311   -1.0000])
%% collect data supblot 3
Tend = 1*fus^-1;
Tsim = 0:dt:Tend;

Theta = [0:pi/30:pi];
dPhi = pi/10;
[phi,theta]=meshgrid(0:dPhi:pi,Theta);
phisc = phi(:); thetasc = theta(:);
POIs=[];
dAngle = pi/(10*resolution);
Angles = 0:dAngle:pi; %Angles = Angles(1:3:end);
iscale_end=1;
Show_sphere={'',''};

VR = [];
VRfluc = [];
maxVRfluc = nan(length(Rvals),length(Angles));
idxOscDip = 1;
PHI = nan(length(Rvals),length(Angles));
THETA = nan(length(Rvals),length(Angles));
Pflag = 0;

tic
for iR = 1:length(Rvals)
    parfor iorien = 1:length(Angles)
        % Select Sources and sinks + select vibrator
        OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
        LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
        CSource = LocDipole+d/2*OrienDipole;
        CSink = LocDipole-d/2*OrienDipole; 
        
        POIs = RPOI.*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];
 
        % Run simulation
        Tsim = 0:dt:Tend;
        VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,...
            'resUS',resUS,'Display',double(iorien==1&&iR==1),'SolutionType',SolutionType,'scale',0,'ShowSphere',Show_sphere{double(iR==1&&iorien==1)+1},'fDependence', fDependence);
        VRfluc = (max(VR,[],2)-min(VR,[],2))/2*10^6;
        
        [maxVRfluc(iR,iorien),idxmaxVRfluc]=max(VRfluc);
        if POIs(idxmaxVRfluc,1) <= 1e-3 && (POIs(idxmaxVRfluc,2)~=0||POIs(idxmaxVRfluc,3)~=0)
            PHI(iR,iorien) = pi/2;
        elseif POIs(idxmaxVRfluc,1) <= 1e-3 && POIs(idxmaxVRfluc,2) <= 1e-3
            PHI(iR,iorien) = pi/2;
        else
            PHI(iR,iorien) = atan(POIs(idxmaxVRfluc,2)/POIs(idxmaxVRfluc,1))+...
                double(POIs(idxmaxVRfluc,1)<0)*pi+double(POIs(idxmaxVRfluc,1)>0&&POIs(idxmaxVRfluc,2)<0)*2*pi; %map values [-pi/2,pi/2] -> [0,2pi]
        end
%         if Angles(iorien) == 0 || Angles(iorien) == pi
%             THETA(iR,iorien) = Angles(iorien);
%         else
%             THETA(iR,iorien) = acos(POIs(idxmaxVRfluc,3)/RPOI_sel);
%         end
        THETA(iR,iorien) = acos(POIs(idxmaxVRfluc,3)/RPOI);
        if iR==1
            if abs(THETA(iR,iorien)-0)<1e-6 && Angles(iorien)>pi/2
                THETA(iR,iorien)= pi;
            elseif abs(THETA(iR,iorien)-pi)<1e-6 && Angles(iorien)<pi/2
                THETA(iR,iorien)= 0;
            end
        end

    end
    PComp = round(iR/length(Rvals)*100);
    if PComp>=Pflag
        disp([num2str(PComp),'%'])
        Pflag=Pflag+1;
    end
end
toc
%% create figure subplot 3

[X,Y] = meshgrid(Angles,Rvals*1000);
Z = abs(THETA-X);
figure()
surf(X,Y,Z,'EdgeColor','none')
xlabel('\theta_{dp}')
ylabel('R [mm]')
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',[0,30,65,87])
colormap(rdpu)
xlim([0,pi]);
ylim([0,RPOI*1000]);
c = colorbar;
c.Ticks = [0:pi/8:pi/4];
c.TickLabels = {'0','\pi/8','\pi/4','3\pi/4','\pi'};
c.Label.String = '|\theta_{POI}-\theta_{dp}|';
set(c,{'fontsize'},{10});
c.Label.FontSize = 12;
view(0,90)
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
c.Label.Position = [2.4889 0.3868 0];
set(gca,'xGrid','off')
%% collect data subplot 4,5,6
Tend = 1*fus^-1;
Tsim = 0:dt:Tend;

Theta = [0:pi/300:pi];
dPhi = pi/4;
[phi,theta]=meshgrid(0:dPhi:pi,Theta);
phisc = phi(:); thetasc = theta(:);
POIs=[];
dAngle = pi/(10*resolution);
Angles = 0:dAngle:pi; 


Rvals = linspace(0,RSphere-0.005,resolution+1);



VR = [];
VRfluc = [];
maxVRfluc = nan(length(Rvals),length(Angles),length(RPOIs));
idxOscDip = 1;
Pflag = 0;
tic
for iscale=1:length(RPOIs)

    RPOI_sel = RPOIs(iscale);
    titles = {'POIs on Brain','POIs on Scalp','POIs 5mm from scalp'};

for iR = 1:length(Rvals)
    for iorien = 1:length(Angles)
        % Select Sources and sinks + select vibrator
        OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
        LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
        CSource = LocDipole+d/2*OrienDipole;
        CSink = LocDipole-d/2*OrienDipole; 
        
        POIs = RPOI_sel.*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];

        % Run simulation
        Tsim = 0:dt:Tend;
        VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,...
            'resUS',resUS,'Display',double(iorien==1&&iR==1),'SolutionType',SolutionType,'scale',0,'ShowSphere','','fDependence', fDependence);
        VRfluc = (max(VR,[],2)-min(VR,[],2))/2*10^6;
        [maxVRfluc(iR,iorien,iscale),idxmaxVRfluc]=max(VRfluc);

    end
    PComp = round(iR/length(Rvals)*100);
    if PComp>=Pflag
        disp([num2str(PComp),'%'])
        Pflag=Pflag+1;
    end
end
end
toc
%save('data_sp456v2.mat')
%%
signalStrenghtattheta0r0 = squeeze(maxVRfluc(1,1,:));
deltaRs = linspace(22,70,5)/1000;
[~,idxpi2] = min(abs(Angles-pi/2));
myinfo = {[],[]};
for i=1:length(RPOIs)
    for x = 1:length(deltaRs)
        R0 = RPOIs(i)-deltaRs(x)
        [~,idxR] = min(abs(Rvals-R0));
        myinfo{1}(i,x) = maxVRfluc(idxR,1,i);
        myinfo{2}(i,x) = maxVRfluc(idxR,idxpi2,i);
    end
end

for myi=1:2
    1-myinfo{myi}([2,3],:)./myinfo{myi}([1,2],:)
end
%% 

CM = rocket;
[X,Y] = meshgrid(Angles,Rvals*1000);
sp={};
c={};
for i = 1:length(RPOIs)
figure()
Z = maxVRfluc(:,:,i);
surf(X,Y,Z,'EdgeColor','none')
% hold on
% cont = contourc(unique(X(:)),unique(Y(:)),Z,[1,1]);
% plotContour(cont,[1,0.3,0.3],0)
% hold off
% hold on
% cont = contourc(unique(X(:)),unique(Y(:)),log10(Z),3);
% plotContour(cont,[0.3,0.3,0.3],1)
% hold off
xlabel('\theta_{dp}')
ylabel('R [mm]')

set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',[0,30,65,RPOIs(i)*1000])
colormap(CM)
xlim([0,pi]);
ylim([0,RPOIs(end)*1000]);
c = colorbar;
c.Label.String = '|\Psi(f=1MHz)| [pV]';
set(gca,'colorscale','log')
c.Ticks = [min(Z(:)),c.Ticks,max(Z(:))];
c.TickLabels = arrayfun(@num2str,round(c.Ticks,3),'UniformOutput',false);
set(c,{'fontsize'},{10});
c.Label.FontSize = 12;
view(0,90)
set(gca,'box','off')
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
%c.Label.Position = double(i==1)*[2.1    1.4986         0]+double(i==2)*[2.4222 0.2412 0]+double(i==3)*[2.1, mean([c.Ticks(1),c.Ticks(end)]),0];
c.Label.Position = [2.3, 10^mean([log10(c.Ticks(1)),log10(c.Ticks(end))]),0];
set(gca,'xGrid','off')

end
%% only contour
CM = rocket;
[X,Y] = meshgrid(Angles,Rvals*1000);
sp={};
c={};
for i = 1:length(RPOIs)
figure()
Z = maxVRfluc(:,:,i);
levels = round(logspace(log10(min(Z(:))),log10(max(Z(:))),7),3,'significant')

contour(X,Y,Z,levels,'showText','off');


xlabel('\theta_{dp}')
ylabel('r_{dp} [mm]')

set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',[0,30,65,RPOIs(i)*1000])
colormap(CM)
xlim([0,pi]);
ylim([0,RPOIs(end)*1000]);
c = colorbar;
c.Label.String = '|\Psi(f=1MHz)| [pV]';
set(gca,'clim',[min(Z(:)),max(Z(:))])
set(gca,'colorscale','log')
c.Ticks = [min(Z(:)),c.Ticks,max(Z(:))];
c.TickLabels = arrayfun(@num2str,round(c.Ticks,3),'UniformOutput',false);
set(c,{'fontsize'},{9});
c.Label.FontSize = 10;
view(0,90)
set(gca,'box','off')
set(findall(gcf,'type','axes'),'fontsize',9)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,6.34,5.21],'centimeters',[1+6.34,3+5.21],'Painters'})
%c.Label.Position = double(i==1)*[2.1    1.4986         0]+double(i==2)*[2.4222 0.2412 0]+double(i==3)*[2.1, mean([c.Ticks(1),c.Ticks(end)]),0];
c.Label.Position = [2.3, 10^mean([log10(c.Ticks(1)),log10(c.Ticks(end))]),0];
set(gca,{'xGrid','yGrid'},{'off','on'})

end
%% surf and contour

contclr = [0.9,0.9,0.9]
CM = rocket;
[X,Y] = meshgrid(Angles,Rvals*1000);
sp={};
c={};
for i = 1:length(RPOIs)
figure()
Z = maxVRfluc(:,:,i);
surf(X,Y,zeros(size(Z)),Z,'EdgeColor','none')
hold on
levels = round(logspace(log10(min(Z(:))),log10(max(Z(:))),7),2,'significant')
[C,h] = contour(X,Y,Z,levels,'showText','on',LineColor=contclr);
clabel(C,h,'fontsize', 9,'Color',contclr)
hold off

xlabel('\theta_{dp}')
ylabel('r_{dp} [mm]')

set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',[0,30,65,RPOIs(i)*1000])
colormap(CM)
xlim([0,pi]);
ylim([0,RPOIs(end)*1000]);
c = colorbar;
c.Label.String = '|\Psi(f=1MHz)| [pV]';
set(gca,'clim',[min(Z(:)),max(Z(:))])
set(gca,'colorscale','log')
c.Ticks = [min(Z(:)),c.Ticks,max(Z(:))];
c.TickLabels = arrayfun(@num2str,round(c.Ticks,3),'UniformOutput',false);
set(c,{'fontsize'},{9});
c.Label.FontSize = 10;
view(0,90)
set(gca,'box','off')
set(findall(gcf,'type','axes'),'fontsize',9)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,6.34,5.21],'centimeters',[1+6.34,3+5.21],'Painters'})
%c.Label.Position = double(i==1)*[2.1    1.4986         0]+double(i==2)*[2.4222 0.2412 0]+double(i==3)*[2.1, mean([c.Ticks(1),c.Ticks(end)]),0];
c.Label.Position = [2.3, 10^mean([log10(c.Ticks(1)),log10(c.Ticks(end))]),0];
set(gca,'xGrid','off')

end
%% Polar

contclr = [0.9,0.9,0.9]
CM = rocket;
[X,Y] = meshgrid(Angles,Rvals*1000);
sp={};
c={};
for i = 1:1%length(RPOIs)
figure()
Z = maxVRfluc(:,:,i);
X2 = Y.*cos(X);
Y2 = Y.*sin(X);
x2 = X(1,:);
y2 = Y(:,1);
polarplot3d(Z,'plottype','surfn','angularrange',x2','radialrange',y2,...
              'polargrid',{1 1},'tickspacing',30,'Colordata',Z,...
              'plotprops',{'Linestyle','none','EdgeColor','None'});
hold on
levels = round(logspace(log10(min(Z(:))),log10(max(Z(:))),7),2,'significant')
[C,h] = contour(X2,Y2,Z,levels,'showText','on',LineColor=contclr, Labelspacing=1000);
clabel(C,h,'fontsize', 9,'Color',contclr)
hold off
view(180,-90)

xlabel('\theta_{dp}')
ylabel('r_{dp} [mm]')


colormap(CM)

c = colorbar;
c.Label.String = '|\Psi(f=1MHz)| [pV]';
set(gca,'clim',[min(Z(:)),max(Z(:))])
set(gca,'colorscale','log')
c.Ticks = [min(Z(:)),c.Ticks,max(Z(:))];
c.TickLabels = arrayfun(@num2str,round(c.Ticks,3),'UniformOutput',false);
set(c,{'fontsize'},{9});
c.Label.FontSize = 10;

set(gca,'box','off')
set(findall(gcf,'type','axes'),'fontsize',9)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,6.34,6.34],'centimeters',[1+6.34,3+6.34],'Painters'})
%c.Label.Position = double(i==1)*[2.1    1.4986         0]+double(i==2)*[2.4222 0.2412 0]+double(i==3)*[2.1, mean([c.Ticks(1),c.Ticks(end)]),0];
c.Label.Position = [2.3, 10^mean([log10(c.Ticks(1)),log10(c.Ticks(end))]),0];
set(gca,'xGrid','off')
pbaspect([2 1 1])

end
%%
CM = rocket
for i = 1:length(RPOIs)
figure()
Z = maxVRfluc(:,:,i);
surf(X,Y,Z,'EdgeColor','none')
xlabel('\theta_{dp}')
ylabel('R [mm]')

set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'})
set(gca,'YTick',[0,30,65,RPOIs(i)*1000])
colormap(CM)
xlim([0,pi]);
ylim([0,RPOIs(end)*1000]);
c = colorbar;
c.Label.String = '|\Psi(f=1MHz)| [pV]';
set(gca,'colorscale','log')
levels = [min(Z(:)),c.Ticks,max(Z(:))];
c.Ticks = levels;
c.TickLabels = arrayfun(@num2str,round(c.Ticks,3),'UniformOutput',false);
set(c,{'fontsize'},{10});
c.Label.FontSize = 12;
view(0,90)
hold on
lend = 2;
levels = interp1(1:lend,levels(1:lend),linspace(1,lend,10));
cont = contourc(unique(X(:)),unique(Y(:)),Z,levels);
plotContour(cont,[0.3,0.3,0.3],0)
hold off
set(findall(gcf,'type','axes'),'fontsize',11)
set(gcf,{'units','color','position','paperunits','papersize','Renderer'},...
    {'centimeters',[1,1,1],[1,3,7.5,6.5],'centimeters',[1+7.5,3+6.5],'Painters'})
%c.Label.Position = double(i==1)*[2.1    1.4986         0]+double(i==2)*[2.4222 0.2412 0]+double(i==3)*[2.1, mean([c.Ticks(1),c.Ticks(end)]),0];
c.Label.Position = [2.3, 10^mean([log10(c.Ticks(1)),log10(c.Ticks(end))]),0];
set(gca,'xGrid','off')

end

%% Compare with Wang et al
Tend = 1*fus^-1;
Tsim = 0:dt:Tend;

Theta = [0,pi/4,pi/2-0.01,pi/2-0.1,0.99*pi/2,pi/2,1.01*pi/2,pi/2+0.01,pi/2+0.1];
%Theta = [0:pi/300:3*pi/4]
dPhi = pi/4;
[phi,theta]=meshgrid(pi/2,Theta)
phisc = phi(:); thetasc = theta(:);
POIs=[];
d = 0.005 

Angles = [0,pi/2-0.01,pi/2-0.1,0.99*pi/2,pi/2,1.01*pi/2,pi/2+0.01,pi/2+0.1];
Rvals = 0;


VR = [];
VRfluc = [];
maxVRfluc = nan(length(Rvals),length(Angles),length(RPOIs));
idxOscDip = 1;
Pflag = 0;
VRflucs = {}
F = 1/(2*pi*1.5e12*Aus*1e-6*I*d*1000);
tic
for iscale=1:length(RPOIs(1))

    RPOI_sel = 1.5*d;
    titles = {'POIs on Brain','POIs on Scalp','POIs 5mm from scalp'};

for iR = 1:length(Rvals)
    for iorien = 1:length(Angles)
        % Select Sources and sinks + select vibrator
        OrienDipole = [0,sin(Angles(iorien)),cos(Angles(iorien))];
        LocDipole = Rvals(iR).*[0,sin(Angles(iorien)),cos(Angles(iorien))];
        CSource = LocDipole+d/2*OrienDipole;
        CSink = LocDipole-d/2*OrienDipole; 
        
        POIs = RPOI_sel.*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];

        % Run simulation
        Tsim = 0:dt:Tend;
        VR=SimDipoleOsc(Tend,CSource,CSink,I,USwave,1/fus,idxOscDip,'POI',POIs,...
            'resUS',resUS,'Display',double(iorien==1&&iR==1),'SolutionType',SolutionType,'scale',0,'ShowSphere','','fDependence', fDependence);
        VRfluc = (max(VR,[],2)-min(VR,[],2))/2*10^6;
        VRflucs{end+1} = VRfluc*F;
        [maxVRfluc(iR,iorien,iscale),idxmaxVRfluc]=max(VRfluc);

    end
    PComp = round(iR/length(Rvals)*100);
    if PComp>=Pflag
        disp([num2str(PComp),'%'])
        Pflag=Pflag+1;
    end
end
end
toc
%%
function plotContour(cont,color,logflag)
i = 1;
while i<size(cont,2)
    zval = cont(1,i);
    if logflag; zval= 10^(zval); end
    N1 = i+cont(2,i);
    plot3(cont(1,i+1:N1),cont(2,i+1:N1),zval*1.1*ones(1,N1-i),'color',color)
    i = N1+1;
end
end
