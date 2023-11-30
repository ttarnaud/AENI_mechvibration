% Run single dipole
addpath(genpath('D:\users\rschoeters\Documents\Imec USEEG\Matlab\Functions'))
if exist('Fignr','var')
    clearvars -except 'Fignr'
    Fignr = Fignr+1;
else
    clear all
    Fignr = 10;
end

GenGIF=false;
if ~GenGIF
    disp('no gif generation')
end
resSphere = 100;
dipoleI = 10;    % [µA]
sigma = 0.33;   %[S/m]
RSphere = 0.058; %[m]
[x,y,z] = sphere(resSphere); %POIs
xsphere = x*RSphere; ysphere = y*RSphere; zsphere = z*RSphere; %Scaled POIs
Inputtype = '';     %Input = exact locations Csource and Csink
SolutionType = '3Sphere';
CSources = [0,0,1;0,1,0;0,1,1;0,1,0];
CSinks = [0,0,1;0,1,0;0,1,1;0,0,1];
for i=1%:size(CSources)
Inputdipole.CSource = 0.01*CSources(i,:)/norm(CSources(i,:));
Inputdipole.CSink = 0.0095*CSinks(i,:)/norm(CSinks(i,:));

VRMAT = PotentialSingleSource(Inputdipole,dipoleI*1e-6,sigma,xsphere,ysphere,zsphere,RSphere,Inputtype,SolutionType,'plot',1);
end
%% Run comparison Bounded unbounded
Inputdipole.CSource = (0.065+0.00025)*CSources(1,:)/norm(CSources(1,:));
Inputdipole.CSink = (0.065-0.00025)*CSinks(1,:)/norm(CSinks(1,:));
figure(Fignr)
ax1 = subplot(1,3,1);
VRMATSphere_tot = PotentialSphere_Multi(Inputdipole.CSource,Inputdipole.CSink,dipoleI,sigma,0.07,resSphere,Fignr,'SolutionType','ClosedBounded');
LIMS(1,:) = get(ax1.Colorbar,'Limits');
title(['Bounded Medium: max = ',num2str(max(VRMATSphere_tot(:)))])
ax2 = subplot(1,3,2);
VRMATSphere_tot = PotentialSphere_Multi(Inputdipole.CSource,Inputdipole.CSink,dipoleI,sigma,0.07,resSphere,Fignr,'SolutionType','unBounded');
LIMS(2,:) = get(ax2.Colorbar,'Limits');
title(['Unbounded Medium: max = ',num2str(max(VRMATSphere_tot(:)))])
ax3 = subplot(1,3,3);
VRMATSphere_tot = PotentialSphere_Multi(Inputdipole.CSource,Inputdipole.CSink,dipoleI,sigma,0.058,resSphere,Fignr,'SolutionType','3Sphere');
LIMS(3,:) = get(ax3.Colorbar,'Limits');
title(['3Sphere model: max = ',num2str(max(VRMATSphere_tot(:)))])
LIMS = [min(LIMS(:,1)),max(LIMS(:,2))];
set(ax1,'CLIM',LIMS);
set(ax1.Colorbar,'Limits',LIMS);
set(ax1.Colorbar,'Visible', 'off')
set(ax2,'CLIM',LIMS);
set(ax2.Colorbar,'Limits',LIMS);
set(ax2.Colorbar,'Visible', 'off')
set(ax3,'CLIM',LIMS);
set(ax3.Colorbar,'Limits',LIMS);
%% Multiple dipoles
clearvars -except Fignr GenGIF
Fignr = Fignr+1;
RSphere = 0.07;
sigma = 0.33;
resSphere = 100;
[xs,ys,zs] = sphere(resSphere);
xsphere = xs*RSphere; ysphere = ys*RSphere; zsphere = zs*RSphere;
Inputtype = '';
SolutionType = 'ClosedBounded';
CM = 'jet';
d = 0.002;
Rd = 0.03;
%% Gen Dipole locations (concentric away from center)
[xd,yd,zd] = sphere(3);
LocDipole = [xd(:),yd(:),zd(:)];
LocDipole = unique(LocDipole,'rows');
CSource = (Rd+d/2)*LocDipole;
CSink = (Rd-d/2)*LocDipole;
I = 1; %[µA]
figure(Fignr)
subplot(1,2,1)
PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,resSphere,Fignr);
Fignr = get(gcf,'number');
title('Orientation: Radial Out')
%% Gen Dipole locations (all up)
subplot(1,2,2)
[xd,yd,zd] = sphere(3);
LocDipole = Rd*[xd(:),yd(:),zd(:)];
LocDipole = unique(LocDipole,'rows');
CSource = LocDipole+[0,0,d/2];
CSink = LocDipole+[0,0,-d/2];
I = 1; %[µA]

PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,resSphere,Fignr);
Fignr = get(gcf,'number');
title('Orientation: Positive Z')
%% Gen Dipole locations (concentric away from center)
Fignr = Fignr+1;
idx=0;
clear('xd','yd','zd')
for i=[0:pi/6:pi]
    circum = Rd*sin(i)*2*pi;
    if (circum-1e-9) <= 0
        Angles = 0;
    else
        Amount = circum*6/(2*pi*Rd);
        interval = 2*pi/Amount;
        Angles = [0:interval:2*pi-interval];
    end
    for j= 1:length(Angles)
        idx=idx+1;
        xd(idx) = cos(Angles(j))*sin(i);
        yd(idx) = sin(Angles(j))*sin(i);
        zd(idx) = cos(i);
    end
end
LocDipole = [xd(:),yd(:),zd(:)];
LocDipole = unique(LocDipole,'rows');
CSource = (Rd+d/2)*LocDipole;
CSink = (Rd-d/2)*LocDipole;
I = 1;%µA

PotentialSphere_Multi(CSource,CSink,I,sigma,RSphere,resSphere,Fignr);

%% Generate Gif
%Gen Dipole locations (concentric away from center) +1 in center
if GenGIF
close all
Fignr = Fignr+1
VRMAT = [];
[xd,yd,zd] = sphere(3);
LocDipole = [xd(:),yd(:),zd(:)];
LocDipole = unique(LocDipole,'rows');
CSource = vertcat((Rd+d/2)*LocDipole,[0,0,d/2]);
CSink = vertcat((Rd-d/2)*LocDipole,[0,0,-d/2]);
idxI = zeros(length(CSource),1);
idxI(end) = 1;
I = 1; %µA
% Variable I in time
Alphafun = @(t,alpha) 2*I*(alpha.*t).*exp(1-alpha.*t);
time = 0:0.05:5;
filename = 'testAnimated.gif';
for n=1:length(time)    
PotentialSphere_Multi(CSource,CSink,Alphafun(time(n),1),sigma,RSphere,resSphere,Fignr+n,'idxI',idxI);
drawnow
h = gcf;
set(h,'position',[-1690,108,1419,851])
LIMS=[-1,2];
ax = gca;
set(ax,'CLIM',LIMS);
set(ax.Colorbar,'Limits',LIMS);
frame = getframe(h);
im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  % Write to the GIF File 
  if n == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
  else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
  end 
  close(h)
  disp([num2str(n/length(time)*100),'%']);
end
end
