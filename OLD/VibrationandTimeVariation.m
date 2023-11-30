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
Alphafun = @(t,Tau) 2*I*(t./Tau).*exp(1-t./Tau);
Tau = 0.001;

% US options
fus = 1e6; %Hz
wus = 2*pi*fus; %rad/s
Aus = 10e-9; %m
Phaseus = 0; %rad
Dirus = [0,0,1]; Dirus = Dirus/norm(Dirus);

% Simulation settings
dt = 1/40*(fus)^-1;
Tend =5*Tau;

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
POIsSphere=[];

    
        
 POIsSphere = RSphere*[cos(phisc).*sin(thetasc),sin(phisc).*sin(thetasc),cos(thetasc)];
    


%% dipole in center and yz-direction
OrienDipole = [0,0,1];
idxOscDip = 1;
POIs = RSphere*[0,0,1;0,1,0;0,1/sqrt(2),1/sqrt(2)];
figure(Fignr)
ax0 = subplot(4,3,1);
PotentialSphere_Multi(d/2*OrienDipole,-d/2*OrienDipole,Alphafun(0.001,Tau),sigma,RSphere,1000,Fignr,'POI',POIs);

% Select Sources and sinks + select vibrator
CSource = d/2*OrienDipole;
CSink = -d/2*OrienDipole;
CSourceus = CSource(idxOscDip,:);
CSinkus = CSink(idxOscDip,:);
% Run simulation
Tsim = 0:dt:Tend;
Pflag = 0;
for i=1:length(Tsim)
    t = Tsim(i);
    CSource(idxOscDip,:) = CSourceus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    CSink(idxOscDip,:) = CSinkus + Dirus.*Aus.*sin(wus.*t+Phaseus);
    [VRMATSphere_tot,VRMATPOI_tot]=PotentialSphere_Multi(CSource,CSink,Alphafun(t,Tau),sigma,RSphere,1,[],'POI',POIs);
    VR(:,i) = VRMATPOI_tot;
    
    if round(i/length(Tsim)*100,0) >= Pflag
    disp([num2str(Pflag),'%'])
    Pflag = Pflag+5;
    end
end

plotinfo(Tsim,Dirus,Aus,wus,Phaseus,[1,2,3],VR,Alphafun,Tau)
Fignr = get(gcf,'number');
Fignr = Fignr+1;
%%
function plotinfo(Tsim,Dirus,Aus,wus,Phaseus,POIs,VR,Alphafun,Tau)


subplot(4,3,2);
plot([Tsim(1:100)*1e6.*ones(3,1)]',[Dirus'.*Aus.*sin(wus.*Tsim(1:100)+Phaseus)]'*1e9)
legend({'x component','y component','z component'})
ylabel('displacement [nm]')
xlabel('Time [µs]')

subplot(4,3,3);
plot(Tsim*1e3,Alphafun(Tsim,Tau))
ylabel('Amplitude dipole [µA]')
xlabel('Time [ms]')

freqs = [0,logspace(0,7,71)];
for ifig =1:3
ax{ifig} = subplot(4,3,3*ifig+1);
nr = POIs(ifig);
plot(Tsim*1e3,VR(nr,:))
if ifig == 3
xlabel('Time [ms]')
end
if ifig == 1
title('Measured potential at POI [µV]')
end
ylabel({['\bf POI = ',num2str(nr)],''})

subplot(4,3,3*ifig+2)
spectrogram(VR(nr,:),128,120,freqs,1/(Tsim(2)-Tsim(1)),'yaxis');
[s,~,t]=spectrogram(VR(nr,:),128,120,freqs,1/(Tsim(2)-Tsim(1)),'yaxis');
set(gca,'Yscale','log')
xlabel('')
if ifig == 1
title('Spectrogram w=128 noverlap=120')
end
if ifig == 3
    xlabel('Time [ms]')
end
subplot(4,3,ifig*3+3)
plot(t*1e3,s(find(freqs==1e6),:))
if ifig == 3
xlabel('Time [ms]')
end
if ifig == 1
title('Short term fourier transform at f = fus')
end
end


end