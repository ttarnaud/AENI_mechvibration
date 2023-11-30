function [VRMATSphere_tot,VRMATPOI]=PotentialSphere_Multi(CSource,CSink,AmpI,sigma,RSphere,ResSphere,Fignr,varargin)
%This function calculates the potential on a sphere due to multiple dipoles
%with sources and sinks in Csource and Csink, respectively. The amplitude
%of the dipole current is given in AmpI either one number (applies to all
%dipoles or row of equal length as Csource/Csink) in µA!!
%OUTPUT: VRMATSphere contains potential at the sphere boundary uV!!
%        VRMATPOI_tot: potential at POIs in uV!!
%Input:  CSource,CSink,AmpI,sigma,RSphere,ResSphere (See
%PotentialSingleSource)
%        Fignr: add fignr if is empty no plot generated
CM = 'jet';    %colormap
SolutionType = '3Sphere';
POIs = [];      %declaration POIs by default no POIs only the potential at sphere is calculated
PlotSS_flag = 0; %plot single source see end of potentialSingleSource
PLOT = true;    %create plot below(contains all Dipoles
DOIs = [];      %dipoles of interest (give different color)
OSCs = [];      % oscillators different collor
scale_flag = 1; %scale position POIs to outer sphere
ndipole_flag = 0; %at number to dipoles
VRMATPOI = [];


% other sphere thicknesses
tSkull = 0.005;   %thickness of skull
tScalp = 0.007;     %thickness of scalp
tAir = 0.005;

% Declare current in dipoles
idxI = ones(size(CSource,1),1);
Istandard = 1;          % standard is 1 µA
I = Istandard*ones(size(CSource,1),1);

if isempty(Fignr)
    PLOT = false;
end
if iscell(varargin) && length(varargin) == 1 && iscell(varargin{1})
    varargin = varargin{1};
end
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'idxI'))
        idxI = varargin{find(strcmpi(varargin,'idxI'))+1};
        end
        if any(strcmpi(varargin,'SolutionType'))
            SolutionType = varargin{find(strcmpi(varargin,'SolutionType'))+1};
        end
        if any(strcmpi(varargin,'POI'))
            POIs = varargin{find(strcmpi(varargin,'POI'))+1};
            VRMATPOI_intm = zeros(size(POIs,1),size(CSource,1));
        end
        if any(strcmpi(varargin,'DOI'))
            DOIs = varargin{find(strcmpi(varargin,'DOI'))+1};            
        end
        if any(strcmpi(varargin,'OSC'))
            OSCs = varargin{find(strcmpi(varargin,'OSC'))+1};
        end
        if any(strcmpi(varargin,'scale'))
            scale_flag = varargin{find(strcmpi(varargin,'scale'))+1};
            idx_scale_flag = find(strcmpi(varargin,'scale'))+1;
        end
        if any(strcmpi(varargin,'tSkull'))
            tSkull = varargin{find(strcmpi(varargin,'tSkull'))+1};            
        end
        if any(strcmpi(varargin,'tScalp'))
            tScalp = varargin{find(strcmpi(varargin,'tScalp'))+1};            
        end
        if any(strcmpi(varargin,'tAir'))
            tAir = varargin{find(strcmpi(varargin,'tAir'))+1};            
        end
        if any(strcmpi(varargin,'ndipole_flag'))
            ndipole_flag = varargin{find(strcmpi(varargin,'ndipole_flag'))+1}; 
        end
    end
else
    error('incorrect input')
end
% Generate points of interest on sphere by default outer sphere
[xs,ys,zs] = sphere(ResSphere);

if strcmpi(SolutionType,'3Sphere')
    xsphere = (RSphere+tSkull+tScalp)*xs; ysphere = (RSphere+tSkull+tScalp)*ys; zsphere = (RSphere+tSkull+tScalp)*zs;
elseif strcmpi(SolutionType,'4Sphere')
    xsphere = (RSphere+tSkull+tScalp+tAir)*xs; ysphere = (RSphere+tSkull+tScalp+tAir)*ys; zsphere = (RSphere+tSkull+tScalp+tAir)*zs;
else
    xsphere = RSphere*xs; ysphere = RSphere*ys; zsphere = RSphere*zs;
end

% Adjust current in dipoles
if length(AmpI)==sum(idxI) || length(AmpI)==1
I(logical(idxI)) = AmpI;
else 
    error('size AmpI and idxI do not match')
end

%Declaration of VRMAT
%3 dimensional array for each CSource the potential at the POIs and sphere
%points is saved. xsphere already in surf format (2D)
VRMAT = zeros(size(xsphere,1),size(xsphere,2));
if ~isempty(POIs)
    if size(POIs,2)==3
        VRMATPOI = zeros(size(POIs,1),1);
    end
end

for ic = 1:size(CSource,1)
    % current is given in absolute values
    if I(ic)<0
    Inputdipole.CSource = CSink(ic,:);
    Inputdipole.CSink = CSource(ic,:);
    else
    Inputdipole.CSource = CSource(ic,:);
    Inputdipole.CSink = CSink(ic,:);
    end
    
    VRMATinterm = PotentialSingleSource(Inputdipole,abs(I(ic)*10^-6),sigma,xsphere,ysphere,zsphere,RSphere,'',SolutionType,horzcat(varargin,{'plot',PlotSS_flag}));
    VRMAT = VRMAT + reshape(VRMATinterm,size(xsphere));
    
    if ~isempty(POIs)
        if size(POIs,2)==3
            VRMATPOIs_intm = PotentialSingleSource(Inputdipole,abs(I(ic)*10^-6),sigma,POIs(:,1),POIs(:,2),POIs(:,3),RSphere,'',SolutionType,horzcat(varargin,{'plot',PlotSS_flag}));
            VRMATPOI = VRMATPOI+VRMATPOIs_intm*1e6;
        else
            error('size POIs incorrect')
        end
    end
    
    if PLOT
        %place dipole arrows, give collor depending on OSC(blue), DOI(red)
        %static (black) see arrow3d for settigns. Volume of arrow depends
        %on strength of current
    figure(Fignr)
    
    d = Inputdipole.CSource-Inputdipole.CSink;
    cdipole = 1/2*(Inputdipole.CSource+Inputdipole.CSink);
    Arrowpoints = [(cdipole-d/norm(d)*RSphere/15)',(cdipole+d/norm(d)*RSphere/15)'];
    Iarrow = I(ic)*5;
    ARROW = arrow3d(Arrowpoints(1,:),Arrowpoints(2,:),Arrowpoints(3,:),0.75,2*(Iarrow/10+1)*RSphere/400,(Iarrow/10+1)*6*RSphere/400);
    set(ARROW,'facecolor',[0,0,0])
    if ic==1; hold on; end
    if ic==size(CSource,1);hold off; end
    if ~isempty(OSCs)
        if any(ic==OSCs)
            set(ARROW,'facecolor',[0 0 1])
        end
    end
    if ~isempty(DOIs)
        if any(ic==DOIs)
            set(ARROW,'facecolor',[1 0 0])
        end
    end
    if ndipole_flag
        text(cdipole(1),cdipole(2),cdipole(3)+0.01,num2str(ic),'FontSize',10,'Color','r');
    end
    end
    
end

VRMATSphere_tot = VRMAT*10^6;%µV

if strcmpi(SolutionType,'3Sphere')
    %model only evaluated at outer sphere
    %xsphere = (RSphere+tSkull+tScalp)*xs; ysphere = (RSphere+tSkull+tScalp)*ys; zsphere = (RSphere+tSkull+tScalp)*zs;
    POIs = POIs./vecnorm(POIs,2,2).*(RSphere+tSkull+tScalp);
elseif strcmpi(SolutionType,'4Sphere')
    %xsphere = (RSphere+tSkull+tScalp+tAir)*xs; ysphere = (RSphere+tSkull+tScalp+tAir)*ys; zsphere = (RSphere+tSkull+tScalp+tAir)*zs;
    if scale_flag
        POIs = POIs./vecnorm(POIs,2,2).*(RSphere+tSkull+tScalp+tAir);
    end
end
% Plot Voltage on sphere
if PLOT
figure(Fignr)
hold on
surf(xsphere,ysphere,zsphere,VRMATSphere_tot,'FaceAlpha',0.2,'EdgeColor','none','SpecularStrength',0);
hold off
view(120,30)
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
c.Label.String = 'Potential at sphere boundary [µV]';
c.Label.FontSize = 14;
colormap(CM)
%LIMS = [floor(min(VRMATSphere_tot(:))),ceil(max(VRMATSphere_tot(:)))];
%LIMS = [-0.4,0.4];
LIMS = [min(VRMATSphere_tot(:)),max(VRMATSphere_tot(:))];
if any(LIMS)
set(hAxis,'CLIM',LIMS);
set(c,'Limits',LIMS);
end
%set(gcf,'position',[-1919,41,1920,963])
end
%Calculate VR in POIs and draw electrodes
if ~isempty(POIs)
    DrawElectrodes(POIs,'radius',0.01)    
end
hold on
h = zeros(3, 1);
h(1) = plot(nan,nan,'sr','MarkerFaceColor','r');
h(2) = plot(nan,nan,'sb','MarkerFaceColor','b');
h(3) = plot(nan,nan,'sk','MarkerFaceColor','k');
hold off
legend(h, {'DOI','Oscillating Dipoles','Static Dipoles'},'location','northwest');
title([num2str(size(CSource,1)),' dipoles']);
end