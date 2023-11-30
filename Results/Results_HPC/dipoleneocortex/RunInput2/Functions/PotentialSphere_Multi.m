function [VRMATSphere_tot,VRMATPOI_tot]=PotentialSphere_Multi(CSource,CSink,AmpI,sigma,RSphere,ResSphere,Fignr,varargin)
%This function calculates the potential on a sphere due to multiple dipoles
%with sources and sinks in Csource and Csink, respectively. The amplitude
%of the dipole current is given in AmpI either one number (applies to all
%dipoles or row of equal length as Csource/Csink) in µA!!

CM = 'jet';    %colormap
SolutionType = '3Sphere';
POIs = [];      %declaration POIs
PlotSS_flag = 0;
PLOT = true;
DOIs = [];      %dipoles of interest (give different color)
scale_flag = 1;

% Generate points of interest on sphere
[xs,ys,zs] = sphere(ResSphere);
xsphere = RSphere*xs; ysphere = RSphere*ys; zsphere = RSphere*zs;
% other sphere thicknesses
tSkull = 0.005;
tScalp = 0.007;
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
            VRMATPOIs = zeros(size(POIs,1),size(CSource,1));
        end
        if any(strcmpi(varargin,'DOI'))
            DOIs = varargin{find(strcmpi(varargin,'DOI'))+1};            
        end
        if any(strcmpi(varargin,'scale'))
            scale_flag = varargin{find(strcmpi(varargin,'scale'))+1};
        if any(strcmpi(varargin,'tSkull'))
            tSkull = varargin{find(strcmpi(varargin,'tSkull'))+1};            
        end
        if any(strcmpi(varargin,'tScalp'))
            tScalp = varargin{find(strcmpi(varargin,'tScalp'))+1};            
        end
        if any(strcmpi(varargin,'tAir'))
            tAir = varargin{find(strcmpi(varargin,'tAir'))+1};            
        end
        end
    end
else
    error('incorrect input')
end

% Adjust current in dipoles
if length(AmpI)==sum(idxI) || length(AmpI)==1
I(logical(idxI)) = AmpI;
else
    error('size AmpI and idxI do not match')
end

%Declaration of VRMAT
VRMAT = zeros(size(xsphere,1),size(xsphere,2),size(CSource,1));

for i = 1:size(CSource,1)
    if I(i)<0
    Inputdipole.CSource = CSink(i,:);
    Inputdipole.CSink = CSource(i,:);
    else
    Inputdipole.CSource = CSource(i,:);
    Inputdipole.CSink = CSink(i,:);
    end
    
    VRMATinterm = PotentialSingleSource(Inputdipole,abs(I(i)*10^-6),sigma,xsphere,ysphere,zsphere,RSphere,'',SolutionType,horzcat(varargin,{'plot',PlotSS_flag}));
    VRMAT(:,:,i) = reshape(VRMATinterm,size(xsphere));
    
    if ~isempty(POIs)
        if size(POIs,2)==3
            VRMATPOIs(:,i) = PotentialSingleSource(Inputdipole,abs(I(i)*10^-6),sigma,POIs(:,1),POIs(:,2),POIs(:,3),RSphere,'',SolutionType,horzcat(varargin,{'plot',PlotSS_flag}));
            
        else
            error('size POIs incorrect')
        end
    end
    
    if PLOT
    figure(Fignr)
    hold on
    d = Inputdipole.CSource-Inputdipole.CSink;
    cdipole = 1/2*(Inputdipole.CSource+Inputdipole.CSink);
    Arrowpoints = [(cdipole-d/norm(d)*RSphere/20)',(cdipole+d/norm(d)*RSphere/20)'];
    
    ARROW = arrow3d(Arrowpoints(1,:),Arrowpoints(2,:),Arrowpoints(3,:),0.75,2*(I(i)/10+1)*RSphere/400,(I(i)/10+1)*6*RSphere/400);
    set(ARROW,'facecolor',[0,0,0])
    hold off 
    if ~isempty(DOIs)
        if any(i==DOIs)
            set(ARROW,'facecolor',[1 0 0])
        end
    end
    end
    
end

VRMATSphere_tot = sum(VRMAT,3)*10^6;%µV
if strcmpi(SolutionType,'3Sphere')
    xsphere = (RSphere+tSkull+tScalp)*xs; ysphere = (RSphere+tSkull+tScalp)*ys; zsphere = (RSphere+tSkull+tScalp)*zs;
    POIs = POIs./RSphere.*(RSphere+tSkull+tScalp);
elseif strcmpi(SolutionType,'4Sphere')
    xsphere = (RSphere+tSkull+tScalp+tAir)*xs; ysphere = (RSphere+tSkull+tScalp+tAir)*ys; zsphere = (RSphere+tSkull+tScalp+tAir)*zs;
    if scale_flag
        POIs = POIs./RSphere.*(RSphere+tSkull+tScalp+tAir);
    end
end
% Plot Voltage on sphere
if PLOT
figure(Fignr)
hold on
surf(xsphere,ysphere,zsphere,VRMATSphere_tot,'FaceAlpha',0.5,'EdgeColor','none','SpecularStrength',0);
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
set(gcf,'position',[-1919,41,1920,963])
end
%Calculate VR in POIs and draw electrodes
if ~isempty(POIs)
    VRMATPOI_tot = sum(VRMATPOIs,2)*10^6;%µV
    if PLOT
        DrawElectrodes(POIs)
    end
else
    VRMATPOI_tot = [];
end
title([num2str(size(CSource,1)),' dipoles']);
end