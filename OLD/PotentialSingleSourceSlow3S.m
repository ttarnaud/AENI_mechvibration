function [VRMAT] = PotentialSingleSourceSlow3S(Inputdipole,dipoleI,sigma,xPOI,yPOI,zPOI,RSphere,Inputtype,SolutionType,varargin)
% this function contains alternative to calculate the 3 sPherical model
% however is apperently slower :(
% default values
RatioSkT = 25;
tSkull = 0.005; %skull thickness [m]
tScalp = 0.007; %Scalp thickness [m]
ShowResult = 0;
CM = 'jet';
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        idxShowResult = find(strcmpi(varargin,'plot'))+1;
        if ~isempty(idxShowResult)
        ShowResult = varargin{idxShowResult};
        end
        idxRatioSkullTissue = find(strcmpi(varargin,'RatioSkT'))+1;
        if ~isempty(idxRatioSkullTissue)
        RatioSkT = varargin{idxRatioSkullTissue};
        end
        idxtSkull = find(strcmpi(varargin,'tSkull'))+1;
        if ~isempty(idxtSkull)
        tSkull = varargin{idxtSkull};
        end        
        idxtScalp = find(strcmpi(varargin,'tScalp'))+1;
        if ~isempty(idxtScalp)
        tScalp = varargin{idxtScalp};
        end
    end
else
    error('Input incorrect')
end

NORMS = vecnorm([xPOI(:),yPOI(:),zPOI(:)],2,2);
if any(round(NORMS-RSphere,-(floor(log10(RSphere))-6)))
    error('POIs not on Sphere')
end

switch lower(Inputtype)
    case 'euler'
        Posdipole = Inputdipole.Posdipole;
        dipoleTheta = Inputdipole.dipoleTheta;
        dipolePhi = Inputdipole.dipolePhi;
        dipoled = Inputdipole.dipoled;
        
        CSource = Posdipole+dipoled*[cos(dipolePhi).*sin(dipoleTheta),sin(dipolePhi).*sin(dipoleTheta),cos(dipoleTheta)];
        CSink = Posdipole-dipoled*[cos(dipolePhi).*sin(dipoleTheta),sin(dipolePhi).*sin(dipoleTheta),cos(dipoleTheta)];
        if norm(CSource,2) >= RSphere || norm(CSink,2) >= RSphere
            error('dipole not in sphere')
        end
    case 'spherical'
        Posdipole = Inputdipole.Posdipole;
        dipoleTheta = Inputdipole.dipoleTheta;
        dipolePhi = Inputdipole.dipolePhi;
        dipoled = Inputdipole.dipoled;
        
        Posdipolexyz = Posdipole(1).*[cos(Posdipole(2)).*sin(Posdipole(3)),sin(Posdipole(2)).*sin(Posdipole(3)),cos(Posdipole(3))];
        CSource = Posdipolexyz+dipoled*[cos(dipolePhi).*sin(dipoleTheta),sin(dipolePhi).*sin(dipoleTheta),cos(dipoleTheta)];
        CSink = Posdipolexyz-dipoled*[cos(dipolePhi).*sin(dipoleTheta),sin(dipolePhi).*sin(dipoleTheta),cos(dipoleTheta)];
        if norm(CSource,2) >= RSphere || norm(CSink,2) >= RSphere
            error('dipole not in sphere')
        end
    otherwise
        CSource = Inputdipole.CSource;
        CSink = Inputdipole.CSink;
        dipoled = norm((CSource-CSink),2);
        if norm(CSource,2) >= RSphere || norm(CSink,2) >= RSphere
            error('dipole not in sphere')
        end
end
if dipoleI<0
    CSource_temp = CSource;
    CSource = CSink;
    CSink = CSource_temp;
    clear('CSource_temp')
end

if ~any(CSource-CSink)
    error('CSource = CSink')
end

switch lower(SolutionType)
    case 'closedbounded'
        ra = @(x,y,z,CSink) sqrt((x-CSink(1)).^2+(y-CSink(2)).^2+(z-CSink(3)).^2);
        rb = @(x,y,z,CSource) sqrt((x-CSource(1)).^2+(y-CSource(2)).^2+(z-CSource(3)).^2);
        cosB = @(x,y,z,CSink) (ra(x,y,z,CSink).^2-(x.^2+y.^2+z.^2)-norm(CSink).^2)./(-2*sqrt((x.^2+y.^2+z.^2))*norm(CSink));
        cosT = @(x,y,z,CSource) (rb(x,y,z,CSource).^2-(x.^2+y.^2+z.^2)-norm(CSource).^2)./(-2*sqrt((x.^2+y.^2+z.^2))*norm(CSource));
        VR = @(x,y,z,CSource,CSink) abs(dipoleI)/(4*pi*sigma)*(2./(rb(x,y,z,CSource))-2./(ra(x,y,z,CSink))+...
            1/RSphere.*log((ra(x,y,z,CSink)+RSphere-norm(CSink).*cosB(x,y,z,CSink))./(rb(x,y,z,CSource)+RSphere-norm(CSource).*cosT(x,y,z,CSource))));
        VRMAT = VR(xPOI(:),yPOI(:),zPOI(:),CSource,CSink);
        
    case 'unbounded'
        ra = @(x,y,z,CSink) sqrt((x-CSink(1)).^2+(y-CSink(2)).^2+(z-CSink(3)).^2);
        rb = @(x,y,z,CSource) sqrt((x-CSource(1)).^2+(y-CSource(2)).^2+(z-CSource(3)).^2);
        VR = @(x,y,z,CSource,CSink) abs(dipoleI)./(4*pi*sigma).*(1./(rb(x,y,z,CSource))-1./ra(x,y,z,CSink));
        VRMAT = VR(xPOI(:),yPOI(:),zPOI(:),CSource,CSink);
    case '3sphere'
        iter_flag = true; % flag for while loop
        istep = 50;
        idx_lim = 1e3;
        PosCD = (CSource+CSink)./2; % Position centre of dipole
        % if PosCD == [0,0,0] ==> ABC turns zero for i>1 thus iteration not
        % necessary
        if ~any(PosCD)
            istep = 2; %if set to 1 would return in array error therefore 2
        end
        iarray = 1:istep;
        DM = dipoleI*(CSource-CSink); %dipole moment
        RScalp = RSphere;  % distance to outerside of head boundary scalp/air
        RSkull = RScalp-tScalp; % distance to boundary skull/scalp
        RBrain = RSkull-tSkull; % distance to boundary brain/skull
        f1 = RBrain/RScalp; % normalized distance brain
        f2 = RSkull/RScalp; % normalized distance skull
        if any(PosCD) 
            nPosCD = PosCD/norm(PosCD); % normalized pos dipole
        else
            nPosCD = CSource/norm(CSource); % normalized pos dipole
        end
        POI = [xPOI(:),yPOI(:),zPOI(:)]; % Vector POIs [Nx3]
        nPOI = POI./vecnorm(POI,2,2); % normalized vector POIs [Nx3]
        CosT = (nPosCD*nPOI')'; % cosine of angle between dipole and POI
        CosT = min(CosT,ones(size(CosT)));
        CosT = max(CosT,-ones(size(CosT)));
        tPosCD = cross(cross(repmat(PosCD,size(POI,1),1),POI,2),repmat(PosCD,size(POI,1),1),2); % tangential of dipole [Nx3]
        ntPosCD = tPosCD./vecnorm(tPosCD,2,2); % normalized tangential
        ntPosCD(isnan(ntPosCD))=0;
        nPosCD = repmat(nPosCD,size(POI,1),1); %from [1x3] to [Nx3]
        gi = @(iSum) ((iSum+1).*RatioSkT+iSum).*(iSum.*RatioSkT./(iSum+1)+1)+(1-RatioSkT).*((iSum+1).*RatioSkT+iSum).*(f1.^(2*iSum+1)-f2.^(2*iSum+1))...
            -iSum.*(1-RatioSkT).^2.*(f1/f2).^(2*iSum+1);
        ABC = @(iSum) RatioSkT*(2*iSum+1).^3./(gi(iSum).*(iSum+1).*iSum).*(norm(PosCD,2)./RScalp).^(iSum-1);
        Gix = @(R,T,iSum,Pi) ABC(iSum).*(iSum.*Pi(:,:,1)'.*R(:,1)-Pi(:,:,2)'.*T(:,1));
        Giy = @(R,T,iSum,Pi) ABC(iSum).*(iSum.*Pi(:,:,1)'.*R(:,2)-Pi(:,:,2)'.*T(:,2));
        Giz = @(R,T,iSum,Pi) ABC(iSum).*(iSum.*Pi(:,:,1)'.*R(:,3)-Pi(:,:,2)'.*T(:,3));
        % start iteration
        
        Pinew = LegendreArray(iarray,CosT,[0,1]);
        Pitest = Pinew(1:end-1,:,:);
        Gitest = 1/(4*pi*sigma*RScalp.^2).*[sum(Gix(nPosCD,ntPosCD,iarray(1:end-1),Pitest),2),...
            sum(Giy(nPosCD,ntPosCD,iarray(1:end-1),Pitest),2),sum(Giz(nPosCD,ntPosCD,iarray(1:end-1),Pitest),2)];
        Ginew = 1/(4*pi*sigma*RScalp.^2).*[sum(Gix(nPosCD,ntPosCD,iarray,Pinew),2),...
            sum(Giy(nPosCD,ntPosCD,iarray,Pinew),2),sum(Giz(nPosCD,ntPosCD,iarray,Pinew),2)];
        VRMATtest = Gitest*DM';
        VRMAT = Ginew*DM';
        if max(abs(VRMAT-VRMATtest)./VRMATtest) < 1e-6 || max(abs(VRMAT-VRMATtest)) == 0
                iter_flag = false;
        end
        while iter_flag && iarray(end)<idx_lim
            iarrayold = iarray;
            iarray = 1:iarrayold(end)+istep;
            Piadd = LegendreArray(iarrayold(end)+1:iarray(end),CosT,[0,1]);
            Pitest = cat(1,Pinew,Piadd(1:end-1,:,:));
            Pinew = cat(1,Pinew,Piadd(1:end,:,:));
            Gitest = 1/(4*pi*sigma*RScalp.^2).*[sum(Gix(nPosCD,ntPosCD,iarray(1:end-1),Pitest),2),...
                sum(Giy(nPosCD,ntPosCD,iarray(1:end-1),Pitest),2),sum(Giz(nPosCD,ntPosCD,iarray(1:end-1),Pitest),2)];
            Ginew = 1/(4*pi*sigma*RScalp.^2).*[sum(Gix(nPosCD,ntPosCD,iarray,Pinew),2),...
                sum(Giy(nPosCD,ntPosCD,iarray,Pinew),2),sum(Giz(nPosCD,ntPosCD,iarray,Pinew),2)];
            VRMATtest = Gitest*DM';
            VRMAT = Ginew*DM';

            if max(abs(VRMAT-VRMATtest)./VRMATtest) < 1e-6 || max(abs(VRMAT-VRMATtest)) == 0
                iter_flag = false;
            end
            
            if iarray(end) == idx_lim
                disp('max iterations reached')
                disp(['max add value',num2str(max(abs(VRMAT-VRMATtest)./VRMATtest))])
                break
            end
            if ~boolean(mod(iarray(end),100))
                disp(num2str(iarray(end)))
            end
                    
%         gi(idx)
%         Ginew
%         Pi(1,:)
%         repmat(nPosCD,size(POI,1),1)*DM'
%         ntPosCD*DM'
%         Pi(2,:)
        end
        
    otherwise
        disp('to be continued')
        VRMAT = 0;
end
if ShowResult
    disp('POIs must create a sphere')
    VRMAT_Plot = reshape(VRMAT,size(xPOI));
    figure
    arrow3d([CSink(1),CSource(1)],[CSink(2),CSource(2)],[CSink(3),CSource(3)],0.75,RSphere/100);
    hold on
    surf(xPOI,yPOI,zPOI,VRMAT_Plot*1e6,'FaceAlpha',0.5,'EdgeColor','none','SpecularStrength',0);
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
    c = colorbar;
    c.Label.String = 'Voltage at sphere boundary [µV]';
    c.Label.FontSize = 18;
    colormap(CM);
    
    %LIMS = [floor(min(VRMAT(:)*1e6)),ceil(max(VRMAT(:)*1e6))];
    LIMS = [min(VRMAT(:)*1e6),max(VRMAT(:)*1e6)];
    if exist('LIMS','var')
        set(hAxis,'CLIM',LIMS);
        set(c,'Limits',LIMS);
    end
end

end
