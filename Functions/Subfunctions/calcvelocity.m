function [veloc, press, loc, fus, maxv,maxv_pos,maxv_dir] = calcvelocity(location,Filename,RBrain,Plot_flag,varargin)
%This function computes the velocity field from the pressure fields
%generated in Sim 4 life. 
%OUTPUT:
%        veloc: structure containing the amplitudes and phases of the
%        vibration of the three euclidean components if database flag:
%        contains also direction of maximum vibration and phase going with
%        this direction
%        press: structure with pressuremap and phase (according to
%        sinusoidal regime)
%        loc: structure containing X,Y,Z grids
%        fus: ultrasonic frequency used in S4l (default 1e6 otherwise
%        specify in varargin)
%        maxv: maximal vibration amplitude
%        maxv_pos: Position of maximal vibration
%        maxv_dir: normalized vibration vector at point of maximal vibration
%INPUT:
%       location: file location to load
%       Filename: what's in a name? That which we call a rose
%       RBrain: radius of modeled brain
%       Plot_flag: create plots or not
%       varargin: varable inputs 
%              'snapshot': pressure field snapshots at  certain time points
%              also add 'totalperiods', initTime
%               'fus: frequency used in sim
%               'downsf': downsample factor = reduce size (not ideal, not
%               via interpolation)
%               'plane': code of plots programmed for XY plane 
%               'scalefactor': linear solver => scale amplitude is possible
%               without creating errors
%               'rotate_flag: depending on position of SEFTS
if ~contains(Filename,'.mat')
    %if .mat included load specified file. if not look if also file with
    %time snapshots is present
    try
        Sim4life_input_t = load([location,Filename,'_l2pt.mat']);
        snapshot_flag = 1;
    catch
        snapshot_flag = 0;
    end
    if any(strcmpi(varargin,'snapshot'))
        snapshot_flag = varargin{find(strcmpi(varargin,'snapshot'))+1};
    end
    
    Sim4life_input_f = load([location,Filename,'_pf.mat']);
    Sim4life_input.p_x_y_z_f_Snapshot0 = Sim4life_input_f.Snapshot0;
    Sim4life_input.Axis0 = Sim4life_input_f.Axis0;
    Sim4life_input.Axis1 = Sim4life_input_f.Axis1;
    Sim4life_input.Axis2 = Sim4life_input_f.Axis2;
    if snapshot_flag
        fnmsSim4Lifet = fieldnames(Sim4life_input_t);
        nrSnapshots = length(fnmsSim4Lifet)-3;  %input('nrSnapshots? = ');
        for ifnms = 1:length(fnmsSim4Lifet)
            if contains(fnmsSim4Lifet{ifnms},'Snapshot')
                nrstr = fnmsSim4Lifet{ifnms}(length('Snapshot')+1:end);
                Sim4life_input.(['p_x_y_z_t_Snapshot',nrstr]) = Sim4life_input_t.(fnmsSim4Lifet{ifnms});
            end
        end
        clear('Sim4life_input_t')
    end
    %clear not to clogg RAM
    clear('Sim4life_input_f')
else
    Sim4life_input = load([location,Filename]);
    nrSnapshots = length(fieldnames(Sim4life_input))-7;
    snapshot_flag = 1;
    if any(strcmpi(varargin,'snapshot'))
        snapshot_flag = varargin{find(strcmpi(varargin,'snapshot'))+1};
    end
end


rhoBrain = 1044.5;
if any(strcmpi(varargin,'fus'))
    fus = varargin{find(strcmpi(varargin,'fus'))+1};
else
fus = 1e6;
end
wus = 2*pi*fus;
SinRegime = 1;
xres = 6;
yres = 20;
if any(strcmpi(varargin,'downsf'))
    %see intro
    downsf = varargin{find(strcmpi(varargin,'downsf'))+1};
else
    downsf = 1;
end

if any(strcmpi(varargin,'Database'))
    %current run based on database
   Database_flag = varargin{find(strcmpi(varargin,'Database'))+1};
else
    Database_flag = 1;
end
if any(strcmpi(varargin,'Plane'))
    Plane = varargin{find(strcmpi(varargin,'Plane'))+1};
else
    Plane = 'XY';
end
if snapshot_flag
    if any(strcmpi(varargin,'totalperiods'))
        totalPeriods = varargin{find(strcmpi(varargin,'totalperiods'))+1};
    else
        totalPeriods = input('totalPeriods? = ');
    end
    if any(strcmpi(varargin,'initTime'))
        initTime = varargin{find(strcmpi(varargin,'inittime'))+1};
    else
        initTime = input('initTime? = ');
    end
    timePoints = linspace(initTime,totalPeriods/fus,nrSnapshots);
end
if any(strcmpi(varargin,'scaleFactor'))
   scaleFactor = varargin{find(strcmpi(varargin,'scaleFactor'))+1};
else
    scaleFactor = 1;
end

disp(['SCALEFACTOR = ',num2str(scaleFactor)])

%change limits of plots based on model (human or mouse)
%shhm=simplifiedHumanHeadModel
if contains(lower(Filename),'human') || contains(lower(Filename),'shhm') 
    vislim = 75;
else
    vislim = 4;
    
end
if any(strcmpi(varargin,'rotate'))
    rotate_flag = varargin{find(strcmpi(varargin,'rotate'))+1};
else
    rotate_flag = input('rotate? \n');
end
quiver_flag = 0 %plot velocity quivers
%%  Step one: determine real and imaganiary part of the pressure phasor
calcReImp = 'pf'; %currently only pf based validated also 2points and allpoints possible, latter two could be more correct if not sinusoidal regime

Pflag = 0;
xaxis = Sim4life_input.Axis0;
yaxis = Sim4life_input.Axis1;
zaxis = Sim4life_input.Axis2;

if ~snapshot_flag
    calcReImp = 'pf';
else
    lastsnapshot = scaleFactor*double(Sim4life_input.(['p_x_y_z_t_Snapshot',num2str(nrSnapshots-1)]))*1e-6;
end

switch lower(calcReImp)
    case '2points'
        %double check before use not validated yet
        selSnapshots = [nrSnapshots-3,nrSnapshots-2];
        seltp = timePoints(selSnapshots);
        Pressuret1 = double(Sim4life_input.(['p_x_y_z_t_Snapshot',num2str(selSnapshots(1))]))*1e-6; %Convert to MPa;
        Pressuret2 = double(Sim4life_input.(['p_x_y_z_t_Snapshot',num2str(selSnapshots(2))]))*1e-6; %Convert to MPa;
        problem.options = optimoptions('fsolve','Display','off');
        problem.solver = 'fsolve';
        
        for iPress=1:length(Pressuret1)
            problem.x0 = [Pressuret1(iPress),0];
            problem.objective = @(p) [p(1).*cos(wus.*seltp(1)+p(2))-Pressuret1(iPress);
                p(1).*cos(wus.*seltp(2)+p(2))-Pressuret2(iPress);];
            x = fsolve(problem);
            phip(iPress) = x(2);
            Ampp(iPress) = x(1);
            Pval = ceil(iPress/length(Pressuret1)*100);
            if Pval>Pflag
                disp([num2str(Pval),'%'])
                Pflag = Pflag + 10;
            end
        end
    case 'allpoints'
        %double check before use not validated yet
        for isnap=1:nrSnapshots
            Pressurei = double(Sim4life_input.(['p_x_y_z_t_Snapshot',num2str(isnap)]))*1e-6; %Convert to MPa
            if isnap == 1
                meanPressure = 1/nrSnapshots.*Pressurei;
                maxPressure = Pressurei;
                minPressure = Pressurei;
            else
                meanPressure = meanPressure + 1/nrSnapshots.*Pressurei;
                maxPressure(Pressurei>maxPressure) = Pressurei(Pressurei>maxPressure);
                minPressure(Pressurei<minPressure) = Pressurei(Pressurei<minPressure);
            end
        end
        if nrSnapshots>=5*(timePoints(end)-timePoints(1))*fus
            Ampp = 1/2*(maxPressure-minPressure);
        else
            Ampp = abs(maxPressure);
        end
        for isnap=1:nrSnapshots
            Pressurei = double(Sim4life_input.(['p_x_y_z_t_Snapshot',num2str(isnap)]))*1e-6; %Convert to MPa
            % adapt acos result
            if isnap == 1
                phip = 1/nrSnapshots*mod(acos(Pressurei./Ampp)-wus*timePoints(isnap),2*pi);
            else
                phip = phip + 1/nrSnapshots*mod(acos(Pressurei./Ampp)-wus*timePoints(isnap),2*pi);
            end
        end
  
    case 'pf'
        Pressuref = double(Sim4life_input.(['p_x_y_z_f_Snapshot0']))*1e-6;
    otherwise
        error('wrong input calcReImp')
end

if strcmpi(calcReImp,'pf')
    realP = scaleFactor*real(Pressuref);
    imagP = scaleFactor*imag(Pressuref);
    clear('Pressuref');
else
    realP = scaleFactor*Ampp.*cos(phip);
    imagP = scaleFactor*Ampp.*sin(phip);
    clear('Ampp','phip')
end
clear('Sim4life_input')
%% step two: determine the gradient
%output of S4L contains corner locations while value calculated middle of
%cubes
xloc = (xaxis(1:end-1)+xaxis(2:end))./2*1e3;% loc in mm
yloc = (yaxis(1:end-1)+yaxis(2:end))./2*1e3;
zloc = (zaxis(1:end-1)+zaxis(2:end))./2*1e3;
RBrain = RBrain*1e3;
% reshape into 3dim matrix yaxis,xaxis,zaxis
if snapshot_flag
lastsnapshot = reshape(lastsnapshot,length(xloc),length(yloc),length(zloc)); % just for figures
lastsnapshot = permute(lastsnapshot,[2,1,3]);
lastsnapshot = lastsnapshot(1:downsf:end,1:downsf:end,1:downsf:end);
end


realPMAT = reshape(realP,length(xloc),length(yloc),length(zloc));
realPMAT = permute(realPMAT,[2,1,3]); %switch x and y 
realPMAT = realPMAT(1:downsf:end,1:downsf:end,1:downsf:end);
imagPMAT = reshape(imagP,length(xloc),length(yloc),length(zloc));
imagPMAT = permute(imagPMAT,[2,1,3]);
imagPMAT = imagPMAT(1:downsf:end,1:downsf:end,1:downsf:end);
% gradient 
[gradrealPx,gradrealPy,gradrealPz] = gradient(realPMAT); % determines gradient for unit step
[gradimagPx,gradimagPy,gradimagPz] = gradient(imagPMAT);
xloc = xloc(1:downsf:end);
yloc = yloc(1:downsf:end);
zloc = zloc(1:downsf:end);
dxloc = [diff(xloc(1:2)),xloc(3:end)-xloc(1:end-2),diff(xloc(end-1:end))]; % determine actual steps
dyloc = [diff(yloc(1:2)),yloc(3:end)-yloc(1:end-2),diff(yloc(end-1:end))];
dzloc = [diff(zloc(1:2)),zloc(3:end)-zloc(1:end-2),diff(zloc(end-1:end))];

[M,N,K] = size(gradrealPx);

gradrealPx = gradrealPx ./ repmat(dxloc,M,1,K); % rescale wrt actual steps
gradrealPy = gradrealPy ./ permute(repmat(dyloc,N,1,K),[2,1,3]);
gradrealPz = gradrealPz ./ permute(repmat(dzloc,M,1,N),[1,3,2]);
gradimagPx = gradimagPx ./ repmat(dxloc,M,1,K);
gradimagPy = gradimagPy ./ permute(repmat(dyloc,N,1,K),[2,1,3]);
gradimagPz = gradimagPz ./ permute(repmat(dzloc,M,1,N),[1,3,2]);


clear('imagP','realP')
%% step3 calculate in sin regime
vxr = complex(gradrealPx,gradimagPx)./complex(0,-wus*rhoBrain)*1e9; %pressure in MPa and displacemanet in mm
vyr = complex(gradrealPy,gradimagPy)./complex(0,-wus*rhoBrain)*1e9; % v in m/s
vzr = complex(gradrealPz,gradimagPz)./complex(0,-wus*rhoBrain)*1e9;
% just to speed up or to decrease matrix sizes 
if snapshot_flag
nelemv = numel(vxr);
parts = 10;
for i=1:parts
    if i~=10
        indices = (i-1)*floor(nelemv/parts)+1:i*floor(nelemv/parts);
    else
        indices = (i-1)*floor(nelemv/parts)+1:nelemv;
    end
    
   ampSnapshotv(indices) = vecnorm(real([vxr(indices)',vyr(indices)',vzr(indices)'].*exp(complex(0,wus*timePoints(end)))),2,2);
end
ampSnapshotv = reshape(ampSnapshotv,length(yloc),length(xloc),length(zloc));
end
clear('gradrealPx','gradrealPy','gradrealPz','gradimagPx','gradimagPy','gradimagPz','indices')
%% Determine polarization v 
angle = 1/2*atan(2*(real(vxr).*imag(vxr)+real(vyr).*imag(vyr)+real(vzr).*imag(vzr))./...
    ((real(vxr).*real(vxr)+real(vyr).*real(vyr)+real(vzr).*real(vzr))-...
    (imag(vxr).*imag(vxr)+imag(vyr).*imag(vyr)+imag(vzr).*imag(vzr)))); % transformation angle to polarization ellipse 
% new vectors
realvarx = real(vxr).*cos(angle)+imag(vxr).*sin(angle);
imagvarx = -real(vxr).*sin(angle)+imag(vxr).*cos(angle);
realvary = real(vyr).*cos(angle)+imag(vyr).*sin(angle);
imagvary = -real(vyr).*sin(angle)+imag(vyr).*cos(angle);
realvarz = real(vzr).*cos(angle)+imag(vzr).*sin(angle);
imagvarz = -real(vzr).*sin(angle)+imag(vzr).*cos(angle);
clear('vxr','vyr','vzr')
% determine longest vector and assign this as maximal amplitude
amprealvar = (realvarx.^2+realvary.^2+realvarz.^2).^(1/2);
ampimagvar = (imagvarx.^2+imagvary.^2+imagvarz.^2).^(1/2);
indicesmax = amprealvar>ampimagvar;
ellipticity = ampimagvar./amprealvar;
ellipticity(ellipticity>1) = ellipticity(ellipticity>1).^(-1);
maxvx = imagvarx;
maxvy = imagvary;
maxvz = imagvarz;
maxvx(indicesmax) = realvarx(indicesmax);
maxvy(indicesmax) = realvary(indicesmax);
maxvz(indicesmax) = realvarz(indicesmax);
maxvMAT = (maxvx.^2 + maxvy.^2 + maxvz.^2).^(1/2);
if Database_flag
       veloc.maxvx = maxvx; veloc.maxvy = maxvy;
       veloc.maxvz = maxvz;
       veloc.phasemaxdir = angle+(1-indicesmax)*pi/2;
end
% gen grid
[X,Y,Z] = meshgrid(xloc,yloc,zloc);
% store grid points in output
loc.X = X; loc.Y = Y; loc.Z = Z;
Dist_O = (X.^2+Y.^2+Z.^2).^(1/2);
idx_DistStRbrain = Dist_O<=RBrain; %search largest vibrator located in brain tissue
maxvMAT_intm = maxvMAT(idx_DistStRbrain);
maxvx_intm = maxvx(idx_DistStRbrain);
maxvy_intm = maxvy(idx_DistStRbrain);
maxvz_intm = maxvz(idx_DistStRbrain);
X_intm = X(idx_DistStRbrain); Y_intm = Y(idx_DistStRbrain);
Z_intm = Z(idx_DistStRbrain);
[maxv, maxv_idx] = max(maxvMAT_intm(:));
maxv_dir = [maxvx_intm(maxv_idx),maxvy_intm(maxv_idx),maxvz_intm(maxv_idx)];
maxv_dir = maxv_dir./norm(maxv_dir);

% store position maxv
maxv_pos = [X_intm(maxv_idx),Y_intm(maxv_idx),Z_intm(maxv_idx)];

clear('indicesmax','maxvMAT_intm','idx_DistStRbrain','Dist_O','X_intm','Y_intm','Z_intm','maxvx_intm','maxvy_intm','maxvz_intm')
% not sure if it has any use
phasevMAT = atan2(ampimagvar,amprealvar);
clear('ampimagvar','amprealvar')
% determine phase and amp of each component for Intensity calculation
phasevarx = atan2(imagvarx,realvarx);
phasevary = atan2(imagvary,realvary);
phasevarz = atan2(imagvarz,realvarz);
ampvarx = (realvarx.^2+imagvarx.^2).^(1/2);
ampvary = (realvary.^2+imagvary.^2).^(1/2);
ampvarz = (realvarz.^2+imagvarz.^2).^(1/2);
%% maximal p in time (equals modulus pf)
maxp = (realPMAT.^2+imagPMAT.^2).^(1/2);
phasep = atan2(imagPMAT,realPMAT);
% store pressure vals in output
press.maxp = maxp; press.phasep = phasep;
%% calculate intensity
Ix = 1/2*ampvarx.*maxp.*cos(angle+phasevarx-phasep)*10^6;
Iy = 1/2*ampvary.*maxp.*cos(angle+phasevary-phasep)*10^6;
Iz = 1/2*ampvarz.*maxp.*cos(angle+phasevarz-phasep)*10^6;
Inorm = (Ix.^2+Iy.^2+Iz.^2).^(1/2);
% store velocitiy vectors into output
veloc.ampx = ampvarx; veloc.phasex = mod(phasevarx+angle,pi); % correction for transformation of velocity vectors
veloc.ampy = ampvary; veloc.phasey = mod(phasevary+angle,pi); % correction for transformation of velocity vectors
veloc.ampz = ampvarz; veloc.phasez = mod(phasevarz+angle,pi); % correction for transformation of velocity vectors


clear('ampvarx','ampvary','ampvarz','phasevarx','phasevary','phasevarz')
clear('Ix','Iy','Iz')
%%

zloc0s = 0;
if Plot_flag
%% maximal pressure and phase
Ysurfaceindices = yloc>-vislim &yloc<vislim;
Permuted_flag = 0;
for i=1:length(zloc0s)
    if strcmpi(Plane,'YZ') && ~Permuted_flag
        % code is programmed for XY plane so switch X and Z for YZ plane
        % (ZY to be more correct)
        xloc0 = zloc0s(i);
        idxz0 = find(abs(xloc-xloc0)==min(abs(xloc-xloc0))); idxz0 = idxz0(1);
        X = permute(X,[1,3,2]);
        Y = permute(Y,[1,3,2]);
        Z = permute(Z,[1,3,2]);
        maxp = permute(maxp,[1,3,2]);
        phasep = permute(phasep,[1,3,2]);
        if snapshot_flag
        ampSnapshotv = permute(ampSnapshotv,[1,3,2]);
        lastsnapshot = permute(lastsnapshot,[1,3,2]);
        end
        maxvMAT = permute(maxvMAT,[1,3,2]);
        ellipticity = permute(ellipticity, [1,3,2]);
        Permuted_flag = 1;
        xloc_interm = xloc;
        xloc = zloc;
        zloc = xloc_interm;
        X_interm = X; X = Z; Z = X_interm;
        Inorm = permute(Inorm, [1,3,2]);
    else
        zloc0 = zloc0s(i);
        idxz0 = find(abs(zloc-zloc0)==min(abs(zloc-zloc0))); idxz0 = idxz0(1);
    end
    xindices = find(xloc>-2 & xloc<2); xindices = xindices(floor(linspace(1,length(xindices),xres)));
    yindices = find(yloc>-6 & yloc<6); yindices = yindices(floor(linspace(1,length(yindices),yres)));
    [XQ,YQ] = meshgrid(xloc(xindices),yloc(yindices));
    %ZQ = PMAT(yindices,xindices,idxz0);
    ZQ = ones(size(XQ));
    UQ = maxvx(yindices,xindices); VQ = maxvy(yindices,xindices);
    WQ = maxvz(yindices,xindices);
    if snapshot_flag
        subplot(2,3,1)
        surf(X(:,:,idxz0),Y(:,:,idxz0),zloc(idxz0)*ones(size(X(:,:,idxz0))),lastsnapshot(:,:,idxz0),'edgecolor','none','facealpha',0.75)
        c = colorbar;
        c.Label.String = 'Pressure [MPa]';
        c.Label.FontSize = 14;
        colormap('hot');
        view(90*(1-double(rotate_flag)),90)
        axis('square')
        xlabel({'x axis [mm]'})
        ylabel({'Pressure','','y axis [mm]'})
        xlims = get(gca,'xlim');
        set(gca,'zlim',xlims);
        title('snapshot')
    end
    subplot(2,2+snapshot_flag,1+snapshot_flag)
    hold on
    surf(X(:,:,idxz0),Y(:,:,idxz0),zloc(idxz0)*ones(size(X(:,:,idxz0))),maxp(:,:,idxz0),'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'pressure [MPa]';
    c.Label.FontSize = 14;
    colormap('hot');
    %quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'zlim',xlims);
    title('max amplitude')
    
    subplot(2,2+snapshot_flag,2+snapshot_flag)
    hold on
    surf(X(:,:,idxz0),Y(:,:,idxz0),zloc(idxz0)*ones(size(X(:,:,idxz0))),phasep(:,:,idxz0),'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'phase';
    c.Label.FontSize = 14;
    c.Ticks = [-pi:pi/2:pi];
    c.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
    colormap('hot');
    %quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'zlim',xlims);
    title('phase')
    
    % velocity plots
    if snapshot_flag
        subplot(2,3,4)
        surf(X(:,:,idxz0),Y(:,:,idxz0),zloc(idxz0)*ones(size(X(:,:,idxz0))),ampSnapshotv(:,:,idxz0),'edgecolor','none','facealpha',0.75)
        c = colorbar;
        c.Label.String = 'velocity [m/s]';
        c.Label.FontSize = 14;
        colormap('hot');
        view(90*(1-double(rotate_flag)),90)
        axis('square')
        xlabel('x axis [mm]')
        ylabel({'Velocity','','y axis [mm]'})
        xlims = get(gca,'xlim');
        set(gca,'zlim',xlims);
        %title('snapshot')
    end
    
    subplot(2,2+snapshot_flag,3+2*snapshot_flag)
    hold on
    surf(X(:,:,idxz0),Y(:,:,idxz0),zloc(idxz0)*ones(size(X(:,:,idxz0))),maxvMAT(:,:,idxz0),'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'velocity [m/s]';
    c.Label.FontSize = 14;
    colormap('hot');
    %quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'zlim',xlims);
    %title('amplitude')
    
    subplot(2,2+snapshot_flag,4+2*snapshot_flag)
    hold on
    surf(X(Ysurfaceindices,:,idxz0),Y(Ysurfaceindices,:,idxz0),zloc(idxz0)*ones(size(X(Ysurfaceindices,:,idxz0))),ellipticity(Ysurfaceindices,:,idxz0),'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'ellipticity';
    c.Label.FontSize = 14;
    caxis('auto')
    %c.Ticks = [-pi:pi/2:pi];
    %c.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
    colormap('hot');
    %quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'zlim',xlims);
    title('ellipticity')
    set(findobj('type','axes'),'fontsize',14)
end
%% maximal pressure and velocity and displacement
figure
Ysurfaceindices = yloc>-vislim &yloc<vislim;
for i=1:length(zloc0s)
    if strcmpi(Plane,'YZ') && ~Permuted_flag
        % code is programmed for XY plane so switch X and Z for YZ plane
        % (ZY to be more correct)
        xloc0 = zloc0s(i);
        idxz0 = find(abs(xloc-xloc0)==min(abs(xloc-xloc0))); idxz0 = idxz0(1);
        X = permute(X,[1,3,2]);
        Y = permute(Y,[1,3,2]);
        Z = permute(Z,[1,3,2]);
        if snapshot_flag
        lastsnapshot = permute(lastsnapshot,[1,3,2]);
        ampSnapshotv = permute(ampSnapshotv,[1,3,2]);
        end
        maxp = permute(maxp,[1,3,2]);
        phasep = permute(phasep,[1,3,2]);
        maxvMAT = permute(maxvMAT,[1,3,2]);
        ellipticity = permute(ellipticity, [1,3,2]);
        Inorm = permute(Inorm, [1,3,2]);
    else
        zloc0 = zloc0s(i);
        idxz0 = find(abs(zloc-zloc0)==min(abs(zloc-zloc0))); idxz0 = idxz0(1);
    end
    xindices = find(xloc>-2 & xloc<2); xindices = xindices(floor(linspace(1,length(xindices),xres)));
    yindices = find(yloc>-6 & yloc<6); yindices = yindices(floor(linspace(1,length(yindices),yres)));
    [XQ,YQ] = meshgrid(xloc(xindices),yloc(yindices));
    %ZQ = PMAT(yindices,xindices,idxz0);
    ZQ = ones(size(XQ));
    UQ = maxvx(yindices,xindices); VQ = maxvy(yindices,xindices);
    WQ = maxvz(yindices,xindices);
    subplot(2,3,1)
    hold on
    surf(X(Ysurfaceindices,:,idxz0),Y(Ysurfaceindices,:,idxz0),...
        zloc(idxz0)*ones(size(X(Ysurfaceindices,:,idxz0))),maxp(Ysurfaceindices,:,idxz0),'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'Amplitude [MPa]';
    c.Label.FontSize = 14;
    colormap('hot');
    if quiver_flag
        quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    end
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'zlim',xlims);
    title('pressure amplitude')
    subplot(2,3,2)
        hold on
    surf(X(Ysurfaceindices,:,idxz0),Y(Ysurfaceindices,:,idxz0),...
        zloc(idxz0)*ones(size(X(Ysurfaceindices,:,idxz0))),maxvMAT(Ysurfaceindices,:,idxz0),'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'velocity [m/s]';
    c.Label.FontSize = 14;
    colormap('hot');
    if quiver_flag
        quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    end
    scatter3(maxv_pos(1),maxv_pos(2),maxv_pos(3),'cX')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'zlim',xlims);
    title('velocity amplitude')
    subplot(2,3,3)
    hold on
    surf(X(Ysurfaceindices,:,idxz0),Y(Ysurfaceindices,:,idxz0),...
        zloc(idxz0)*ones(size(X(Ysurfaceindices,:,idxz0))),maxvMAT(Ysurfaceindices,:,idxz0)/wus*1e9,'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'Displacement [nm]';
    c.Label.FontSize = 14;
    colormap('hot');
    if quiver_flag
        quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    end
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'zlim',xlims);
    title('displacement amplitude')
    if ~rotate_flag
    subplot(2,3,[4,5,6])
    xloc0 = 0;
    idxx0 = find(abs(xloc-xloc0)==min(abs(xloc-xloc0)));
    yyaxis left
    plot(Y(Ysurfaceindices,idxx0,idxz0),maxp(Ysurfaceindices,idxx0,idxz0))
    ylabel('Pressure Amplitude [MPa]')
    xlabel('y axis [mm]')
    yyaxis right
    plot(Y(Ysurfaceindices,idxx0,idxz0),maxvMAT(Ysurfaceindices,idxx0,idxz0)/wus*1e9)
    ylabel('Displacement Amplitude [nm]')
    title('x \approx 0 mm, z \approx 0 mm')
    else
    subplot(2,3,[4,5,6])
    yloc0 = 0;
    Xsurfaceindices = xloc>-vislim &xloc<vislim;
    idxx0 = find(abs(yloc-yloc0)==min(abs(yloc-yloc0)));
    yyaxis left
    plot(X(idxx0(1),Xsurfaceindices,idxz0(1)),maxp(idxx0(1),Xsurfaceindices,idxz0(1)))
    ylabel('Pressure Amplitude [MPa]')
    xlabel('y axis [mm]')
    yyaxis right
    plot(X(idxx0(1),Xsurfaceindices,idxz0(1)),maxvMAT(idxx0(1),Xsurfaceindices,idxz0(1))/wus*1e9)
    ylabel('Displacement Amplitude [nm]')
    title('x \approx 0 mm, z \approx 0 mm')
    end
end

%% maximal pressure and velocity and intensity
figure
Ysurfaceindices = yloc>-vislim &yloc<vislim;
Xsurfaceindices = xloc>-vislim &xloc<vislim;
for i=1:length(zloc0s)
    zloc0 = zloc0s(i);
    idxz0 = find(abs(zloc-zloc0)==min(abs(zloc-zloc0))); idxz0 = idxz0(1);
    xindices = find(xloc>-2 & xloc<2); xindices = xindices(floor(linspace(1,length(xindices),xres)));
    yindices = find(yloc>-6 & yloc<6); yindices = yindices(floor(linspace(1,length(yindices),yres)));
    [XQ,YQ] = meshgrid(xloc(xindices),yloc(yindices));
    %ZQ = PMAT(yindices,xindices,idxz0);
    ZQ = ones(size(XQ));
    UQ = maxvx(yindices,xindices); VQ = maxvy(yindices,xindices);
    WQ = maxvz(yindices,xindices);
    subplot(1,3,1)
    hold on
    surf(X(Ysurfaceindices,Xsurfaceindices,idxz0),Y(Ysurfaceindices,Xsurfaceindices,idxz0),...
        zloc(idxz0)*ones(size(X(Ysurfaceindices,Xsurfaceindices,idxz0))),maxp(Ysurfaceindices,Xsurfaceindices,idxz0),'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'Amplitude [MPa]';
    c.Label.FontSize = 16;
    colormap('hot');
    %quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    set(gca,'xlim',[-vislim,vislim],'ylim',[-vislim,vislim],'zlim',[-vislim,vislim])
    title('pressure amplitude')
    subplot(1,3,2)
        hold on
    surf(X(Ysurfaceindices,Xsurfaceindices,idxz0),Y(Ysurfaceindices,Xsurfaceindices,idxz0),...
        zloc(idxz0)*ones(size(X(Ysurfaceindices,Xsurfaceindices,idxz0))),maxvMAT(Ysurfaceindices,Xsurfaceindices,idxz0)/wus*1e9,'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'Displacement [nm]';
    c.Label.FontSize = 16;
    colormap('hot');
    %quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    set(gca,'xlim',[-vislim,vislim],'ylim',[-vislim,vislim],'zlim',[-vislim,vislim])
    title('displacement amplitude')
    subplot(1,3,3)
    hold on
    surf(X(Ysurfaceindices,Xsurfaceindices,idxz0),Y(Ysurfaceindices,Xsurfaceindices,idxz0),...
        zloc(idxz0)*ones(size(X(Ysurfaceindices,Xsurfaceindices,idxz0))),Inorm(Ysurfaceindices,Xsurfaceindices,idxz0)*1e-4,'edgecolor','none','facealpha',0.75)
    c = colorbar;
    c.Label.String = 'Intensity [W/cm^2]';
    c.Label.FontSize = 16;
    colormap('hot');
    %quiver3(XQ,YQ,zloc(idxz0)*ones(size(XQ)),UQ,VQ,WQ,'b')
    hold off
    view(90*(1-double(rotate_flag)),90)
    axis('square')
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    xlims = get(gca,'xlim');
    set(gca,'xlim',[-vislim,vislim],'ylim',[-vislim,vislim],'zlim',[-vislim,vislim])
    
    title('Intensity amplitude')
    
end

set(findobj('type','axes'),'fontsize',16)
end
end


