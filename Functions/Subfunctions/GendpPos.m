function [CSource,CSink] = GendpPos(type,OrienType,d,RSphere,N,Display,varargin)
% declare Csource and Csink of N dipoles according to distribution type
% 'type' and Orientation 'OrienType'.
%other Inputs
%       d: distance CSource and CSink
%       RSphere: Max radial position
%       N number of dipoles to distribute
%       display progress
%       varargin: variable input according to string value pair 
%                  -nLayers for fibonaccy multiple
%                  -'Resolution' on concentricsingle
%                  -'S4l'Input
if Display
fprintf('SETTINGS Dipoles:\n')
fprintf('\t Dipole distribution type: %s\n',type);
fprintf('\t Dipole orientation type: %s\n',OrienType);
fprintf('\t Distance CSource-CSink: %3.2g\n',d);
fprintf('\t Maximum R value: %3.2g\n',RSphere);
end

% Create centrum dipole locations
switch lower(type)
    case 'random1'
        %random distribution along sphere
        rvals = 2*rand(N,1)-1;
        elevation = asin(rvals);
        azimuth = 2*pi*rand(N,1);
        radii = RSphere*(rand(N,1).^(1/3));
        [x,y,z] = sph2cart(azimuth,elevation,radii);
        Dipoles = [x,y,z];
        fprintf('\t Number of Dipoles: %d\n',N);
    case 'random2'
        %randm distribution along sphere (similar to 1)
        u = rand(N,1);
        v = rand(N,1);
        theta = 2*pi*u;
        phi = acos(2*v-1);
        radii = RSphere*(rand(N,1)).^(1/3);
        Dipoles = [radii.*cos(theta).*sin(phi),radii.*sin(theta).*sin(phi),radii.*cos(phi)];
        fprintf('\t Number of Dipoles: %d\n',N);
    case 'random3'
        %random distribution along sphere special case! source?
        X = randn(N,3);
        s2 = sum(X.^2,2);
        Dipoles = X.*repmat(RSphere*(gammainc(s2/2,3/2).^(1/3))./sqrt(s2),1,3);
        fprintf('\t Number of Dipoles: %d\n',N);
    case 'concentricsingle'
        % dipoles on a sphere but not uniform (like a rose)
        
        Rd = RSphere; % concentric sphere on which dipoles are situated
        Resolution = 10;
        if mod(length(varargin)+1,2)
            if ~isempty(varargin)
                if any(strcmpi(varargin,'Resolution'))
                    Resolution = varargin{find(strcmpi(varargin,'Resolution'))+1};
                end
            end
        end        
        dTheta = pi/Resolution;
        dPhi0 = pi/Resolution;
        idx = 0;
        Dipoles=[];
        for i=[0:dTheta:pi]
            if i==0 || i==pi
                Angles = 0;
            else
                dPhi = dPhi0;%*(1+Resolution/4*4/pi^2*(pi/2-i)^2);
                Angles = [dPhi:dPhi:2*pi];
            end
            for j= 1:length(Angles)
                idx=idx+1;
                Dipoles(idx,:) = Rd*[cos(Angles(j))*sin(i),sin(Angles(j))*sin(i),cos(i)];
            end
        end
        fprintf('\t Number of Dipoles: %d\n',size(Dipoles,1));
        fprintf('\t Resolution: %d\n',Resolution);
    case 'fibonaccisingle'
        %(nearly)unifrom distribution along certain sphere Lookup for math
        %behind it
        indices = [0:N-1]+0.5;
        theta = acos(1-2*indices'/N);
        phi = pi*(1+5^(1/2))*indices';
        Dipoles = RSphere*[cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta)];
    case 'fibonaccimultiple'
        %muliple layers of uniform distribtuions. distribute total # of
        %dipoles (N) according to surface size each sphere
        nLayers = 10;
        if mod(length(varargin)+1,2)
            if ~isempty(varargin)
                if any(strcmpi(varargin,'nLayers'))
                    nLayers = varargin{find(strcmpi(varargin,'nLayers'))+1};
                end
            end
        end
        Rvals = linspace(0,RSphere,nLayers+1); Rvals = Rvals(2:end);
        surfaceRvals = 4*pi*Rvals.^2;
        fracN = floor(surfaceRvals/sum(surfaceRvals)*N);
        fracN(end) = N-sum(fracN(1:end-1));
        for i = 1:length(Rvals)
        indices = [0:fracN(i)-1]+0.5;
        theta = acos(1-2*indices'/fracN(i));
        phi = pi*(1+5^(1/2))*indices';
        if i == 1
            Dipoles = Rvals(i)*[cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta)];
        else
            Dipoles = vertcat(Dipoles,Rvals(i)*[cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta)]);
        end
        end
    case 's4l_hotspot'
        %Determine dipole positions based on maximal values of simulations
        %(hotspots=position high vibration) vicinity is limited by
        %resolution

        if any(strcmpi(varargin,'S4l_Input'))
            S4l_Input = varargin{find(strcmpi(varargin,'S4l_Input'))+1};
        else
            error('no S4l_data provided')
        end
        if isfield(S4l_Input,'resolution')
            Res = S4l_Input.resolution;
            if numel(Res)==1
                exclMethod = 'sphere';
            fprintf('\ninput resolution of %3.2f mm is used\n',Res*1e3)
            else
                exclMethod = 'cuboid';
            fprintf('\ninput resolution of [x,y,z]~[%3.2f,%3.2f,%3.2f] mm is used\n',Res*1e3)
            end
        else
            Res = 0.003;
            fprintf('\nstandard resolution of %3.2f mm is used\n',Res*1e3)
        end
        
        switch S4l_Input.loc.unit
            case 'm'
                f_loc = 1;
            case 'mm'
                f_loc = 1e-3;
            otherwise
                error('unit not implemented')
        end
        
        S4l_Input.loc.X = S4l_Input.loc.X*f_loc;
        S4l_Input.loc.Y = S4l_Input.loc.Y*f_loc;
        S4l_Input.loc.Z = S4l_Input.loc.Z*f_loc;
        %     S4l_maxX = max(S4l_Input.loc.X(:)); S4l_minX = min(S4l_Input.loc.X(:));
        %     S4l_maxY = max(S4l_Input.loc.Y(:)); S4l_minY = min(S4l_Input.loc.Y(:));
        %     S4l_maxZ = max(S4l_Input.loc.Z(:)); S4l_minZ = min(S4l_Input.loc.Z(:));
        %     S4l_maxloc = [S4l_maxX, S4l_maxY, S4l_maxZ];
        %     S4l_minloc = [S4l_minX, S4l_minY, S4l_minZ];
        
        locx = S4l_Input.loc.X(:);
        locy = S4l_Input.loc.Y(:);
        locz = S4l_Input.loc.Z(:);
        vamp = (S4l_Input.veloc.ampx(:).^2+S4l_Input.veloc.ampy(:).^2+S4l_Input.veloc.ampz(:).^2).^(1/2);
        Dp_pos_avail = ones(size(vamp));
        Dp_pos_select = zeros(size(vamp));
        idx_zero = ((locx.^2+locy.^2+locz.^2).^(1/2))>RSphere; %exclode positions not inside brain
        Dp_pos_avail(idx_zero) = 0;
        
        iter = 0;
        itermax = N*2;
        Pflag = 0;
        
        
        if sum(Dp_pos_avail)<N
            dpsearch_flag = 0;
            Dp_pos_select(Dp_pos_avail) = 1;
            fprintf('\n less than requested dipole positions available\n')
        else
            dpsearch_flag = 1;
        end
        while dpsearch_flag && any(Dp_pos_avail) && iter<=itermax
            iter=iter+1;
            vamp(~Dp_pos_avail) = 0; %not available set vamp = 0=> will not be considered as hotspot
            idx_max_avail = find(vamp==max(vamp)); 
            if iter==1
                idx_absmax = idx_max_avail(1);
            end
            while ~isempty(idx_max_avail)
                %multiple positions with same amplitude possible
                imaxs = 1;
                %check if current max is avail and not yet selected
                if Dp_pos_avail(idx_max_avail(imaxs)) && ~Dp_pos_select(idx_max_avail(imaxs))
                    Dp_pos_select(idx_max_avail(imaxs)) = 1; % found maximum turn into dipole
                    xpos_sel = locx(idx_max_avail(imaxs)); ypos_sel = locy(idx_max_avail(imaxs));
                    zpos_sel = locz(idx_max_avail(imaxs));
                    switch exclMethod
                        case 'sphere'
                            dist_otherp = ((locx-xpos_sel).^2+(locy-ypos_sel).^2+(locz-zpos_sel).^2).^(1/2); %distance other dps
                            idx_cl_res = dist_otherp<Res; %idx of positions closer than resolution
                        case 'cuboid'
                            idx_cl_resx = abs((locx-xpos_sel))<Res(1);
                            idx_cl_resy = abs((locy-ypos_sel))<Res(2);
                            idx_cl_resz = abs((locz-zpos_sel))<Res(3);
                            idx_cl_res = idx_cl_resx & idx_cl_resy & idx_cl_resz;
                            
                        otherwise
                            error('method not encoded or wrong');
                    end
                    Dp_pos_avail(idx_cl_res)=0; %set positions closer than resolution to zero idx_cl_res is array siwe of Dp_pos_avail with 0 and 1, Dp_posAvail start all ones but changes to zeros, keeps having same size
                    idx_max_avail(imaxs)=[]; %position is added so delete from search
                    % delete other max that where removed from Dp_pos_avail
                    % this iteration. idx_max_avail, small array with few
                    % numbers(positions) => create array of 0and1 size
                    % idx_cl_res if overlap should be excluded  % '(maxlist_temp&idx_cl_res)'
                    maxlist_temp = zeros(size(idx_cl_res)); maxlist_temp(idx_max_avail)=1;
                    idx_max_avail = find(maxlist_temp-(maxlist_temp&idx_cl_res));
                    
                else
                    error('problem with algorithm, position not avail or already selected')
                end
                
                if sum(Dp_pos_select)==N
                    break
                elseif sum(Dp_pos_select)>N
                    error('too much dipoles selected')
                end
                
            end
            if Display
                Pval = round(sum(Dp_pos_select)/N*100,0);
                if Pval >= Pflag
                    disp(['Dp position selection: ',num2str(Pval),'%'])
                    Pflag = Pflag+5;
                end
            end
            dpsearch_flag = sum(Dp_pos_select)<N;
            
        end
        Dp_pos_select(idx_absmax) = 0;
        %DOI is by default first dipole in this case biggest oscillator
        Dipoles = [locx(idx_absmax),locy(idx_absmax),locz(idx_absmax);
            locx(logical(Dp_pos_select)),locy(logical(Dp_pos_select)),locz(logical(Dp_pos_select))];
    otherwise
        error('false input')
end
% Create CSource and CSink positions
switch lower(OrienType)
    case 'radial'
        CSource = (1+d/2./vecnorm(Dipoles,2,2)).*Dipoles;
        CSink = (1-d/2./vecnorm(Dipoles,2,2)).*Dipoles;
    case 'random'
        Orientation =  rand(size(Dipoles,1),3)*2-1; %rand numbers between [-1 and 1]
        Orientation = Orientation./vecnorm(Orientation,2,2); %normalize
        CSource = Dipoles+d/2*Orientation;
        CSink = Dipoles-d/2*Orientation;
    otherwise
        error('valse orien type input')
end
        
end