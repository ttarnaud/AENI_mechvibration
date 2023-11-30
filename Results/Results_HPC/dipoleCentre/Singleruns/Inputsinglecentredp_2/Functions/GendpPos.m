function [CSource,CSink] = GendpPos(type,OrienType,d,RSphere,N,varargin)
fprintf('SETTINGS Dipoles:\n')
fprintf('\t Dipole distribution type: %s\n',type);
fprintf('\t Dipole orientation type: %s\n',OrienType);
fprintf('\t Distance CSource-CSink: %3.2g\n',d);
fprintf('\t Maximum R value: %3.2g\n',RSphere);
% Create centrum dipole locations
switch lower(type)
    case 'random1'
        rvals = 2*rand(N,1)-1;
        elevation = asin(rvals);
        azimuth = 2*pi*rand(N,1);
        radii = RSphere*(rand(N,1).^(1/3));
        [x,y,z] = sph2cart(azimuth,elevation,radii);
        Dipoles = [x,y,z];
        fprintf('\t Number of Dipoles: %d\n',N);
    case 'random2'
        u = rand(N,1);
        v = rand(N,1);
        theta = 2*pi*u;
        phi = acos(2*v-1);
        radii = RSphere*(rand(N,1)).^(1/3);
        Dipoles = [radii.*cos(theta).*sin(phi),radii.*sin(theta).*sin(phi),radii.*cos(phi)];
        fprintf('\t Number of Dipoles: %d\n',N);
    case 'random3'
        X = randn(N,3);
        s2 = sum(X.^2,2);
        Dipoles = X.*repmat(RSphere*(gammainc(s2/2,3/2).^(1/3))./sqrt(s2),1,3);
        fprintf('\t Number of Dipoles: %d\n',N);
    case 'concentricsingle'
        % dipoles on a sphere
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
        indices = [0:N-1]+0.5;
        theta = acos(1-2*indices'/N);
        phi = pi*(1+5^(1/2))*indices';
        Dipoles = RSphere*[cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta)];
    case 'fibonaccimultiple'
        nLayers = 4;
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
    otherwise
        error('false input')
end
% Create CSource and CSink positions
switch lower(OrienType)
    case 'radial'
        CSource = (1+d/2./vecnorm(Dipoles,2,2)).*Dipoles;
        CSink = (1-d/2./vecnorm(Dipoles,2,2)).*Dipoles;
    case 'random'
        Orientation =  rand(size(Dipoles,1),3)*2-1;
        Orientation = Orientation./vecnorm(Orientation,2,2);
        CSource = Dipoles+d/2*Orientation;
        CSink = Dipoles-d/2*Orientation;
    otherwise
        error('valse orien type input')
end
        
end