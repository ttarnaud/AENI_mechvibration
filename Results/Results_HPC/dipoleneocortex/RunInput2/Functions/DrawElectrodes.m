function DrawElectrodes(POIs,varargin)

axes(gca)
hold on

R = 2*10^-3;
L = 1*10^-3;
Color = 'k';
if mod(length(varargin)+1,2)
    if ~isempty(varargin)
        if any(strcmpi(varargin,'axes'))
             axes(varargin{find(strcmpi(varargin,'axes'))+1});
        end
        if any(strcmpi(varargin,'length'))
            L = varargin{find(strcmpi(varargin,'Size'))+1};
        end
        if any(strcmpi(varargin,'radius'))
            R = varargin{find(strcmpi(varargin,'Size'))+1};
        end
    end
else
    error('incorrect input')
end



for i=1:size(POIs,1)
% calculating rotation matrix and rotate
[X,Y,Z] = cylinder(R,100);
Z(2,:) = L.*Z(2,:);
RPOI = norm(POIs(i,:));
nPOI = [POIs(i,:)]/RPOI;
v = cross([0,0,1],nPOI);
s = norm(v);
c = dot([0,0,1],nPOI);
skewv = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
RM = eye(3)+skewv+skewv^2*(1-c)/(s^2+double(s==0));
Cyl1 = [X(1,:);Y(1,:);Z(1,:)];
Cyl2 = [X(2,:);Y(2,:);Z(2,:)];
Cyl1 = RM*Cyl1+POIs(i,:)';
Cyl2 = RM*Cyl2+POIs(i,:)';
X = [Cyl1(1,:);Cyl2(1,:)]; Y = [Cyl1(2,:);Cyl2(2,:)]; Z = [Cyl1(3,:);Cyl2(3,:)];

% Plot and fill
surf(X,Y,Z,'facecolor',Color,'LineStyle','none');
fill3(X(1,:),Y(1,:),Z(1,:),Color)
fill3(X(2,:),Y(2,:),Z(2,:),Color)
text((RPOI+8*L)*nPOI(1),(RPOI+8*L)*nPOI(2),(RPOI+8*L)*nPOI(3),num2str(i),'FontSize',10,'Color','r');
end

hold off
end