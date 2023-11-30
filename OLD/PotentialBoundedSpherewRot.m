function [VRMAT,VRMAT2] = PotentialBoundedSpherewRot(Inputdipole,dipoleI,sigma,xsphere,ysphere,zsphere,Inputtype,SolutionType)

ShowRotCorrect = 1;
RSphere = norm([xsphere(1,1),ysphere(1,1),zsphere(1,1)]);

switch Inputtype
    case 'Euler'
    Posdipole = Inputdipole.Posdipole;
    dipoleTheta = Inputdipole.dipoleTheta;
    dipolePhi = Inputdipole.dipolePhi;
    dipoled = Inputdipole.dipoled;
 
    CSource = Posdipole+dipoled*[cos(dipolePhi).*sin(dipoleTheta),sin(dipolePhi).*sin(dipoleTheta),cos(dipoleTheta)];
    CSink = Posdipole-dipoled*[cos(dipolePhi).*sin(dipoleTheta),sin(dipolePhi).*sin(dipoleTheta),cos(dipoleTheta)];
    if norm(CSource,2) >= RSphere || norm(CSink,2) >= RSphere
        error('dipole not in sphere')
    end
    case 'Spherical'
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
if any(CSource(1:2))
    unitCsource = CSource/norm(CSource);         % calculate unit vector PoscurrentSource
    v = cross(unitCsource,[0,0,1]);
    s = norm(v);
    c = dot(unitCsource,[0,0,1]);
    % rotation matrix transfrom CSource to z-axis
    RCSource = eye(3) + [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0]+[0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0]^2.*(1-c)/s^2;
    
else
    RCSource = eye(3);
end
CSourceInit = CSource; 
CSinkInit = CSink;
CSourceInter = (RCSource*CSourceInit')';
CSinkInter = (RCSource*CSinkInit')';

if CSink(2)~=0    
    rotphi = -(atan(CSinkInter (2)./CSinkInter(1))+double(CSinkInter(1)<0)*pi);
    RCSink = [cos(rotphi) -sin(rotphi) 0;sin(rotphi) cos(rotphi) 0;0 0 1];
else
    RCSink = eye(3);
end
Rtot = RCSink*RCSource;
CSource = (Rtot*CSourceInit')';   
CSink = (Rtot*CSinkInit')';
   
if ShowRotCorrect
    figure
    subplot(1,2,1)
    plot3([0,CSourceInit(1)],[0,CSourceInit(2)],[0,CSourceInit(3)],'b--','LineWidth',2)
    hold on
    plot3([0,CSinkInit(1)],[0,CSinkInit(2)],[0,CSinkInit(3)],'r--','LineWidth',2)
    plot3([0,CSourceInter(1)],[0,CSourceInter(2)],[0,CSourceInter(3)],'b-.')
    plot3([0,CSinkInter(1)],[0,CSinkInter(2)],[0,CSinkInter(3)],'r-.')
    plot3([0,CSource(1)],[0,CSource(2)],[0,CSource(3)],'b-')
    plot3([0,CSink(1)],[0,CSink(2)],[0,CSink(3)],'r-')
    %surf(xs,ys,zs,'FaceAlpha',0.1,'EdgeColor','none','SpecularStrength',0);
    axis equal
    hold off
    view(120,30)
    hAxis=gca;
    hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
    hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
    hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
    hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
    xlabel('xaxis')
    ylabel('yaxis')
    zlabel('zaxis')
end

switch SolutionType
    case 'ClosedBounded'
        ra = @(x,y,z,CSink) sqrt((x-CSink(1)).^2+(y-CSink(2)).^2+(z-CSink(3)).^2);
        rb = @(x,y,z,CSource) sqrt((x-CSource(1)).^2+(y-CSource(2)).^2+(z-CSource(3)).^2);
        cosB = @(x,y,z,CSink) (ra(x,y,z,CSink).^2-(x.^2+y.^2+z.^2)-norm(CSink).^2)./(-2*sqrt((x.^2+y.^2+z.^2))*norm(CSink));
        cosT = @(x,y,z,CSource) (rb(x,y,z,CSource).^2-(x.^2+y.^2+z.^2)-norm(CSource).^2)./(-2*sqrt((x.^2+y.^2+z.^2))*norm(CSource));
        cosT2 = @(x,y,z,CSource) (z./sqrt(x.^2+y.^2+z.^2));
        NewCoorSphere = (Rtot*[xsphere(:),ysphere(:),zsphere(:)]')';
        xspherenew = reshape(NewCoorSphere(:,1),size(xsphere));
        yspherenew = reshape(NewCoorSphere(:,2),size(xsphere));
        zspherenew = reshape(NewCoorSphere(:,3),size(xsphere));        
        
        VR = @(x,y,z,CSource,CSink,cosT) dipoleI/(4*pi*sigma)*(2./(rb(x,y,z,CSource))-2./(ra(x,y,z,CSink))+...
            1/RSphere.*log((ra(x,y,z,CSink)+RSphere-norm(CSink).*cosB(x,y,z,CSink))./(rb(x,y,z,CSource)+RSphere-norm(CSource).*cosT(x,y,z,CSource))));
        VRMAT = VR(xspherenew(:),yspherenew(:),zspherenew(:),CSource,CSink,cosT2);
        VRMAT2 = VR(xsphere(:),ysphere(:),zsphere(:),CSourceInit,CSinkInit,cosT);
        VRMAT = reshape(VRMAT,size(xsphere));
        VRMAT2 = reshape(VRMAT2,size(xsphere));
        
    otherwise
        disp('to be continued')
        VRMAT = 0;
end
if ShowRotCorrect
    figure
    subplot(1,3,1)
    arrow3d([CSinkInit(1),CSourceInit(1)],[CSinkInit(2),CSourceInit(2)],[CSinkInit(3),CSourceInit(3)],0.75,RSphere/100)
    hold on
    surf(xsphere,ysphere,zsphere,VRMAT2,'FaceAlpha',0.5,'EdgeColor','none','SpecularStrength',0);
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
    subplot(1,3,2)
    arrow3d([CSinkInit(1),CSourceInit(1)],[CSinkInit(2),CSourceInit(2)],[CSinkInit(3),CSourceInit(3)],0.75,RSphere/100)
    hold on    
    surf(xsphere,ysphere,zsphere,VRMAT,'FaceAlpha',0.5,'EdgeColor','none','SpecularStrength',0);
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
    subplot(1,3,3)
    arrow3d([CSink(1),CSource(1)],[CSink(2),CSource(2)],[CSink(3),CSource(3)],0.75,RSphere/100)
    hold on    
    surf(xspherenew,yspherenew,zspherenew,VRMAT,'FaceAlpha',0.5,'EdgeColor','none','SpecularStrength',0);
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
end

end
