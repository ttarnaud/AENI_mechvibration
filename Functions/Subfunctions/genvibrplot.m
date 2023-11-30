function genvibrplot(Csource,Csink,dAmpx,dAmpy,dAmpz,phasex,phasey,phasez,nr)
% see course of electromagnetism: here we transform the directions to the
% dominant and recessive of polarization ellips. this is done prior to determination
% of ellepticity. a plot is created with polarization ellips displayed as
% wel as both dominant and recessive directions and position of rotor at
% t=0
x = Csink(1);
y = Csink(2);
z = Csink(3);
u = Csource(1)-Csink(1);
v = Csource(2)-Csink(2);
w = Csource(3)-Csink(3);

x2 = (Csource(1)+Csink(1))/2;
y2 = (Csource(2)+Csink(2))/2;
z2 = (Csource(3)+Csink(3))/2;

vxr = dAmpx.*exp(complex(0,phasex-pi/2)); % -pi/2 because phase is of sin not cosine
vyr = dAmpy.*exp(complex(0,phasey-pi/2));
vzr = dAmpz.*exp(complex(0,phasez-pi/2));

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

realvec = [realvarx,realvary,realvarz];
amprealvec = norm(realvec);
if amprealvec == 0 % check not zero avoid error trown
    valforplot = realvec;
else
    valforplot = realvec/amprealvec;
end
strforplot = sprintf('%5.2e, ',valforplot);
strrealvec = ['real vector: ',sprintf('%5.2e',amprealvec),'\cdot [',strforplot(1:end-2),']'];
imagvec = [imagvarx,imagvary,imagvarz];
ampimagvec = norm(imagvec);
if ampimagvec == 0
    valforplot = imagvec;
else
    valforplot = imagvec/ampimagvec;
end
strforplot = sprintf('%5.2e, ',valforplot);
strimagvec = ['imaginary vector: ',sprintf('%5.2e',ampimagvec),'\cdot [',strforplot(1:end-2),']'];

%at t=0s
t0vec = [real(vxr),real(vyr),real(vzr)];
ampt0vec = norm(t0vec);
valforplot = t0vec/ampt0vec;
strforplot = sprintf('%5.2e, ',valforplot);
strt0vec = ['vibration @ t=0: ',sprintf('%5.2e',ampt0vec),'\cdot  [',strforplot(1:end-2),']'];
%quiver3(x,y,z,u,v,w,'k')

ellipticity = ampimagvec./amprealvec;
ellipticity(ellipticity>1) = ellipticity(ellipticity>1).^(-1);
ellipticity(isnan(ellipticity)) = 0;
if ellipticity==1
    strpolarplane = 'circular';
elseif ellipticity<0.001
    strpolarplane = 'linear';
else
    strpolarplane = 'ellipsoidal';
end

quiver3(x2,y2,z2,realvarx,realvary,realvarz,'b-','DisplayName',strrealvec)
hold on
quiver3(x2,y2,z2,imagvarx,imagvary,imagvarz,'b--','DisplayName', strimagvec)
quiver3(x2,y2,z2,real(vxr),real(vyr),real(vzr),'r','DisplayName', strt0vec)
thetavals = 0:pi/10:2*pi;
xval_ellips = realvarx*cos(thetavals)-imagvarx*sin(thetavals)+x2;
yval_ellips = realvary*cos(thetavals)-imagvary*sin(thetavals)+y2;
zval_ellips = realvarz*cos(thetavals)-imagvarz*sin(thetavals)+z2;
plot3(xval_ellips,yval_ellips,zval_ellips,'k.','displayname',['polarization plane: ',strpolarplane]);
hold off
legend('show')
axis square
vecprod = cross([realvarx,realvary,realvarz],[imagvarx,imagvary,imagvarz]);
if norm(vecprod)~=0
vecprod = vecprod/norm(vecprod);
view(vecprod);
end
title([num2str(nr),':  dipole @ [',num2str([x2,y2,z2]),']'])
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
pause(0.1)
% if strcmpi(strpolarplane,'linear')
% Lims = [get(gca,'Xlim');get(gca,'Ylim');get(gca,'Zlim')];
% set(gca,'Xlim',[min(Lims(:,1)),max(Lims(:,2))]);
% set(gca,'Ylim',[min(Lims(:,1)),max(Lims(:,2))]);
% set(gca,'Zlim',[min(Lims(:,1)),max(Lims(:,2))]);
% 
% end



end