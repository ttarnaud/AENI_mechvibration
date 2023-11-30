% 
theta_test = [pi/4,pi/3]';
phi_test = [pi/3.5,pi/4.1]';
points = [cos(phi_test).*sin(theta_test),sin(phi_test).*sin(theta_test),cos(theta_test)];
[xsphere,ysphere,zsphere] = sphere(100);
angle1 = acos(dot(points(1,:),points(2,:)));
angle2 = 2*pi-angle1
q=@(t) (points(1,:).*sin((1-t').*angle) + points(2,:).*sin(t'*angle))./sin(angle);
t=[0:0.001:1];
qvals = q(t);
surf(xsphere,ysphere,zsphere)
hold on
plot3(qvals(:,1),qvals(:,2),qvals(:,3));
scatter3(points(:,1),points(:,2),points(:,3))
hold off
%%

theta_test = [pi/4,pi/3.5,pi/3,pi/2.7]';
phi_test = [pi/3.5,pi/3.7,pi/3.8,pi/4.1]';
points = [cos(phi_test).*sin(theta_test),sin(phi_test).*sin(theta_test),cos(theta_test)];
v1 = cross(points(1,:),points(2,:));
v2 = cross(points(3,:),points(4,:));
v1n = v1/norm(v1); v2n = v2/norm(v2);
if norm(v1n-v2n)<1e-8
    disp('v1, v2 are identical')
end
D = cross(v1n,v2n);
D1n = D/norm(D);
D2n = -D1n;
%%

angle1 = acos(dot(points(1,:),points(2,:)));
angle2 = 2*pi-angle1
q=@(t,angle) (points(1,:).*sin((1-t').*angle) + points(2,:).*sin(t'*angle))./sin(angle);
t=[0:0.001:1];
qvals1 = [q(t,angle1);q(t,angle2)];

angle3 = acos(dot(points(3,:),points(4,:)));
angle4 = 2*pi-angle3
q=@(t,angle) (points(3,:).*sin((1-t').*angle) + points(4,:).*sin(t'*angle))./sin(angle);
t=[0:0.001:1];
qvals2 = [q(t,angle3);q(t,angle4)];



surf(xsphere,ysphere,zsphere,'FaceColor','y','EdgeColor','none','facealpha',0.5)
%s.EdgeColor = 'none';
hold on
scatter3(points(:,1),points(:,2),points(:,3));
hold on
scatter3(D1n(1),D1n(2),D1n(3),'rX')
scatter3(D2n(1),D2n(2),D2n(3),'kX')
plot3(qvals1(:,1),qvals1(:,2),qvals1(:,3));
plot3(qvals2(:,1),qvals2(:,2),qvals2(:,3));
hold off
axis square