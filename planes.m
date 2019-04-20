function [ ] = planes( Thetax,Thetay,Thetaz,T )
[y,x]=ndgrid( -0.0716:0.0004:0.078,-0.0716:0.0004:0.078);
z = 0.19*ones(375,375);
plot3(x,y,z,'b.')
x = reshape(x,[1,375^2]);
y = reshape(y,[1,375^2]);
z = reshape(z,[1,375^2]);
R = RotationMat(Thetax,Thetay,Thetaz);
T = T*ones(1,375^2);
A = R*[x;y;z]+T;
xp = A(1,:);
yp = A(2,:);
zp = A(3,:);
xp = reshape(xp,[375,375]);
yp = reshape(yp,[375,375]);
zp = reshape(zp,[375,375]);
hold all
plot3(xp,yp,zp,'r.')

end

