



function R = RotationMat(theta_x,theta_y,theta_z)

theta_x = theta_x/180*pi;
Rx = [ 1     0        0 
        0  cos(theta_x) -sin(theta_x)
        0  sin(theta_x) cos(theta_x)];
theta_y = theta_y/180*pi;
Ry = [ cos(theta_y)  0  sin(theta_y)
         0       1    0
       -sin(theta_y)  0  cos(theta_y)];

theta_z = theta_z/180*pi;
Rz = [ cos(theta_z) -sin(theta_z) 0
        sin(theta_z)  cos(theta_z) 0
          0       0      1];

R = Rx*Ry*Rz;

end
