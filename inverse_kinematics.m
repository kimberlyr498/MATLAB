clc

% Input theta r(1)..r(6)
theta  = [0; 0; 0; 0; 0; 0];
%theta = [pi; -pi/2; 0; 0; pi/2; 0];
%theta = [3.94; -1.0472; 1.0472; -.5236; pi/4; 5.2023];

% Denavit-Hartenberg table
pi = 3.1415926535;
%a     = [0; 0; -9.59*25.4; -8.39*25.4; 0; 0];
%alpha = [0; pi/2; 0; 0; pi/2; -pi/2]; 
%d     = [5.98*25.4; 4.75*24.5; 0; -3.66*25.4; 3.28*25.4; 3.28*25.4];

theta = [0+theta(1);pi+theta(2);0+theta(3);0+theta(4);(-pi/2)+theta(5);0+theta(6)];
% Test data from UR site (doesnt even work??)
a = [0; 0; 243.65; 213.25; 0; 82.4];
alpha = [0; pi/2; 0; 0; -pi/2; 0]; 
d = [151.85; 119.88; -92.964; 83.3; 83.4; 0];

% Define and calculate each Transformation Matrix
T = cell(6,1);

for i = 1:numel(theta)
    T{i} = [cos(theta(i)),               -sin(theta(i)),              0,              a(i);
            sin(theta(i))*cos(alpha(i)), cos(theta(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i);
            sin(theta(i))*sin(alpha(i)), cos(theta(i))*sin(alpha(i)), cos(alpha(i)),  cos(alpha(i))*d(i);
            0,                           0,                           0,              1];
end

% Multiply all Trasformation Matrices in cell "T" 
% Obtain new transformation matrix "DH"
DH = T{1};
for n = 2:6
  DH = DH * T{n};
end

% Rotation angles from fixed axis equations
% atan2 returns arctangent of (a,b)=(a/b)
rY = atan2(-DH(3,1),sqrt((DH(1,1)^2.+DH(2,1)^2.)));
rZ = atan2((DH(2,1)/cos(rY)),(DH(1,1)/cos(rY)));
rX = atan2((DH(3,2)/cos(rY)),(DH(3,3)/cos(rY)));

% Assign headers + create table
VarN = {'X_mm', 'Y_mm', 'Z_mm', 'rX_rad', 'rY_rad', 'rZ_rad'};
Table = table(DH(1,4), DH(2,4), DH(3,4), rX, rY, rZ, 'VariableNames' , VarN)
