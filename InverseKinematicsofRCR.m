clc
clear all
t1 = 1*pi/180;
t2 = 1*pi/180;
t3 = 1*pi/180;
L1 = 6;
L2 = 6;
L3 = 3;
d1 = 1; 
c = sqrt(L2^2+L3^2);

x = cos(t1)*L1 + cos(t1)*cos(t2)*d1 - d1*sin(t1)*sin(t2) + cos(t3)*L2*(cos(t1)*cos(t2) - sin(t1)*sin(t2));
y = L1*sin(t1) + cos(t1)*d1*sin(t2) + cos(t2)*d1*sin(t1) + cos(t3)*L2*(cos(t1)*sin(t2) + cos(t2)*sin(t1));
z = L2*sin(t3);

q3 = asin(z/L2);
q2 = acos(((x^2)+(y^2)-(L1^2)-(L2*cos(t3)+d1)^2)/(2*L1*(L2*cos(t3)+d1)));
q1 = atan(y/x)-atan(((L2*cos(t3)+d1)*sin(t2))/(L1+(L2*cos(t3)+d1)*cos(t2)));