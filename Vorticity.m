%Homework 1_Advanced Fluid Dynamics
%No. 19
C = .5; %constant
v = 1; %kinematic viscosity in m^2/s
t = linspace(10,50,5); %take time at 10, 20, 30, 40, and 50 seconds
r = linspace(0,20); %take radius values from 0 to 20 meters
v_theta = zeros(1,100); %initialize the theta component of velocity
vorticity = zeros(1,100); %initialize the vorticity
for i=1:5 %outer for loop to plot at the five sample times between 0 and 50 seconds
for j=1:100 %inner for loop to plot the theta component of velocity and the vorticity versus r
v_theta(j) = (C/r(j))*(1-exp(-r(j)^2/(4*v*t(i))));
vorticity(j) = 2*C*exp(-r(j)^2/(4*v*t(i)))/(4*v*t(i));
end
%plot the theta component of velocity
figure(1)
plot(r,v_theta)
hold on
%plot the vorticity 
figure(2)
plot(r,vorticity)
hold on
end
figure(1)
title('Velocity Profile')
xlabel('Radius (m)')
ylabel('Angular Component of Velocity (m/s)')
legend('t=10', 't=20', 't=30', 't=40', 't=50')
title(legend,'Time (s)')
figure(2)
title('Vorticity Profile')
xlabel('Radius (m)')
ylabel('Vorticity (1/s)')
legend('t=10', 't=20', 't=30', 't=40', 't=50')
title(legend, 'Time (s)')