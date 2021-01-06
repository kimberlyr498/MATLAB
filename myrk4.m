% Kimberly Rodriguez
% EGR5110 S20 
% March 31, 2020
 
                        % Assignment 3 RK4
                 % Use Runge-Kutta 4th Order Method
       % Plot Trajectory of a Satelite in Orbit Around Barycenter

function [x,y,vx,vy,t] = myrk4(t0,tf,dt,x0,y0,vx0,vy0,fd,eps)

%% Define Known Variables
mu = 1/82.45;        % Earth/Moon Mass -(location of earth on x-axis)
mup = 1-mu;          % Location of Moon on x-axis
N = (tf-t0)/dt;      % Number of steps needed to calcualte motion over time-span

%% Ordinary Differential Equations
% Motion of a body in orbit about two large bodies
 % Equations are rewriten as a first-order differential equation
  % f1=x, f2=x'=vx, f3=y, f4=y'=vy

fprime1 = @(f1, f2, f3, f4) f2;
fprime2 = @(f1, f2, f3, f4) 2*f4 + f1 - (mup*(f1 + mu))/(sqrt(((f1+mu)^2) + (f3^2))^3) - (mu*(f1 - mup))/(sqrt(((f1-mup)^2) +(f3^2))^3) - fd*f2;
fprime3 = @(f1, f2, f3, f4) f4;
fprime4 = @(f1, f2, f3, f4) f3 - 2*f2 - (mup*f3)/(sqrt(((f1+mu)^2) + (f3^2))^3) - (mu*f3)/(sqrt(((f1 - mup)^2) + (f3^2))^3) - fd*f4;

f = {fprime1 fprime2 fprime3 fprime4};        % An array of ODEs listed above for calling

%% Runge-Kutta Method

% Fill arrays with zeros
k1 = zeros(4,1); k2 = zeros(4,1);        % Partial Slopes
k3 = zeros(4,1); k4 = zeros(4,1);
phi = zeros(4,1);                        % phi = dt*(k1+2k2+2k3+k4)/6
r1 = zeros(10,1);                        % Distance from earth; used in error calcualtion
x = zeros(N+1,1);                        % Horizontal distance from Barycenter (center of mass between Earth and the moon)
y = zeros(N+1,1);                        % Vertical Distance from Barycenter
vx = zeros(N+1,1);                       % Velocity in the x-direction
vy = zeros(N+1,1);                       %        "        y-direction
t = zeros(N+1,1);                        % Independent variable (time,s)

for n = 1:5                              % Loop for changes in step-size  
    for i = 1:N                          % Calculate values for N amount of steps {N = (tf-t0)/dt}
        x(1) = x0; y(1) = y0;            % Initial Conditions
        vx(1) = vx0; vy(1) = vy0; 
        t(1)=t0;
                                         % Partial Slopes calculate next
                                         %  value; using Ftemporary (Ft)
        for j = 1:4
            k1(j) = f{j}(x(i),                vx(i),               y(i),               vy(i));
        end 
        for j = 1:4                      % Ft(respsective) = 'variable' + k1(respective)*dt/2
            k2(j) = f{j}((x(i)+k1(1)*dt/2) , (vx(i)+k1(2)*dt/2) , (y(i)+k1(3)*dt/2) , (vy(i)+k1(4)*dt/2));
        end 
        for j = 1:4                      % Ft. = 'variable' + k2.*dt/2
            k3(j) = f{j}((x(i)+k2(1)*dt/2) , (vx(i)+k2(2)*dt/2) , (y(i)+k2(3)*dt/2) , (vy(i)+k2(4)*dt/2));
        end 
        for j = 1:4                      % Ft. = 'variable' + k3.*dt
            k4(j) = f{j}((x(i)+k3(1)*dt)   , (vx(i)+k3(2)*dt)   , (y(i)+k3(3)*dt)   , (vy(i)+k3(4)*dt));
        end 
        for j = 1:4
            phi(j) = dt*(k1(j) + 2*k2(j) + 2*k3(j) + k4(j))/6;
        end
        
        x(i+1) = x(i) + phi(1);          % Predict Following Values  
        vx(i+1) = vx(i) + phi(2);        % using current value and partial slopes
        y(i+1) = y(i) + phi(3);
        vy(i+1) = vy(i) + phi(4);        
        t(i+1) = t(i) + dt;
    end
    
    r1(n) = sqrt((x(end)+mu)^2 + y(end)^2);                  % Distance from earth at end of tracking data     

    if n > 1 && abs((r1(n) - (r1(n-1))) / r1(n)) < eps       % Compares new step size to previous for major shifts
        break                                                % Break loop if current step size lies within error margin
    end
    
    N = N*2;                             % Double the number of steps to improve accuracy of results
    dt = dt/2;                           % Respectively half the step-size
end

%% Plot
figure(1)
% Trajectory
subplot(2,1,1)
plot(x,y)
legend('384,400^{km}/_{unit}')
xlabel('X-direction (unit)')
ylabel('Y-direction (unit)')
text(-mu,0,'Earth')
text(mup,0,'Moon')
title('Trajectory of Spaceship about Barycenter')
grid on 
%Phase Space
subplot(2,1,2)
plot(sqrt((x+mu).^2+y.^2),sqrt(vx.^2+vy.^2))
legend('^{384,400^{km}/_{unit}}/_{3.75e5 ^{s}/_{unit}}')
xlabel('Distance (unit)')
ylabel('Velocity (unit(^{km}/_{s}))')
title('Phase Space')
grid on 

%% Simulation; using 'drawnow' command
% Animate motion of spaceship
figure(2)
h = animatedline('Color','b');
title('Spaceship Trajectory');
axis([-1.5 1.5 -1 1]);
for i = 1:length(t)
    addpoints(h,x(i),y(i))          
    drawnow
end 