% Kimberly Rodriguez
% EGR 5110 Numerical Methods
% May 2, 2020

                         % Assignment 6
                % 2D Unconstrained Optimization
       % Use Golden Ration and Gradient for steepest ascend
       
function [xypos1,numsteps1] = calcMaxStudent(f, xi, yi, eps)

% Initiate Position Array
x(1) = xi; y(1) = yi;        % Initial x and y coordinates
numsteps1 = 0;               % Initial number of steps
delta = 0.0001;              % Arbitrary value to calculate derivatives
        
for i = 1:5e5
    numsteps1 = numsteps1 + 1;        % Step counter
    
    % Calculate derivative df/dx and df/dy
    % Uses center-difference method
    derivx(i) = (f(x(i)+delta,y(i)) - f(x(i)-delta,y(i)))/(2*delta);
    derivy(i) = (f(x(i),delta+y(i)) - f(x(i),y(i)-delta))/(2*delta);
    
    % Reference Golden Ratio Method in function (end of code)
    % function solves for optimal (maximum) height(h) to travel
    h = goldenratio(x(i),y(i),derivx(i),derivy(i),f,eps);
    
    % Terminate loop when flat area is reached; f'(x,y)=0
    if sqrt(derivx(i)^2+derivy(i)^2) < eps
        x = x.';                      % Transpose to column vector
        y = y.';
        xypos1 = [x,y];               % Combine x y for output
        break
    end
    
    % Approximate position coordinates using exact line search method
    x(i+1) = x(i) + derivx(i)*h;
    y(i+1) = y(i) + derivy(i)*h;
end

%% Plot path along countour graph
xval = -10:10;                        % Range assigned for plotting
yval = -10:10;
[n,p] = meshgrid(xval,yval);
z = f(n,p);
contour(n,p,z)                        % Contour plot of function (f)
title('Path Traveled')
colorbar 
xlabel('x-direction')
ylabel('y-direction')
hold on
plot(x,y,'b-o')                             % Hold until path is plotted
hold off                       % includes trajectory on same figure

%% Golden-ratio solver for g'(h)=0
% Optimize g(h)=f(x,y) to find local maximum (h value for max g(h))
function h = goldenratio(x,y,derivx,derivy,f,eps)
g = @(h) f(x+derivx*h,y+derivy*h);              % 1D function from 2D
phi = (-1+sqrt(5))/2;               % L2/L1 = L1/L0 i.e. golden ratio
xa = 0; xb = 100;                 % Arbitrary values for search space
d = phi*(xb-xa);           % Variable distance to generate new values      

x1 = xa+d;                       % New values generated
x2 = xb-d;

for n = 1:5e5
    f1 = g(x1);                  % Calculate g(h) at each point
    f2 = g(x2);
    
    % Recreate search space based on the following criteria
    if f1 > f2                      % If f1>f2; xa=x2, x2=x1
        xa = x2;
        x2 = x1;
        d = phi*(xb-xa);            % find new d using golden ration
        x1 = xa+d;                    % for calculating missing values
    elseif f1 < f2                  % If f2>f1; xb=x1, x1=x2
        xb = x1;
        x1 = x2;
        d = phi*(xb-xa);
        x2 = xb-d;
    end
    
    % Terminate loop if it meets the criteria below
    if (abs(xb-xa)/(abs(x1)+abs(x2))) < eps
        if f1 > f2
            h = x1;
        else
            h = x2;
        end
        break
    end
end
end
end 