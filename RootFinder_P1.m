% Kimberly Rodriguez
% EGR 5110 Numerical Methods 
% March 4, 2020

                        % Homework Assignment 2 - Part 1
                 % Solving Sets of Non-Linear Algbraic Equations
              % Using Bisecant, Newton-Raphson and Secant Methods

function [fBis, iterBis, fNR, iterNR, fSec, iterSec] = RootFinder_P1(y, dy, fa, fb, f0, f1, f2, eps)

%% Bisection Method
iter = 1;                            % Counter keeps track of number of iterations that the following loop progresses
llim = 1e-10; ulim = 1e-4;           % Lower and Upper limits for epsilon values in plot
while iter < 100                     % Run loop while iterations are still within allowed parameters
    fmid = (fa + fb) / 2;            % Generate new vaule between fmin and fmax 
                                                        
    if y(fa)*y(fmid) > 0             % These values are reassigned if they cross over the axis
        fa = fmid;                   %... to allow convergence on the root of the function
    elseif y(fa)*y(fmid) < 0 
        fb = fmid;
    end

   
    if iter > 1                                 % After first iteration
       tolerance = abs(((fmid-test)/fmid));
       if llim <= tolerance && tolerance <= ulim 
          epsval(iter) = tolerance;
          iterplot(iter) = iter;
       end 
       if tolerance < eps        % Check termination criterion to break loop if necessary
          fBis = fmid;                        % Assign current f value as the solution
          iterBis = iter;                     % Export number of iterations
          break                               % 'break' loop if termination criterion is met
       end
    end
    

test = fmid;                          % Save current "f" to check future termination criteria 
iter = iter + 1;                      % Continue iteration loop if termination criteria is not met
end
 
%% Newton-Raphson Method
iter = 1;                            % Counter keeps track of number of iterations that the following loop progresses

while iter < 100                     % Run loop while iterations are still within allowed parameters
    fnew = f0 - (y(f0)/dy(f0));      % Generate next value using taylor expansion 
       
        if abs(((fnew-f0)/fnew)) < eps        % Check termination criterion to break loop if necessary
            fNR = fnew;                       % Assign current f value as the solution
            iterNR = iter;                    % Export number of iterations
            break                             % 'break' loop if termination criterion is met
        end
        
f0 = fnew;                            % Save current "f" to check future termination criteria
iter = iter + 1;                      % Continue iteration loop if termination criteria is not met
end

%% Secant Method
iter = 1;                            % Counter keeps track of number of iterations that the following loop progresses

while iter < 100                                        % Run loop while iterations are still within allowed parameters
    f3 = f2 - ( y(f2)*(f2 - f1)) / (y(f2)-y(f1) );      % Generate next value using Secant Method 
    tolerance = abs((f3-f2)/f3);                        % Define Tolerance

        if tolerance < eps                % Check termination criterion to break loop if necessary
            fSec = f3;                    % Assign current f value as the solution
            iterSec = iter;               % Export number of iterations
            break                         % 'break' loop if termination criterion is met
        end
        
f1 = f2;
f2 = f3;                              % Save current "f" to check future termination criteria       
iter = iter + 1;                      % Continue iteration loop if termination criteria is not met
end 

%% Plot Secant Method Convergence
figure;
semilogx(epsval,iterplot,'-b','MarkerSize',15)                   
title('Iteration vs. Tolerance')
xlabel('\epsilon')
ylabel('Iteration')