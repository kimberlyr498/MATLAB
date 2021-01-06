% Kimberly Rodriguez
% EGR 5110 Numerical Methods 
% March 4, 2020

                        % Homework Assignment 2 - Part 2
                 % Solving Sets of Non-Linear Algbraic Equations
         % Use Modified Secant Method to Calculate Flowrate in Pipe Network
           
function [Flowrates] = RootFinder_P2(dens, visc, rough, eps)

%% Problem Specifics and Parameters
d = [.1234 .1596 .1234 .1234 .1762 .0968 .1234 .1098];        % Diameter of each link
L = [200 100 260 200 100 200 300 450];                        % Length of each pipe respectively
N = 8;                                                        % Number of links

%% Plug in Delta to Jacobian
Q = [1.1 0.9 0.2 0.9 0.1 0.4 0.5 0.3];        % Initial Guesses for each Pipe
f = zeros(1,N);                               % Initiate friction values' matrix
J = zeros(1,N);                               % Initiate Jacobian matrix functions
dt = 0.01;                                    % Set theta for derivatives definition 

iter = 1;                                     % Counter for loop iterations

for t = 1: 1000000
   
    for n = 1:N
        Re = 4*abs(Q(n))*dens/(pi*d(n)*visc);                                                           % Reynolds Number
    	f(n) =((-1.8)*log10((((rough/(d(n)*3.7))^1.11) + (6.9/Re))))^(-2.);                             % Friction of each branch/pipe
        J(n) =8*f(n)*L(n)*((Q(n)+dt)*abs(Q(n)+dt)-(Q(n)*abs(Q(n)))) / ((pi^2)*9.81*(d(n)^5)*dt);        % Jacobian matrix equations
    end 

B =[(-2000/3600)+Q(1)+Q(2);                   % Solution matrix, where calcualted flowrates are input 
    (200/3600)+Q(4)-Q(3)-Q(2);                % Convegence when all elements in B are close or equal to zero
    (-1000/3600)+Q(7)+Q(3)+Q(8);
    (1300/3600)-Q(7)-Q(4)+Q(5);
    (800/3600)-Q(5)-Q(1)+Q(6);
    (8/(9.81*(pi^2))) * (f(7)*L(7)*Q(7)*abs(Q(7))/(d(7)^5) + f(5)*L(5)*Q(5)*abs(Q(5))/(d(5)^5) - f(8)*L(8)*Q(8)*abs(Q(8))/(d(8)^5) + f(6)*L(6)*Q(6)*abs(Q(6))/(d(6)^5));
    (8/(9.81*(pi^2))) * ((-f(7))*L(7)*Q(7)*abs(Q(7))/(d(7)^5) + f(3)*L(3)*Q(3)*abs(Q(3))/(d(3)^5) + f(4)*L(4)*Q(4)*abs(Q(4))/(d(4)^5));
    (8/(9.81*(pi^2))) * (f(1)*L(1)*Q(1)*abs(Q(1))/(d(1)^5) - f(5)*L(5)*Q(5)*abs(Q(5))/(d(5)^5) - f(2)*L(2)*Q(2)*abs(Q(2))/(d(2)^5) - f(4)*L(4)*Q(4)*abs(Q(4))/(d(4)^5))];
    
    if sqrt(B(1)^2 + B(2)^2 + B(3)^2 + B(4)^2 + B(5)^2 + B(6)^2 + B(7)^2 + B(8)^2) < eps
       Flowrates = [Q(1) Q(2) Q(3) Q(4) Q(5) Q(6) Q(7) Q(8)];        % Output flowrates if termination criteria is met
       break                                                         % Stop the code if convergence occurs at roots of equations
    end
    
A = [1 1 0 0 0 0 0 0;                        % Jacobian Matrix
    0 -1 -1 1 0 0 0 0;                       % dF/dy = [d(xi-delta)- d(xi)] / delta
    0 0 1 0 0 0 1 1;
    0 0 0 -1 1 0 -1 0;
    -1 0 0 0 -1 1 0 0;           
    0 0 0 0 J(5) J(6) J(7) -J(8);            % Input values of derivatives with flowrate input values
    0 0 J(3) J(4) 0 0 -J(7) 0;
    J(1) -J(2) 0 -J(4) -J(5) 0 0 0];        

negB = B*(-1);                               % -F = Jacobian*(Xdelta), in our case -B = A*Q
delta = GaussElim(A,negB);                   % Use Gaussian Elimination to solve Xdelta
Q = Q + delta;                               % Add delta to previous flowrate, Q, and repeat loop
iter = iter + 1;
end

function [X] = GaussElim(A,B)
%% Augmented matrix [A B]
ab = [A B]; 
[n,m]= size(ab);                     % [rows, colums] of [A B]

%% Gaussian Elimination 
for i = 1:n
        
    if i == n                        % Terminate loop on last row
        ab(i,:) = ab(i,:)/ab(i,i);   % Divide entire row by first value
    end
        
    pivot = ab(i,i);              % Use first value on column as pivot
    abmax = abs(ab(i:n,i));       % Convert column to positive integers
    max = abmax(1);
                            % Compare rows for highest value by magnitude
    for c = 2:(n-i+1)
        if abmax(c) > max         
             pivot = ab(i+c-1,i);  % Set as pivot if larger value exists
             max = abmax(c);       % Assign as new value for indexing
        end 
    end
                                  % Swap pivot row and current row
    for w = i:n
         if ab(w,i) == pivot        
              old = ab(i,:);      % Divide entire row by 1st real value:
              ab(i,:) = ab(w,:);  % to equate diagonals to 1
              ab(w,:) = old;
              ab(i,:) = ab(i,:)/ab(i,i);
         end
    end
                                  % Subtract rows by modified first row
    for z = i+1:n                 % R_rem-R_cur*(R_rem(i,i))
         ab(z,:) = ab(z,:) - ab(i,:)*(ab(z,i));
    end 
end

%% Backward Substitution
% Solve for each variable x(1) -> x(m)

% Seperate matrix: variable multipliers (U) and constants (C)
U = ab(:,1:m-1);                   % All except last column
C = ab(:,m);                       % Last column
X(n) = C(n)/U(n,n);                % Define known X(final) variable

for j = n-1:-1:1                   % Decrease index for variable analysis 
    for g = j+1
        sig = (-1)*U(j,g).*X(g);   % Sigma(u(i,j)*x(j,j))
        if j <= n-2                 
        for w = j+2:n
            sig = sig - U(j,w).*X(w);
        end 
        end 
    end
    X(j) = (C(j) + sig)./U(j,j);   % Backward substituion 
end

end                                % End GaussElim.m  

end                                % End RootFinder_P2.m

%% Graveyard 
% f1 = (1/(-1.8*log10(((rough/d(1))/3.7)^1.11 + 6.9/(4*abs(Q(1))*(dens/(pi*d(1)*visc))))))^2;
% f2 = (1/(-1.8*log10(((rough/d(2))/3.7)^1.11 + 6.9/(4*abs(Q(2))*(dens/(pi*d(2)*visc))))))^2;
% f3 = (1/(-1.8*log10(((rough/d(3))/3.7)^1.11 + 6.9/(4*abs(Q(3))*(dens/(pi*d(3)*visc))))))^2;
% f4 = (1/(-1.8*log10(((rough/d(4))/3.7)^1.11 + 6.9/(4*abs(Q(4))*(dens/(pi*d(4)*visc))))))^2;
% f5 = (1/(-1.8*log10(((rough/d(5))/3.7)^1.11 + 6.9/(4*abs(Q(5))*(dens/(pi*d(5)*visc))))))^2;
% f6 = (1/(-1.8*log10(((rough/d(6))/3.7)^1.11 + 6.9/(4*abs(Q(6))*(dens/(pi*d(6)*visc))))))^2;
% f7 = (1/(-1.8*log10(((rough/d(7))/3.7)^1.11 + 6.9/(4*abs(Q(7))*(dens/(pi*d(7)*visc))))))^2;
% f8 = (1/(-1.8*log10(((rough/d(8))/3.7)^1.11 + 6.9/(4*abs(Q(8))*(dens/(pi*d(8)*visc))))))^2;

% a = (8*f(5)*L(5)/((pi^2)*9.81)*(d(5)^5))*((Q(5)+dt)*abs(Q(5)+dt)-(Q(5))*abs(Q(5)));
% b = -(8*f(6)*L(6)/((pi^2)*9.81)*(d(6)^5))*((Q(6)+dt)*abs(Q(6)+dt)-(Q(6))*abs(Q(6)));
% c = (8*f(7)*L(7)/((pi^2)*9.81)*(d(7)^5))*((Q(7)+dt)*abs(Q(7)+dt)-(Q(7))*abs(Q(7)));
% d = -(8*f(8)*L(8)/((pi^2)*9.81)*(d(8)^5))*((Q(8)+dt)*abs(Q(8)+dt)-(Q(8))*abs(Q(8)));
% e = (8*f(3)*L(3)/((pi^2)*9.81)*(d(3)^5))*((Q(3)+dt)*abs(Q(3)+dt)-(Q(3))*abs(Q(3)));
% f = -(8*f(4)*L(4)/((pi^2)*9.81)*(d(4)^5))*((Q(4)+dt)*abs(Q(4)+dt)-(Q(4))*abs(Q(4)));
% w = (8*f(1)*L(1)/((pi^2)*9.81)*(d(1)^5))*((Q(1)+dt)*abs(Q(1)+dt)-(Q(1))*abs(Q(1)));
% v = -(8*f(2)*L(2)/((pi^2)*9.81)*(d(2)^5))*((Q(2)+dt)*abs(Q(2)+dt)-(Q(2))*abs(Q(2)));