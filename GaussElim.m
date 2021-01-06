% By: Kimberly Rodriguez
% ENGR5110 Numerical Methods
% Updated: Feb. 12, 2020
              % Gaussian Elimination for Solving Systems of Lin. Equations
                   % Function calls matrices A(nxn) and B(nx1)
                        % Outputs solution vector X(nx1)

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