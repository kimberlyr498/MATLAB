% Author: Kimberly Rodriguez
% ME5711 Fundamentals of Robotics
% 8-chain connectivity
% Outputs the chain directory of the shape
% Up equals 2 moving counter-clockwise for upper right to equate 1
% November 17, 2019
% Lab 4

clc
clear all

mg = imread('3.png');
imshow(mg)
    t = 1;
    chain(t) = 0;
new = zeros(size(mg)+2);
for i = 1:size(mg)
    for n = 1:size(mg)
        new(i+1,n+1) = mg(i,n);
    end
end

mg = new;

for x = 2:(size(mg,1)-1)
    for j = 2:(size(mg,2)-1)
        A = [mg(x-1,j-1) mg(x-1,j) mg(x-1,j+1) mg(x,j-1) mg(x,j+1) mg(x+1,j-1) mg(x+1,j) mg(x+1,j+1)];
        chck = mean(A);
        if chck < 1 
        end
        if mg(x-1,j-1) > 0
        chain(t) = 3; 
        x=x-1;
        j=j-1; 
        t = t+1;  
        end 
       if mg(x-1,j) > 0
        chain(t) = 2; 
        x = x-1;
        t = t+1;
       end 
       if mg(x-1,j+1) > 0
        chain(t) = 1;
        x = x-1;
        j = j+1;
        t = t+1;
        end 
        if mg(x,j-1) > 0
        chain(t) = 4;
        j = j-1;
        t = t+1;
        end 
            if mg(x,j+1) > 0
        chain(t) = 0;
        j=j+1;
        t = t+1;
            end
        if mg(x+1,j-1) > 0
        chain(t) = 5;
        x = x+1;
        j= j-1;
        t = t+1;
        end 
        if mg(x+1,j) > 0
        chain(t) = 6;
        x = x+1;
        t = t+1;
        end 
        if mg(x+1,j+1) > 0
        chain(t) = 7;
        x = x+1;
        j = j+1;
        t = t+1;
        end 
    end 
end 