clc
clear 

% Author: Kimberly Rodriguez
% ME5711 Fundamentals of Robotics
% Filtering Includes: Histogram, Median, Mean, Sobel-X/Y, 
%  and Laplacian Filters respectively
% November 17, 2019
% Lab 4

%% Import Image
% Change .png/file directory for each case
img = imread('3.png');

%% Histogram

% Derive values form histogram to modify plot for clarity
[pixelCounts, grayLevels] = imhist(img,256);
indexvalues = pixelCounts;
avg = .005*max(pixelCounts) - mean(pixelCounts);

[values] = find(pixelCounts > avg);

for n = 1:numel(pixelCounts)
    if pixelCounts(n) > avg
        pixelCounts(n) = avg;
    end
end

figure, 
bar(grayLevels, pixelCounts);
% Add note regarding reduced bins that are shortened
x = mean(values);
y = max(pixelCounts)-2;
txt = transpose(values);
txt = regexprep(num2str(txt),'\s+',', ');
txt = {'pixel values of' txt 'exceed avg of' avg};
text(x,y,txt);

%% Binary Image
%Modify Rows and Columns
    %Create 256+(2*(p-1)) matrix with zeros
    %This allows calculation of edges without predefining
n = 3; %3x3 matrix (8-bit) being calculated
%empty = uint8(zeros(size(img)+(n-1)));
empty = zeros(size(img)+(n-1));

bnr = empty;
%Fill new array with previous values from "img"
for x = 1:size(img,1)
    for y = 1:size(img,2)
        bnr(x+n-2,y+n-2) = img(x,y);
    end
end 

for r = 1:length(bnr)
    for c = 1:length(bnr)
    if bnr(r,c) < avg
        bnr(r,c) = 0;
    elseif bnr(r,c) > avg
        bnr(r,c) = 256;
    end
    end 
end
subplot(2,3,1);
image(bnr);
title('Binary Image');

%% Median Filter
mdn = empty;
%Fill new array with values from "img"
for x=1:size(img,1)
    for y=1:size(img,2)
        mdn(x+n-2,y+n-2)=img(x,y);
    end
end 

%Calculate only the values starting at row and comlumn (2,2)
%And ending at (257,257) to avoid edges missing surrounding values
 for i=2:(size(mdn,1)-1)
    for j=2:(size(mdn,2)-1)
        surround = [mdn(i-1,j-1) mdn(i-1,j) mdn(i-1,j+1) mdn(i,j-1) mdn(i,j) mdn(i,j+1) mdn(i+1,j-1) mdn(i+1,j) mdn(i+1,j+1)];
        arrng = sort(surround); %Ascending order
        mdn(i-1,j-1) = arrng(ceil((n)^2/2)); %Median assigned as new value
    end
 end
subplot(2,3,2);
imshow(mdn);
title('Median Filter');

%% Mean Filter
meanfilt = empty;
%Fill new array with values from "img"
for x=1:size(img,1)
    for y=1:size(img,2)
        meanfilt(x+n-2,y+n-2)=img(x,y);
    end
end 

 for i=2:(size(meanfilt,1)-1)
    for j=2:(size(meanfilt,2)-1)
        surround = [meanfilt(i-1,j-1) meanfilt(i-1,j) meanfilt(i-1,j+1) meanfilt(i,j-1) meanfilt(i,j) meanfilt(i,j+1) meanfilt(i+1,j-1) meanfilt(i+1,j) meanfilt(i+1,j+1)];
        val = mean(surround); %Average of 3x3 matrix
        meanfilt(i-1,j-1) = val; %Mean assigned as new value
    end
 end
subplot(2,3,3);
imshow(meanfilt);
title('Mean Filter');

%% Sobel Operator
%Horizontal operation
sobelx = empty;
%Fill new array with values from "img"
for x=1:size(img,1)
    for y=1:size(img,2)
        sobelx(x+n-2,y+n-2)=img(x,y);
    end
end 

for i=2:(size(sobelx,1)-1)
    for j=2:(size(sobelx,2)-1)
        operator = [sobelx(i-1,j-1) sobelx(i-1,j) sobelx(i-1,j+1) sobelx(i,j-1) sobelx(i,j) sobelx(i,j+1) sobelx(i+1,j-1) sobelx(i+1,j) sobelx(i+1,j+1)];
        hsobel= transpose([-1 -2 -1 0 0 0 1 2 1]); %Sobel matrix horizontal, Wx
        val = operator*hsobel; %Multiply our matrix by sobel's
        sobelx(i-1,j-1) = sum(val); %Sum of product assigned as new value
    end
 end
subplot(2,3,4);
imshow(sobelx);
title('Sobel-X Filter');

%Sobel vertical operation
sobely = empty;
for x=1:size(img,1)
    for y=1:size(img,2)
        sobely(x+n-2,y+n-2)=img(x,y);
    end
end 

 for i=2:(size(sobely,1)-1)
    for j=2:(size(sobely,2)-1)
        operator = [sobely(i-1,j-1) sobely(i-1,j) sobely(i-1,j+1) sobely(i,j-1) sobely(i,j) sobely(i,j+1) sobely(i+1,j-1) sobely(i+1,j) sobely(i+1,j+1)];
        vsobel= transpose([-1 0 1 -2 0 2 -1 0 1]); %Sobel matrix vertical, Wy
        val = operator*vsobel; %Multiply our matrix by sobel's
        sobely(i-1,j-1) = sum(val); %Sum of product assigned as new value
    end
 end
subplot(2,3,5);
imshow(sobely);
title('Sobel-Y Filter');

%% Laplacian Filter
lap = empty;
for x=1:size(img,1)
    for y=1:size(img,2)
        lap(x+n-2,y+n-2)=img(x,y);
    end
end 

%Calculate only the values starting at row and comlumn (2,2)
%And ending at (257,257) to avoid edges missing surrounding values
 for i=2:(size(lap,1)-1)
    for j=2:(size(lap,2)-1)
        operator = [lap(i-1,j-1) lap(i-1,j) lap(i-1,j+1) lap(i,j-1) lap(i,j) lap(i,j+1) lap(i+1,j-1) lap(i+1,j) lap(i+1,j+1)];
        laplacem = transpose([0 -1 0 -1 4 -1 0 -1 0]); %Sobel matrix horizontal, Wx
        val = operator*laplacem; %Multiply our matrix by sobel's
        lap(i-1,j-1) = sum(val); %Sum of product assigned as new value
    end
 end
subplot(2,3,6);
imshow(lap);
title('Laplacian Filter');