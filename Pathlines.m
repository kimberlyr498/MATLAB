%Define arrays x, y, u, and v.
[x,y] = meshgrid(0:0.1:1,0:0.1:1);
u = x;
v = -y;
% Create a quiver plot of the data.
%Plot streamlines that start at different
%points along the line y=1
figure
quiver(x,y,u,v)
startx = 0.1:0.1:1;
starty = ones(size(startx));
streamline(x,y,u,v,startx,starty);
title('Pathline')