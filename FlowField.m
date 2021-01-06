clear all 
[x,y]=meshgrid(0:10);
u=2*y;
v=x;
figure
quiver(x,y,u,v)
title('Flow Field')