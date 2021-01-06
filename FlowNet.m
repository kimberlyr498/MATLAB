clear all 
c=.5;
syms x y
phi = 3/(2*pi)*log((x^2)+(y^2))+6/pi*atan(y/x);
p = fcontour(phi, [-4,4],'r');
hold on
psi = 3/(pi)*atan(y/x)-3/pi*log((x^2)+(y^2));
h = fcontour(psi,[-4,4], 'k');
[x,y] = meshgrid(-4:4,-4:4);
u= 3*(x-2*y)./((x.^2+y.^2)*pi);
v= 3*(y+2*x)./((x.^2+y.^2)*pi);
quiver(x,y,u,v,'k');
hold off
title('Flow Net')
legend([p h],'Phi','Psi')
xlabel('x')
ylabel('y')