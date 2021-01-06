% Kimberly Rodriguez
% EGR 5110 Numerical Methods
% April 23, 2020

                             % Assignment 4
                    % Partial Differential Equations
                % Solving Heat Diffusion in Elliptic PDE

function [T,Ttipsim,Qfinsim] = calcTvstime(T,Nx,Ny,Nt,lam,kcond,h,dx,dt,Lx,Ly,Lz,Bi,Tb,Tinf)

% Insert the base temperature to first column
T(1,:,:) = Tb;
sc = 5;                 % Scenario number to complete table
%% Calculate the temperature distribution at each time step
for t = 1:Nt-1
    
    %Interior Nodes
    for y = 2:Ny-1
        for x = 2:Nx-1
            T(x,y,t+1) = lam*(T(x-1,y,t)+T(x+1,y,t)+T(x,y+1,t)+T(x,y-1,t))+(1-4*lam)*T(x,y,t);
        end 
    end 
    
    % Top and Bottom Nodes
    for x = 2:Nx-1
        T(x,Ny,t+1) = lam*(2*T(x,Ny-1,t)+T(x+1,Ny,t)+T(x-1,Ny,t)+2*Bi*Tinf)+(1-4*lam-2*Bi*lam)*T(x,Ny,t);
        T(x,1,t+1) = lam*(2*T(x,2,t)+T(x+1,1,t)+T(x-1,1,t)+2*Bi*Tinf)+(1-4*lam-2*Bi*lam)*T(x,1,t);
    end
    
    % Right Corner Nodes
    for x = Nx
       T(x,Ny,t+1) = 2*lam*(T(x,Ny-1,t)+T(x-1,Ny,t)+2*Bi*Tinf)+(1-4*lam-4*Bi*lam)*T(x,Ny,t);
       T(x,1,t+1) = 2*lam*(T(x,2,t)+T(x-1,1,t)+2*Bi*Tinf)+(1-4*lam-4*Bi*lam)*T(x,1,t); 
    end
        
    % Outer Edge: Tip of Fin
    for y = 2:Ny-1
       T(Nx,y,t+1) = lam*(2*T(Nx-1,y,t)+T(Nx,y-1,t)+T(Nx,y+1,t)+2*Bi*Tinf)+(1-4*lam-2*Bi*lam)*T(Nx,y,t);
    end
    
end

%% Calculate the average temperature at the tip at the last time step

Ttipsim=0;

for i = 1:Ny
    Ttipsim = Ttipsim + T(Nx,i,Nt);
end 

Ttipsim = Ttipsim/Ny;

%% Calculate the heat rate into the fin at the last time step
Qfinsim = 0;

for i = 2:Ny-1
    Qcond = kcond*dx*Lz*(T(1,i,Nt)-T(2,i,Nt))/dx;
    Qfinsim = Qfinsim + Qcond;
end

% Bottom and Top Nodes
Qcond = kcond*dx*Lz*(T(1,1,Nt)-T(2,1,Nt))/dx;
Qfinsim = Qfinsim + Qcond;

%% Create an animation of the temperature distribution (only plot 100 time steps)
figure(1)
times = dt*Nt/60;
sx = Lx/Nx;
sy = Ly/Ny;
yval = 0.0:sx:(Lx-sx);
xval = 0.0:sy:(Ly-sy);
[y,x] = meshgrid(xval,yval);
z = T(:,:,Nt);
contour(x,y,z,'ShowText','on')
xlabel('length (m)')
ylabel('length (m)')
title(['Temperature Distribution (C), Scenario:' num2str(sc)])
text(.025,.005,['time =' (num2str(times)) 'minutes']);
colorbar

% Animation video
figure(2)
% Create File
myvideo = VideoWriter('calTvstime.mp4','MPEG-4'); 
open(myvideo);

frame = getframe(gcf);            % Obtains the current figure as movie frame
writeVideo(myvideo,frame);        % Begin File Writing

ts = 100;                         % Time steps for video

for t = 1:(Nt/ts):Nt
    yval = 0.0:sx:(Lx-sx);
    xval = 0.0:sy:(Ly-sy);
    [y,x] = meshgrid(xval,yval);
    z = T(:,:,t);
    contour(x,y,z,'ShowText','on')
    xlabel('length (m)')
    ylabel('length (m)')
    title('Temperature Distribution (C)')
    text(.025,.005,['time =' (num2str(dt*t/60)) 'minutes']);
    colorbar
    frame = getframe(gcf);
    writeVideo(myvideo,frame);
end
close(myvideo)

% End of code