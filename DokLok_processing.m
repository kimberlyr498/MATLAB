clear all;

%!          Enter username,and copy and paste file name into trial 
%                -File must begin with Oil name 1st letter; J, E, X, or A
%                -File name must include yyyymmdd in the last 8 characters

    username= ('Kimberly');
    trial=    ('AW46_20180815_30C_20180815'); 

%!         Enter temperature in Celsius 
    temperatureC=  30;
    
%!         Verify path i.e. change "desktop" to "dowloads" etc., followed by folder name         
    read_path=      ['C:\Users\'  username '\Desktop\PetroCan Tests Streamdata\' trial '.csv'];
    write_path=     ['C:\Users\' username '\Desktop\PetroCan Analysis\' trial '.csv'];
    [num,text,raw]= xlsread(read_path);
    
%% Processing

    %Temp Celsius to Kelvin
    tempK=[temperatureC+273];
   
    %Number of samples
    n=numel(raw(:,1));
    sample=[0:n,1];
    
    %define columns
    time=raw(:,1);
    torque=raw(:,2);
    speed=raw(:,3);
    pressure1=raw(:,4);
    pressure2=raw(:,5);
    temp=raw(:,6);
    locked=raw(:,7);
    unlocked=raw(:,8);
    
%DATE- obtained from file name
    day=[];    
    for i = numel(trial)
        datematrix=trial;
        date=datematrix(i-7:i);
        day= [date(5:6) '/' date(7:8) '/' date(1:4)];
    end

%Hydraulic Fluid name and Walther constants
    if trial(1) == 'J'
            fluid= 'Enviroguard';
            a=0;
            b=15;
    elseif trial(1)=='E'
            fluid= 'Extreme';
            a=1.71465;
            b=4.46553;
    elseif trial(1)=='X'
            fluid= 'XV';
            a=2.98699;
            b=7.67571;
    elseif trial(1)=='A'
            fluid= 'AW46';
            a=3.66598;
            b=9.37336;
    end
    
%L+UL
lul(:)=0;
    for f=2:n 
    lock= raw(2:n,7);   lockedd=cell2mat(lock);
    unlock= raw(2:n,8); unlockedd=cell2mat(unlock);
    lul(f)= lockedd(f-1)+unlockedd(f-1); 
    end 
    
%RUN
runn(:)=1;
    for z=(2:n)
        speedindex(1)=1;
        speedindex(2:n)=cell2mat(raw(2:n,3));
     if speedindex(z-1)> 3 && speedindex(z)<3
         runn(z)=runn(z-1)+1;
     else
         runn(z)=runn(z-1);
     end
    end
    
%position
   position{:}='unlocked';
   t= rem(runn,2);
    for l=2:n
     if t(l)==1 
        position{l}=('unlocked');
     elseif t(l)==0 
        position{l}=('locked');
     end
    end

%DELAY
delay(1:n)=0;    
    for l=2:n
        speedind=cell2mat(raw(2:n,3));
     if t(l)==1 && speedind(l-1)> 3
        position{l}=('up');
        delay(l)=0.1;
     elseif t(l)==0 && speedind(l-1)> 3
        position{l}=('down');
        delay(l)=0.1;
     end
    end
    
%KinVis
%kinvis(:)=0;
%kinv(:)=0;
%for h=2:n
%    tempC=cell2mat(raw(2:n,6));
%    kinv(h)=10.^(a+(b*(log10(tempC(h-1)+273.15))));
%    kinvis(h)=(10.^(kinv(h)))-0.7;
%end

%Power
power(:)=0;
for p=2:n
    torquev=cell2mat(raw(2:n,2));
    rpm=cell2mat(raw(2:n,3));
    power(p)=((torquev(p-1)*rpm(p-1))/63025);
end 
    
% Place names & headers in new sheet
    preservecolumns=        [1 2 3 4 5 6 7 8];
    newcolumns=             [6 7 8 9 10 11 13 14];
    headers(1,newcolumns) = raw(1,preservecolumns);
    headers{1} = 'Date';       headers{2} = 'Fluid';
    headers{3} = 'Run';        headers{4}= 'Temp [C]';
    headers{5}='Sample';       headers{12}='Power [hp]';  
    headers{15} = 'L+UL';      headers{16} = 'Delay [s]'; 
    headers{17}='Position';    %headers{18}='KinVis';
    
%Place new columns in new sheet
    cell(1,:)=headers;
        for j=2:n
            cell{j,1}=day; 
            cell{j,2}=fluid;
            cell{j,3}=runn(j);
            cell{j,4}=temperatureC;
            cell{j,5}=sample(j);
            cell{j,6}=time{j};
            cell{j,7}=torque{j};
            cell{j,8}=speed{j}; 
            cell{j,9}=pressure1{j};
            cell{j,10}=pressure2{j};
            cell{j,11}=temp{j}; 
            cell{j,12}=power(j);
            cell{j,13}=locked{j};
            cell{j,14}=unlocked{j};
            cell{j,15}=lul(j);
            cell{j,16}=delay(j);
            cell{j,17}=position{j};
            %cell{j,18}=kinvis(j);
        end

%% Place new file in folder (cannot locate filewrite, change file path on line 16) 

filewrite=[write_path];
sheet = 1;
xlRange = 'A1';
xlswrite(filewrite,cell,sheet,xlRange)
disp(['File ''' trial '.xls'' Sucessfully Processed'])