%% PHASE SPACE RECONSTRUCTION: TIME DELAY

% Author: Santiago Duarte (Escuela Colombiana de Ingenieria Julio Garavito), santiagoduare09@gmail.com
% March 2017

%INPUTS:
%=============
% PRECIPITATION  : 
%
% OUTPUT:
%=============
% LAG : 

function [LAG]=Delay(DATA,FOLDER)

if nargin < 2
    FOLDER = '';
end

%AUTOCORRELATION FUNCTION
da=1;dm=1;
[Auto] = autocorr(DATA,length(DATA)-1);
while Auto(da)*Auto(da+1)>0;
    da=da+1;
end
da=da+1;
if da~=0
    LAG(1,1)=da;
else
    LAG(1,1)=1;
end

%MUTUAL INFORMATION FUNCTION

save Precipitation.dat DATA -ascii
system([mfilename('fullpath'),'\..\Auxiliar\cygwin64\mutual.exe -V0 ',pwd,'\Precipitation.dat -b20 -D',num2str(length(DATA)),' -o',pwd,'\Mutua.dat']);
fileID = fopen('Mutua.dat');
Mut = textscan(fileID, '%s %f32');
Mutua=cell2mat(Mut(2));
Mutua=Mutua(2:end);
fclose(fileID);
while Mutua(dm+1)<Mutua(dm);
    dm=dm+1;
end
if dm~=0
    LAG(2,1)=dm; 
else
    LAG(2,1)=1; 
end

%SAVE DELAY FIGURES (OPTIONAL)
if ~strcmp(FOLDER,'')
    %PLOT MUTUAL INFORMATION AND AUTOCORRELATION
    figure('Visible','off');
    hold all
    plot(Mutua,'-b')
    plot(Auto,'-r.');
    title('Mutual Information and Autocorrelation'); xlabel('Time Delay(dias)');ylabel('Mutual Information and Autocorrelation');xlim([0 max(da,dm)]+5);legend('Mutual Information','Autocorrelation'); grid on;
    saveas(gcf,strcat(FOLDER,'Time_delay'),'jpg')  
    hold off
    
    %PLOT RECONSTRUCTION PHASE SPACE IN 3D
    Del=min(LAG); Dimension=3;
    Delay=NaN(length(DATA)-((Dimension-1)*Del),Dimension);
    for i=1:Dimension
        Delay(:,i)=DATA(1+Del*(Dimension-i):end-Del*(i-1));
    end
    plot3(Delay(:,1),Delay(:,2),Delay(:,3),'-b');view([135 45])
    title('Phase Space Reconstruction in 3 Dimensions');xlabel('t (days)');ylabel('t- \tau (days)');zlabel('t- 2\tau (days)');grid on;
    saveas(gcf,strcat(FOLDER,'Phase_space'),'jpg')  
    hold off
    
    %SPACE-TIME SEPARATION PLOT
    system([mfilename('fullpath'),'\..\Auxiliar\cygwin64\stp.exe -V0 ',pwd,'\Precipitation.dat -d',num2str(Del),' -m',num2str(Dimension),' -o',pwd,'\STP.dat']); 
    STP=load('STP.dat');pl=0;   
    for n=1:100
        pi=pl+n; pl=100*n;
        plot(STP(pi:pl,1),STP(pi:pl,2),'-b');hold all
    end
    title('Space-time Separation Plot');xlabel('Relative Time');ylabel('Distance');
    saveas(gcf,strcat(FOLDER,'STP'),'jpg')  
    hold off
    delete('STP.dat','Mutua.dat')   
end
end