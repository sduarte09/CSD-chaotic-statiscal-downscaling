%% CHAOS_PROPERTIES_COMPARATION

% Author: Santiago Duarte (Escuela Colombiana de Ingenieria Julio Garavito), santiagoduare09@gmail.com
% November 2016

%INPUTS:
%=============
% STATION                    :
% CHAOS_DATA_STATION         :
% CLIMATE_MODEL              :
% CHAOS_DATA_CLIMATE_MODEL   :
% INTERVAL_ACCUMULATION_POS  : 
% OUTPUT_FILE_NAME           :
%

function []=Chaos_Comparation(STATION,GCM,CHAOS_STATION,CHAOS_GCM,IAP,FOLDER)

INI=find(GCM(:,1)==STATION(1,1));
FIN=find(GCM(:,1)==STATION(end,1));

%ACUMMULATED STATION
k=1;Info=NaN(length(STATION(:,2)),1);p=CHAOS_STATION{1,IAP};
while k<=(length(STATION(1:end,2))-p)
    if isnan(STATION(k:k+p-1,2))
        Info(k)=NaN;
    else
        Info(k)=nansum(STATION(k:k+p-1,2));% Transforming data to NaN and Numbers
    end
    k=k+p;
end
STATION_CUM=Info(1:p:end); % Accumlated data extraction
if isnan(STATION_CUM(end))
    STATION_CUM(end)=[];
end

%ACUMMULATED GCM
k=1;Info=NaN(length(GCM(INI:FIN,2)),1);p=CHAOS_GCM{1,IAP};
while k<=(length(GCM(INI:FIN,2))-p)
    if isnan(GCM(INI+k:INI+k+p-1,2))
        Info(k)=NaN;
    else
        Info(k)=nansum(GCM(INI+k:INI+k+p-1,2));% Transforming data to NaN and Numbers
    end
    k=k+p;
end
GCM_CUM=Info(1:p:end); % Accumlated data extraction
if isnan(GCM_CUM(end))
    GCM_CUM(end)=[];
end

%ACUMMULATED TIME
TIME_CUM=STATION(1:CHAOS_STATION{1,IAP}:end,1); % Accumlated data extraction
if length(TIME_CUM)>length(STATION_CUM)
    TIME_CUM(end)=[];
end

%% POWER SPECTRAL DENSITY
figure('Visible','off');
subplot(2,2,1); 
[pxx,f] = periodogram(TIME_CUM,STATION_CUM,'onesided',512,1000);
[pxx1,f1] = periodogram(TIME_CUM,GCM_CUM,'onesided',512,1000);
plot(f,10*log10(pxx),'b');hold on
plot(f1,10*log10(pxx1),'r-.')
xlabel('Frequency (Hz)');ylabel('Power/Frequency (dB/Hz)');grid on;
title('Periodogram Power Spectral Density Estimate');
legend('Serie 1(Blue): STATION','Serie 2(Red): GCM-RCP');
hold off

%% HURST EXPONENT
delete(subplot(2,2,4))
NS = length(STATION_CUM);NG=length(GCM_CUM);
Rs=NaN(1,20);Rg=NaN(1,20);Ss=NaN(1,20);Sg=NaN(1,20);
tau=NaN(1,20);RSs=NaN(1,20);RSg=NaN(1,20);
for i = 1:20;
    m = floor(NS/(2*i));n=floor(NG/(2*i));
    for j=1:i;
        rs = STATION_CUM(1+(j-1)*m:j*m);rg = GCM_CUM(1+(j-1)*n:j*n);
        Ms = mean(rs);Mg = mean(rg);
        xs = (rs-Ms);xg = (rg-Mg);
        Vs = cumsum(xs);Vg = cumsum(xg);
        Rs(1,j) = max(Vs)-min(Vs);Rg(1,j) = max(Vg)-min(Vg);
        Ss(1,j) = std(rs); Sg(1,j) = std(rg);
    end
    tau(1,i) = m;
    RSs(1,i) = nanmean(Rs(1,:)./Ss(1,:),2);RSg(1,i) = nanmean(Rg(1,:)./Sg(1,:),2);
end
subplot(2,2,4);
title(['Hurst Exponent, Interval Accumulation: ' num2str(CHAOS_STATION{1,IAP}) ' days'])
xlabel('log(\tau)','FontSize',12)
ylabel('log(R/S)','FontSize',12)
hold on
qs = polyfit(log10(tau),log10(RSs),1);qg = polyfit(log10(tau),log10(RSg),1);
t = 1:.01:5;
ys = qs(1)*t+qs(2);yg = qg(1)*t+qg(2);
plot(t,ys,'b','LineWidth',2)
text(1.5,2,['y sta =' num2str(qs(1)),' x + ' num2str(qs(2))],'FontSize',12)
plot(t,yg,'r','LineWidth',1)
text(3,1,['y gcm =' num2str(qg(1)),' x + ' num2str(qg(2))],'FontSize',12)
legend('Serie 1 (Azul): Station','Serie 2(Roja): GCM-RCP')
hold off

%% POINCARE MAP
save Precipitation1.dat STATION_CUM -ascii
system([pwd,'\Chaos_Analysis\Auxiliar\cygwin64\poincare.exe ',pwd,'\Precipitation1.dat -m',num2str(cell2mat(CHAOS_STATION(20,IAP))),' -d',num2str(cell2mat(CHAOS_STATION(15,IAP))),' -C1 -o',pwd,'\Poincare1.dat']);
Poincare1=load('Poincare1.dat');
save Precipitation2.dat GCM_CUM -ascii
system([pwd,'\Chaos_Analysis\Auxiliar\cygwin64\poincare.exe ',pwd,'\Precipitation2.dat -m',num2str(cell2mat(CHAOS_STATION(20,IAP))),' -d',num2str(cell2mat(CHAOS_GCM(15,IAP))),' -C1 -o',pwd,'\Poincare2.dat']);
Poincare2=load('Poincare2.dat');
subplot(2,2,3);
plot(Poincare1(:,1),Poincare1(:,2),'b');hold on
plot(Poincare2(:,1),Poincare2(:,2),'r');
title(['Poincaré, Interval Accumulation: ' num2str(CHAOS_STATION{1,IAP}) ' days']);xlabel('t (dias)');ylabel('t- \tau (dias)');
legend(strcat('Serie 1 (Azul): Station \tau :',num2str(cell2mat(CHAOS_STATION(15,IAP))),' dias - m:',num2str(cell2mat(CHAOS_STATION(20,IAP)))),...
strcat('Serie 2(Roja): GCM-RCP \tau :',num2str(cell2mat(CHAOS_GCM(15,IAP))),' dias - m:',num2str(cell2mat(CHAOS_GCM(20,IAP)))));grid on;
hold off

%% CHAOS CHART
h2 = subplot(2, 2, 2);
xl = xlim(h2); xPos = xl(1) + diff(xl)/8 ;yl = ylim(h2);yPos = yl(1) + diff(yl) / 2; 
t = text(xPos, yPos, sprintf('%s\n%s%i%s%i\n%s%i%s%i\n%s%i%s%i\n%s%.2f%s%.2f\n%s%.2f%s%.2f\n%s%.2f%s%.2f',...
    'DINAMYC SYSTEM PROPERTIES         STATION      GCM',...
    'Intervalo de Acumulacion (dias)                    ',CHAOS_STATION{1,IAP},'               ', CHAOS_GCM{1,IAP},...    
    '             Time Delay                               ',CHAOS_STATION{15,IAP},'               ', CHAOS_GCM{15,IAP},...
    '        Embedding Dimension                    ',CHAOS_STATION{20,IAP},'              ', CHAOS_GCM{20,IAP},...
    '    Sum of Lyapunov Exponent          ',CHAOS_STATION{22,IAP},'          ', CHAOS_GCM{22,IAP},...
    '      Maximum Lyapunov Exponent             ',CHAOS_STATION{23,IAP},'           ', CHAOS_GCM{23,IAP},...
    '              Kaplan-Yorke Dimension                   ',CHAOS_STATION{25,IAP},'           ', CHAOS_GCM{25,IAP},...
    ''), 'Parent', h2);
set(t, 'HorizontalAlignment', 'left');
set (h2, 'visible', 'off');
hold off

delete('Poincare1.dat','Poincare2.dat','Precipitation1.dat','Precipitation2.dat')
saveas(gcf,strcat(FOLDER,'\CHAOS_COMP'),'fig');
end