%% STATISTICAL_PROPERTIES_COMPARATION

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

function []=Statistical_Comparation(STATION,GCM, CHAOS_STATION,CHAOS_GCM,CORD_STA,CORD_GCM,IAP,SHAPE,FOLDER)
 
LAT_EST=CORD_STA(1,1);
LONG_EST=CORD_STA(1,2);
LAT_GCM=CORD_GCM(1,:);
LONG_GCM=CORD_GCM(2,:);

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

%TIEMPO VS PRECIPITACION
figure('Visible','off');
subplot(2,2,4)
plot(TIME_CUM,STATION_CUM,'b');hold on
plot(TIME_CUM,GCM_CUM,'r');
title(sprintf('%s%s%s','Time serie for rainfall accumulation of ',mat2str(CHAOS_STATION{1,IAP}), ' dias'));xlabel('Time(days)');ylabel('Precipitation(mm)');
datetick('x','yyyy');legend('Serie 1: Station(blue)','Serie 2: GCM-RCP(red)');grid on;axis([min(TIME_CUM) max(TIME_CUM) 0 max(max(STATION_CUM),max(GCM_CUM))]);
hold off

%LOCATION
subplot(2,2,[1 3])
geoshow(shaperead(SHAPE), 'FaceColor', [1 1 0.8]);
lat=[LAT_GCM(1) LAT_GCM(2) LAT_GCM(2)  LAT_GCM(1)];
lon=[LONG_GCM(1) LONG_GCM(1) LONG_GCM(2)  LONG_GCM(2)];
geoshow(lat,lon, 'DisplayType','polygon','FaceColor', [0.9 0.2 0.2] );
geoshow(LAT_EST, LONG_EST, 'DisplayType', 'Point', 'Marker', '*', 'MarkerEdgeColor','b');
title('Station Location');xlabel('Longitude(°)');ylabel('Latitude(°)');legend('Serie 2: GCM-RCP(red)','Serie 1: Station(blue)')
hold off

%STATISTICAL CHART
h2 = subplot(2, 2, 2);
xl = xlim(h2); xPos = xl(1) + diff(xl) / 4;yl = ylim(h2);yPos = yl(1) + diff(yl) / 3; 
t = text(xPos, yPos, sprintf(...
    '%s\n%s%i%s%i\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f\n%s%.1f%s%.1f',...
    'STATISTICAL PROPERTIES       STATION    GCM',...
    'Interval Accumulation (days)                               :',CHAOS_STATION{1,IAP},',     ', CHAOS_GCM{1,IAP},...
    'N° of data                                     :',CHAOS_STATION{2,IAP},',     ', CHAOS_GCM{2,IAP},...
    'Mean(mm)                                         :',CHAOS_STATION{3,IAP},',     ', CHAOS_GCM{3,IAP},...
    'Median(mm)                                     :',CHAOS_STATION{4,IAP},',     ', CHAOS_GCM{4,IAP},...
    'Maximum(mm)                                     :',CHAOS_STATION{6,IAP},',     ', CHAOS_GCM{6,IAP},...
    'Standard Deviation(mm)                   :',CHAOS_STATION{8,IAP},',     ', CHAOS_GCM{8,IAP},...
    'Variance                                          :',CHAOS_STATION{9,IAP},',     ', CHAOS_GCM{9,IAP},...
    'Coefficient of variation                      :',CHAOS_STATION{10,IAP},',     ', CHAOS_GCM{10,IAP},...
    'Coefficient of asymmetry                      :',CHAOS_STATION{11,IAP},',     ', CHAOS_GCM{11,IAP},...
    'Kurtosis                                               :',CHAOS_STATION{12,IAP},',     ', CHAOS_GCM{12,IAP},...
    ''), 'Parent', h2);
set(t, 'HorizontalAlignment', 'left');
set (h2, 'visible', 'off');
hold off

saveas(gcf,strcat(FOLDER,'\STAT_COMP'),'fig');
end