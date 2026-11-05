%% Calibration

% Author: Santiago Duarte (Escuela Colombiana de Ingenieria Julio Garavito), santiagoduare09@gmail.com
% May 2017

%INPUTS:
%=============
% EST                :
% GCM                :
% CHAOS_EST          :
% CHAOS_GCM          :
% INTER_ACUM_POS     : 
% LIM_INF            :
% LIM_SUP            :
%
%% OUTPUT:
%=============
% DS : 

function [CALB]=Calibration(STATION,GCM,CHAOS_STATION,CHAOS_GCM,IAP,CAL_DATES,FOLDER)

CAL_END=datenum(CAL_DATES(2));
INI=find(GCM(:,1)==STATION(1,1));
step=CHAOS_GCM{1,IAP};

%% DELAY VECTORS
Dimension_STA=cell2mat(CHAOS_STATION(20,IAP));
Dimension_GCM=cell2mat(CHAOS_GCM(20,IAP));
Delay_STA=cell2mat(CHAOS_STATION(15,IAP));
Delay_GCM=cell2mat(CHAOS_GCM(15,IAP));

%ACUMMULATED STATION CALIBRATION
if CAL_END>STATION(end,1)
   CAL_END=floor(length(STATION(:,1))*0.8)+STATION(1,1);
end
STA_CAL=STATION(1:find(STATION(:,1)==CAL_END),2);
k=1;Info_CAL=NaN(length(STA_CAL),1);
while k<=(length(STA_CAL)-step)
    if isnan(STA_CAL(k:k+step-1))
        Info_CAL(k)=NaN;
    else
        Info_CAL(k)=nansum(STA_CAL(k:k+step-1));
    end
    k=k+step;
end
STA_CAL_CUM=Info_CAL(1:step:end);
if isnan(STA_CAL_CUM(end))
    STA_CAL_CUM(end)=[];
end

%ACUMMULATED GCM CALIBRATION
GCM_CAL=GCM(INI:find(GCM(:,1)==CAL_END),2);
k=1;Info=NaN(length(GCM_CAL),1);
while k<=(length(GCM_CAL)-step)
    if isnan(GCM_CAL(k:k+step-1))
        Info(k)=NaN;
    else
        Info(k)=nansum(GCM_CAL(k:k+step-1));
    end
    k=k+step;
end
GCM_CAL_CUM=Info(1:step:end); 
if isnan(GCM_CAL_CUM(end))
    GCM_CAL_CUM(end)=[];
end

%CREATE DELAY VECTORS
STACALVEC=NaN(length(STA_CAL_CUM)-((Dimension_STA-1)*Delay_STA),Dimension_STA);
GCMCALVEC=NaN(length(GCM_CAL_CUM)-((Dimension_GCM-1)*Delay_GCM),Dimension_GCM);
for i=1:Dimension_STA
    STACALVEC(:,i)=STA_CAL_CUM(1+Delay_STA*(Dimension_STA-i):length(STA_CAL_CUM)-Delay_STA*(i-1));
end
for i=1:Dimension_GCM
    GCMCALVEC(:,i)=GCM_CAL_CUM(1+Delay_GCM*(Dimension_GCM-i):length(GCM_CAL_CUM)-Delay_GCM*(i-1));
end

%% PARAMETERS: LAG AND AMPLITUDE (Generalized Similarity Function)
fprintf('%s\n','CALIBRATION: LAG AND AMPLITUDE');
LEN=min(length(STACALVEC),length(GCMCALVEC));
lag_c=1;Sg_mat=NaN(100,(floor(30/step)));
AM(1)=max(STA_CAL_CUM)/max(GCM_CAL_CUM);
AM(2)=mean(STA_CAL_CUM)/mean(GCM_CAL_CUM);
LAG_MAT=0:(floor(30/step)-1); AMP_MAT=0.1:(max(AM)-min(AM))/99:max(AM);

% COMPUTE AMP AND LAG
for LAG=LAG_MAT(1:end)
    am_c=1;
    for AMP=AMP_MAT(1:end)         
        STACALVNLAG=STACALVEC(1+LAG:LEN)-repmat(mean(STACALVEC(1+LAG: LEN)),1,LEN-LAG); % Zero mean + LAG
        STACALVN=STACALVEC(1:LEN-LAG)-repmat(mean(STACALVEC(1:LEN-LAG)),1,LEN-LAG); % Zero mean
        GCMCALVN=GCMCALVEC(1:LEN-LAG)-repmat(mean(GCMCALVEC(1:LEN-LAG)),1,LEN-LAG); % Zero mean
        Sg_mat(am_c,lag_c)=(((STACALVNLAG-AMP*GCMCALVN).^2/(LEN-LAG))/(((STACALVN.^2)/(LEN-LAG).*(AMP*GCMCALVN.^2)/(LEN-LAG)).^(1/2)))^(1/2);       
        am_c=am_c+1;
    end
    lag_c=lag_c+1;
end
[row,col,Sg] = find(min(Sg_mat(:))==Sg_mat);
LAG=LAG_MAT(col); AMP=AMP_MAT(row);
        
%LAG CORRECTION
STACALV=STACALVEC(1+LAG:end,:); GCMCALV=GCMCALVEC(1:end-LAG,:);
LENCALV=min(length(GCMCALV),length(STACALV));
STACALV=STACALV(1:LENCALV,:); GCMCALV=GCMCALVEC(1:LENCALV,:);

%AMPLITUDE CORRECTION
GCMCALV=GCMCALV*AMP;

%% PARAMETERS: U AND SS  Mutual False Nearest Neighbors
fprintf('%s\n','CALIBRATION: U');
ss=NaN(LENCALV,1);u=NaN(LENCALV,1);

% FIND U AND EPSILON
for t=1:LENCALV
    % DISTANCE IN STATION
    distanceA = pdist2(GCMCALV(1:LENCALV,:),GCMCALV(t,:));
    [nearA, nposA]=sort(distanceA);
    l=find(nearA~=0,1,'first');
    Dist_GCM_GCM=nearA(l);
    pas(1,:)=STACALV(t,:);pas(2,:)=STACALV(nposA(l),:);
    Dist_Est_GCM=pdist(pas);
    % DISTANCE IN GCM
    distanceB = pdist2(STACALV,STACALV(t,:));
    [nearB, nposB]=sort(distanceB);
    k=find(nearB~=0,1,'first');
    Dist_Est_Est=nearB(k);
    pes(1,:)=GCMCALV(t,:);pes(2,:)=GCMCALV(nposB(k),:);
    Dist_GCM_Est=pdist(pes);
    % COMPUTE SS, U AND EPSILON
    ss(t,1)=(Dist_Est_GCM*Dist_GCM_Est)/(Dist_GCM_GCM*Dist_Est_Est);
    u(t,1)=sum(ss(1:t))/t;
    clear('pas','pes')   
end

%WRITE VARIABLES
CALB{1}=LAG;
CALB{2}=AMP;
CALB{3}=ss;
CALB{4}=u;
CALB{5}=Sg;

%PLOT MFNN
figure('Visible','off');
plot(u,'r');hold on
title(sprintf('%s%s%s','Mutual False Nearest Neighbors, rainfall accumulation: ',mat2str(CHAOS_STATION{1,IAP}), ' days'));
xlabel('Time(days)');ylabel('Parameter \mu');
grid on;hold off
saveas(gcf,strcat(FOLDER,'\MFNN'),'jpg');

%PLOT Sg
AMP_MAT1=repmat(AMP_MAT',1,length(LAG_MAT));
LAG_MAT1=repmat(LAG_MAT,length(AMP_MAT),1);
surf(AMP_MAT1,LAG_MAT1,Sg_mat);hold on
title(sprintf('%s%s%s','Generalized Similarity Function, rainfall accumulation: ',mat2str(CHAOS_STATION{1,IAP}), ' days'));
xlabel('AMPLITUDE');ylabel('LAG');zlabel('Sg');
grid on;hold off
saveas(gcf,strcat(FOLDER,'\Similarity'),'jpg');
fprintf('%s\n','CALIBRATION COMPLETE');

end