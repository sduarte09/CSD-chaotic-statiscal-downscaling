%% Validation

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

function [ERROR]=Validation(STATION,GCM,CHAOS_STATION,CHAOS_GCM,IAP,APP_DATES,CALB,FOLDER)

APP_INI=datenum(APP_DATES(1));APP_FIN=datenum(APP_DATES(2));
INI=find(GCM(:,1)==STATION(1,1));
step=CHAOS_GCM{1,IAP};

%% DELAY VECTORS
Dimension_STA=cell2mat(CHAOS_STATION(20,IAP));
Dimension_GCM=cell2mat(CHAOS_GCM(20,IAP));
Delay_STA=cell2mat(CHAOS_STATION(15,IAP));
Delay_GCM=cell2mat(CHAOS_GCM(15,IAP));

%ACUMMULATED STATION CALIBRATION
if APP_INI>STATION(end,1)
   APP_INI=floor(length(STATION(:,1))*0.8)+STATION(1,1)+1;
end
STA_CAL=STATION(1:(find(STATION(:,1)==APP_INI)-1),2);
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

%ACUMMULATED STATION APPLICATION
DIF=sum(isnan(Info_CAL(end-step:end)))+1-step;
STA_VAL=STATION((find(STATION(:,1)==APP_INI)-DIF):end,2);
k=1;Info_VAL=NaN(length(STA_VAL),1);
while k<=(length(STA_VAL)-step)
    if isnan(STA_VAL(k:k+step-1))
        Info_VAL(k)=NaN;
    else
        Info_VAL(k)=nansum(STA_VAL(k:k+step-1));
    end
    k=k+step;
end
STA_VAL_CUM=Info_VAL(1:step:end);
if isnan(STA_VAL_CUM(end))
    STA_VAL_CUM(end)=[];
end

%ACUMMULATED GCM CALIBRATION
GCM_CAL=GCM(INI:(find(GCM(:,1)==APP_INI)-1),2);
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

%ACUMMULATED GCM APPLICATION
GCM_VAL=GCM(find(GCM(:,1)==APP_INI)-DIF:find(GCM(:,1)==APP_FIN),2);%VALIDATION GCM TIME SERIES
k=1;Info=NaN(length(GCM_VAL),1);
while k<=(length(GCM_VAL)-step)
    if isnan(GCM_VAL(k:k+step-1))
        Info(k)=NaN;
    else
        Info(k)=nansum(GCM_VAL(k:k+step-1));
    end
    k=k+step;
end
GCM_VAL_CUM=Info(1:step:end); 
if isnan(GCM_VAL_CUM(end))
    GCM_VAL_CUM(end)=[];
end

%ACUMMULATED TIME
TIME_CUM=GCM(find(GCM(:,1)==APP_INI):CHAOS_STATION{1,IAP}:find(GCM(:,1)==APP_FIN),1); 
if length(TIME_CUM)>length(GCM_CAL_CUM)+length(GCM_VAL_CUM)
    TIME_CUM(end)=[];
end

%CREATE DELAY VECTORS
STACALVEC=NaN(length(STA_CAL_CUM)-((Dimension_STA-1)*Delay_STA),Dimension_STA);
GCMCALVEC=NaN(length(GCM_CAL_CUM)-((Dimension_GCM-1)*Delay_GCM),Dimension_GCM);
GCMVALVEC=NaN(length(GCM_VAL_CUM)-((Dimension_GCM-1)*Delay_GCM),Dimension_GCM);
for i=1:Dimension_STA
    STACALVEC(:,i)=STA_CAL_CUM(1+Delay_STA*(Dimension_STA-i):length(STA_CAL_CUM)-Delay_STA*(i-1));
end
for i=1:Dimension_GCM
    GCMCALVEC(:,i)=GCM_CAL_CUM(1+Delay_GCM*(Dimension_GCM-i):length(GCM_CAL_CUM)-Delay_GCM*(i-1));
    GCMVALVEC(:,i)=GCM_VAL_CUM(1+Delay_GCM*(Dimension_GCM-i):length(GCM_VAL_CUM)-Delay_GCM*(i-1));
end

%% VALIDATION
fprintf('%s\n','VALIDATION IN PROGRESS');

%PARAMETERS:
LAG=CALB{1};
AMP=CALB{2};
ss=CALB{3};
u=CALB{4};

%LAG CORRECTION
STACALV=STACALVEC(1+LAG:end,:);
GCMCALV=GCMCALVEC(1:end-LAG,:);
if LAG==0
    GCMVALV=GCMVALVEC;
else
    GCMVALV(1:LAG,:)=GCMCALVEC(end-LAG+1:end,:);
    GCMVALV(LAG+1:length(GCMVALVEC),:)=GCMVALVEC(1:end-LAG,:);
end
LENCALV=min(length(GCMCALV),length(STACALV));
STACALV=STACALV(1:LENCALV,:); GCMCALV=GCMCALVEC(1:LENCALV,:);

% AMPLITUDE CORRECTION
GCMCALV=GCMCALV*AMP;

% MUTUAL FALSE NEAREST NEIGHBORS
t=length(ss);
N=30;
LENGCMCALV=length(GCMVALV);

% FIND N AND VERIFY CALIBRATION
for v=1:LENGCMCALV    
    %FIND Xn, XnNND AND YnNND (Direct)
    GCMCALV((end+1),:)=GCMVALV(v,:);
    EE=repmat(GCMVALV(v,:),length(GCMCALV),1);
    distanceC= sqrt(sum((EE-GCMCALV).^2,2));            
    [nearC, nposC]=sort(distanceC);
    lm=find(nearC~=0,1,'first');
    Xn(1,:)=GCMVALV(v,:);
    XnNND=GCMCALV(nposC(lm),:);
    YnNND=STACALV(nposC(lm),:); 
    %FIND Yn, YnNNR AND XnNNR (Iteration)
    DISX=pdist2(Xn,XnNND)./pdist2(Xn,GCMCALV(1:(t+v-1),:));
    mulYnNND=repmat(YnNND,length(STACALV),1);
    C(1,2:Dimension_STA)=STACALV(length(STACALV)-Delay_STA+1,1:(Dimension_STA-1));
    C=repmat(C,length(STACALV),1); C1=C; C2=C;
    c2d=-2*sum((C(:,2:end).*mulYnNND(:,2:end)),2)+ sum(mulYnNND.^2,2);      
    %FUNCTION U
    DIST_END=100000;
    sso=nanmean(ss(floor(length(ss)/2):end));
    epsilon=nanstd(ss(floor(length(ss)/2):end)); 
    if sso-epsilon>0
        FSmin=sso-epsilon;
    else
        FSmin=0;
    end   
    FSmax=sso+epsilon;
    FS_EPS=FSmin:(FSmax-FSmin)/N:FSmax;
    for FS=FS_EPS(1:end)       
        FACTOR=(FS*DISX').^2;       
        a2=1-FACTOR;
        b2= 2*(FACTOR).*STACALV(:,1) -2*mulYnNND(:,1);    
        c2n=(1-FACTOR).*sum((C(:,2:end)).^2,2);   
        c2r=(2*FACTOR).*sum((C(:,2:end).*STACALV(:,2:end)),2)-(FACTOR).*sum(STACALV.^2,2);
        c2=c2d+c2r+c2n;      
        c31=(-b2+sqrt(b2.^2-4*a2.*c2))./(2*a2);
        c31(imag(c31)~=0)=NaN;c31(c31<0)=NaN;C1(:,1)=c31;
        distanceD1= sqrt(sum((STACALV-C1).^2,2));
        mindist1=min(distanceD1);
        c32=(-b2-sqrt(b2.^2-4*a2.*c2))./(2*a2);
        c32(imag(c32)~=0)=NaN;c32(c32<0)=NaN;C2(:,1)=c32;
        distanceD2= sqrt(sum((STACALV-C2).^2,2));  
        mindist2=min(distanceD2);
        mindist=min([mindist1 mindist2]);
        if mindist1<mindist2
            distanceD=distanceD1;
            CEND=c31;
        else
            distanceD=distanceD2;
            CEND=c32;
        end
        if mindist<DIST_END 
           EPS=distanceD;
           C(:,1)=CEND;
           DIST_END =mindist;
        end
    end  
    %WRITE Yn, YnNNR AND XnNNR
    if ~exist('EPS','var')
        C(:,1)=STACALV(:,1);
        EPS=sqrt(sum((STACALV-C).^2,2));
    end
    [nearD, nposD]=sort(EPS);
    m=find(nearD~=0,1,'first');
    Yn=C(nposD(m),:);
    XnNNR=GCMCALV(nposD(m),:);
    YnNNR=STACALV(nposD(m),:);
    ss(t+v)=(pdist2(Yn,YnNND)*pdist2(Xn,XnNNR))/(pdist2(Yn,YnNNR)*pdist2(Xn,XnNND));
    u(t+v)=sum(ss(1:(t+v)))/(t+v);
    STACALV((end+1),:)=Yn;
    clear C C1 C2
end

%SERIE PRECIPITATION RECONSTRUCTION
DSC(:,1)=STACALV(:,end);
DSC(end+1:end+Delay_STA*(Dimension_STA-1),1)=STACALV(end-Delay_STA*(Dimension_STA-1)+1:end,1);
VAL=DSC(LENCALV+1:end,1);

%COMPUTE ERROR
delta1=sqrt(sum((STA_CAL_CUM(1:length(VAL))-VAL).^2)/length(VAL));
delta2=abs(max(STA_CAL_CUM(1:length(VAL)))-max(VAL));       
VALIDATION(:,1)=TIME_CUM(1:length(VAL));
VALIDATION(:,2)=VAL;

%WRITE VARIABLES
ERROR(1)=delta1;
ERROR(2)=delta2;

%PLOT TIME SERIES
figure('Visible','off');
plot(VALIDATION(:,1),VALIDATION(:,2),'b');hold on
plot(TIME_CUM(1:length(STA_VAL_CUM),1),STA_VAL_CUM,'r');
title(sprintf('%s%s%s','TIME SERIES, RAINFALL ACCUMULATION: ',mat2str(CHAOS_STATION{1,IAP}), ' days'));
xlabel('Time(days)');ylabel('Precipitation(mm)');
datetick('x','yyyy');legend('Serie 1: STATION (MODEL)(blue)','Serie 2: STATION (OBSEVATIONS)(red)');
grid on;axis([min(TIME_CUM) max(TIME_CUM) 0 max(max(STA_VAL_CUM),max(VALIDATION(:,2)))]);
saveas(gcf,strcat(FOLDER,'\VALIDATION'),'jpg');
fprintf('%s\n','VALIDATION COMPLETE');
end