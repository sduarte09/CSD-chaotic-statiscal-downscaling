%% CHAOS ANALYSIS OF TIME SERIES
% Author: Santiago Duarte (Escuela Colombiana de Ingenieria Julio Garavito), santiagoduare09@gmail.com
% May 2017
% Version 1
%INPUTS:
%=============
% DATA                   :
% INTERVAL_ACCUMULATION  :
%
% OUTPUT:
%=============
% CHAOS_INFO: 

function [CHAOS_INFO]=Chaos_Analysis(DATA ,INTERVAL_ACCUMULATION,FOLDER)
%% Create the output variable of chaos analysis results
CHAOS_INFO=cell(27,length(INTERVAL_ACCUMULATION)+1);  
CHAOS_INFO{1,1}='Interval Accumulation';
CHAOS_INFO{2,1}='N° of data accumulated';
CHAOS_INFO{3,1}='Mean';
CHAOS_INFO{4,1}='Median';
CHAOS_INFO{5,1}='Minimum';
CHAOS_INFO{6,1}='Maximum';
CHAOS_INFO{7,1}='Range';
CHAOS_INFO{8,1}='Standard Deviation';
CHAOS_INFO{9,1}='Variance';
CHAOS_INFO{10,1}='Coefficient of variation'; 
CHAOS_INFO{11,1}='Coefficient of asymmetry';
CHAOS_INFO{12,1}='Kurtosis';
CHAOS_INFO{13,1}='Time Delay Autocorrelation';
CHAOS_INFO{14,1}='Time delay Mutual Information';
CHAOS_INFO{15,1}='Time Delay';
CHAOS_INFO{16,1}='Correlation Dimension';
CHAOS_INFO{17,1}='False Nearest Neighbors Dimension';
CHAOS_INFO{18,1}='Cao''s Dimension';
CHAOS_INFO{19,1}='Embedding Dimension';
CHAOS_INFO{20,1}='Embedding Dimension Verified';
CHAOS_INFO{21,1}='Lyapunov Exponents';
CHAOS_INFO{22,1}='Sum of Lyapunov Exponent';
CHAOS_INFO{23,1}='Maximum Lyapunov Exponent';
CHAOS_INFO{24,1}='Type of motion of Dynamical System';
CHAOS_INFO{25,1}='Kaplan-Yorke Dimension';
CHAOS_INFO{26,1}='Correlation Entropy';
CHAOS_INFO{27,1}='Pesin''s identity';

for p=INTERVAL_ACCUMULATION(1:end)    
    %% ANALYSIS OF DATA
    FOLDER_STA=strcat(FOLDER,'\',num2str(p),'_day_');
    k=1;Info=NaN(length(DATA),1);
    while k<=(length(DATA)-p)
        if isnan(DATA(k:k+p-1))
            Info(k)=NaN;
        else
            Info(k)=nansum(DATA(k:k+p-1));% Transforming data to NaN and Numbers
        end
        k=k+p;
    end
    a=isnan(Info);b=sum(sum(a));c=length(a)-b;
    Precipitation=Info(1:p:end); % Accumlated data extraction
    if isnan(Precipitation(end))
        Precipitation(end)=[];
    end

    %% PHASE SPACE RECONSTRUCTION: TIME DELAY
    [lag]=Delay(Precipitation,FOLDER_STA);
    
    %% PHASE SPACE RECONSTRUCTION: EMBEDDING DIMENSION
    [dim, entropy]=Dimension(Precipitation, min(lag),FOLDER_STA);
    
    %% LYAPUNOV SPECTRUM
    [Lambda, ky,dimv]=Lyapunov(Precipitation, min(lag),min(dim),FOLDER_STA);
    if sum(Lambda)<0 && max(Lambda)>0
        TYPE_MOTION='Deterministic_Chaos';
    elseif sum(Lambda)>0 && max(Lambda)>0
        TYPE_MOTION='Noise_(Random)';
    elseif sum(Lambda)<0 && max(Lambda)==0
        TYPE_MOTION='Stable_Fixed_Point';
    else
        TYPE_MOTION='Stable_Limit_Cycle';
    end
    %% WRITE RESULTS
    CHAOS_INFO{1,find(INTERVAL_ACCUMULATION==p)+1}=p;
    CHAOS_INFO{2,find(INTERVAL_ACCUMULATION==p)+1}=c;
    CHAOS_INFO{3,find(INTERVAL_ACCUMULATION==p)+1}=mean(Precipitation);
    CHAOS_INFO{4,find(INTERVAL_ACCUMULATION==p)+1}=median(Precipitation);
    CHAOS_INFO{5,find(INTERVAL_ACCUMULATION==p)+1}=min(Precipitation);
    CHAOS_INFO{6,find(INTERVAL_ACCUMULATION==p)+1}=max(Precipitation);
    CHAOS_INFO{7,find(INTERVAL_ACCUMULATION==p)+1}=range(Precipitation);
    CHAOS_INFO{8,find(INTERVAL_ACCUMULATION==p)+1}=std(Precipitation);
    CHAOS_INFO{9,find(INTERVAL_ACCUMULATION==p)+1}=var(Precipitation);
    CHAOS_INFO{10,find(INTERVAL_ACCUMULATION==p)+1}=std(Precipitation)/mean(Precipitation);
    CHAOS_INFO{11,find(INTERVAL_ACCUMULATION==p)+1}=skewness(Precipitation);
    CHAOS_INFO{12,find(INTERVAL_ACCUMULATION==p)+1}=kurtosis(Precipitation);
    CHAOS_INFO{13,find(INTERVAL_ACCUMULATION==p)+1}=lag(1);
    CHAOS_INFO{14,find(INTERVAL_ACCUMULATION==p)+1}=lag(2);
    CHAOS_INFO{15,find(INTERVAL_ACCUMULATION==p)+1}=min(lag);
    CHAOS_INFO{16,find(INTERVAL_ACCUMULATION==p)+1}=dim(1);
    CHAOS_INFO{17,find(INTERVAL_ACCUMULATION==p)+1}=dim(2);
    CHAOS_INFO{18,find(INTERVAL_ACCUMULATION==p)+1}=dim(3);
    CHAOS_INFO{19,find(INTERVAL_ACCUMULATION==p)+1}=min(dim);
    CHAOS_INFO{20,find(INTERVAL_ACCUMULATION==p)+1}=dimv;
    CHAOS_INFO{21,find(INTERVAL_ACCUMULATION==p)+1}=Lambda;
    CHAOS_INFO{22,find(INTERVAL_ACCUMULATION==p)+1}=sum(Lambda);
    CHAOS_INFO{23,find(INTERVAL_ACCUMULATION==p)+1}=max(Lambda);
    CHAOS_INFO{24,find(INTERVAL_ACCUMULATION==p)+1}=TYPE_MOTION;
    CHAOS_INFO{25,find(INTERVAL_ACCUMULATION==p)+1}=ky;
    CHAOS_INFO{26,find(INTERVAL_ACCUMULATION==p)+1}=entropy;
    CHAOS_INFO{27,find(INTERVAL_ACCUMULATION==p)+1}=sum(Lambda(Lambda>0));  
end
end