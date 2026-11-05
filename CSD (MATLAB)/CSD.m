%% CHAOTIC STATISTICAL DOWNSCALING
% Ejecutar el archivo en la carpeta CSD
% Author: Santiago Duarte (Escuela Colombiana de Ingenieria Julio Garavito), santiagoduare09@gmail.com
% May 2017
% Version 1
% INPUTS:
%=============
% Precipitation_Local     :
% Catalog_File_Path       :
% GCM_Hist_File_Path      :
% RCP_Fut_File_Path       :
% Output_Folder           : 
% CSD(EST, CAT, GCM, RCP, OUT)
function [Downscaling,ERROR]=CSD(Precipitation_Station,Station_Catalog,Model_Historic,Model_RCP,Output_Folder)
%% PARAMETERS
% TIME SERIES ACCUMULATION INTERVALS
Int_Acum=[3,5,7,10,15,30];

% CALIBRATION, VALIDATION AND APPLICATION DATES
Cal_dates= datetime([1850,2005],[1,12],[1,31]);
Val_dates= datetime([2006,2016],[1,8],[1,31]);
App_dates= datetime([2006,2100],[1,12],[1,31]);

% SHAPE
SHAPE='Zonificacion_hidrografica_2013.shp';

%% LOAD MODEL DATASET (Historic and RCP) (.nc)
setup_nctoolbox()
gh= ncdataset(Model_Historic);
rcp= ncdataset(Model_RCP);
GCM_hist_time = gh.time('time')-0.5;% Read GCM time data and convert to MatLab format
RCP_fut_time = rcp.time('time')-0.5;% Read RCP time data and convert to MatLab format
model_latitude_array= gh.data('lat_bnds');% Read GCM latitude data box
model_longitude_array = gh.data('lon_bnds');% Read GCM longitude data box

%% LOAD OBSERVATION DATASET (.xlsx)
[~,~,raw]=xlsread(Precipitation_Station);  % Read the xls file
NS=length(raw(1,:));  % Number of stations to read (all default)
Station=NaN(length(raw(:,1))-1,length(raw(1,:)));
Station(:,1)=datenum(raw(2:end,1),'dd/mm/yyyy');
Station(:,2:end)=cell2mat(raw(2:end,2:end));
resume=NaN(NS-1,3);

%% COORDINATES:LOCAL VARIABLE (.xlsx)
[~,~,loc]=xlsread(Station_Catalog);
[rowStation,~] = find(strcmp(loc(2:end,1), char(raw(1,2:NS))));
Location=cell2mat(loc(rowStation+1,2:4));
station_latitude=Location(:,1); 
station_longitude=Location(:,2);
lon_change= station_longitude<0;
station_longitude(lon_change)=station_longitude+360;

for z=2:2 
    %% CREATE MAIN DIRECTORY(STATION)
    SFolder=strcat(Output_Folder,char(raw(1,z)));   
    mkdir(SFolder);
    
    %% FILL DATA
    [Data,LM,LE]=Fill_Data(Station,Location,z);     
    fprintf(strcat(' STATION: \t',char(raw(1,z)),' \n Missing data:\t',num2str(LM),'\n Existent data:\t',num2str(LE),'\n Total data:\t',num2str(LM+LE)))
    save(strcat(SFolder,'\Data.mat'),'Data')
    
    %% CHAOS ANALYSIS:LOCAL VARIABLE   
    if exist(strcat(SFolder,'\Chaos.mat'),'file')==0
        [Station_Chaos]=Chaos_Analysis(Data(:,2),Int_Acum,SFolder);
        save(strcat(SFolder,'\Chaos.mat'),'Station_Chaos')
    else
        load(strcat(SFolder,'\Chaos.mat'))
    end
    
    %% GET GCM CELL
    GCM_latitude=1;GCM_longitude=1;
    while GCM_latitude<=length(model_latitude_array)
        if model_latitude_array(GCM_latitude,1)<=station_latitude(z-1) && station_latitude(z-1)<=model_latitude_array(GCM_latitude,2)
            break
        end
        GCM_latitude=GCM_latitude+1;
    end
    while GCM_longitude<=length(model_longitude_array)
        if model_longitude_array(GCM_longitude,1)<=station_longitude(z-1) && station_longitude(z-1)<=model_longitude_array(GCM_longitude,2)
            break
        end
        GCM_longitude=GCM_longitude+1;
    end
    
    %% CREATE MAIN DIRECTORY(MODEL)
    MFolder=strcat(Output_Folder,'GCM_Lat_',num2str(model_latitude_array(GCM_latitude,1)),'_Long_',num2str(model_longitude_array(GCM_longitude,1)));   
    mkdir(MFolder);
   
    %% CHAOS ANALYSIS:CLIMATE MODEL   
    if exist(strcat(MFolder,'\Chaos.mat'),'file')==0
        % GET GCM DATA
        historic_model=gh.data('pr', [1 GCM_latitude GCM_longitude], [length(GCM_hist_time) GCM_latitude GCM_longitude])*86400; % yyy*86400 = Conversion from kg/ m^2 / s to mm/day
        future_model=rcp.data('pr', [1 GCM_latitude GCM_longitude], [length(RCP_fut_time) GCM_latitude GCM_longitude])*86400; % yyy*86400 = Conversion from kg/ m^2 / s to mm/day
        
        % ENSAMBLE MODEL DATA
        Model(1:length(historic_model),1)=GCM_hist_time;
        Model(1:length(historic_model),2)=historic_model;
        Model((length(historic_model)+1):(length(historic_model)+length(future_model)),1)=RCP_fut_time;
        Model((length(historic_model)+1):(length(historic_model)+length(future_model)),2)=future_model;
        
        [Model_Chaos]=Chaos_Analysis(Model(:,2),Int_Acum,MFolder);
        save(strcat(MFolder,'\Chaos.mat'),'Model_Chaos')
        save(strcat(MFolder,'\Data.mat'),'Model')
    else
        load(strcat(MFolder,'\Chaos.mat'))
        load(strcat(MFolder,'\Data.mat'))
    end
    
    %% APPROPIATE INTERVAL ACCUMULATION
    IA=NaN(1,1);iap=2;
    while isnan(IA)&& iap<=length(Int_Acum)+1
        if strcmp(Station_Chaos{24,iap},'Deterministic_Chaos') && strcmp(Model_Chaos{24,iap},'Deterministic_Chaos')
            IA=iap;
        end    
        iap=iap+1;
    end
    if isnan(IA)
        display('Deterministic chaos was not found in the predictor and predictand for the same rainfall accumulation)')     
    else
        if exist(strcat(Output_Folder,'resume.mat'),'file')==1
            load(strcat(Output_Folder,'resume.mat'),'resume')
        end
        resume(z-1,1)=IA; % APPROPIATE INTERVAL ACCUMULATION
        resume(z-1,2)=Station_Chaos{15,IA}; % KAPLAN-YORKE DIMENSION STATION
        resume(z-1,3)=Station_Chaos{20,IA}; % EMBEDDING DIMENSION STATION
        resume(z-1,4)=Station_Chaos{25,IA}; % KAPLAN-YORKE DIMENSION STATION
        resume(z-1,5)=Model_Chaos{15,IA};  % TIME DELAY DIMENSION GCM
        resume(z-1,6)=Model_Chaos{20,IA}; % EMBEDDING DIMENSION STATION
        resume(z-1,7)=Model_Chaos{25,IA}; % KAPLAN-YORKE DIMENSION STATION
        save(strcat(Output_Folder,'\resume.mat'),'resume')
    end
    
    %% COMPARE HISTORIC STATISTICAL PROPERTIES
    CORD_STA=[station_latitude,station_longitude];    
    CORD_GCM(1,:)=model_latitude_array(GCM_latitude,:);
    CORD_GCM(2,:)=model_longitude_array(GCM_longitude,:);   
    Statistical_Comparation(Data,Model,Station_Chaos,Model_Chaos,CORD_STA,CORD_GCM,IA,SHAPE ,SFolder)
 
    %% COMPARE HISTORIC CHAOTIC PROPERTIES
    Chaos_Comparation(Data,Model,Station_Chaos,Model_Chaos,IA,SFolder)      

    %% CREATE RESULT VARIABLE
    if ~exist('Downscaling','var')
        INI=find(Model(:,1)==datenum(Val_dates(2))+1);
        Downscaling=NaN(length(Model(INI:end,1)),NS);
        Downscaling(:,1)=Model(INI:end,1);   
    end
    
    %% CALIBRATION
    [CALB]=Calibration(Data,Model,Station_Chaos,Model_Chaos,IA,Cal_dates,SFolder);

    %% VALIDATION
    [ERROR]=Validation(Data,Model,Station_Chaos,Model_Chaos,IA,Val_dates,CALB,SFolder);
    
    %% SINCHRONYZATION
    [SYNC]=Synchronization(Data,Model,Station_Chaos,Model_Chaos,IA,App_dates,CALB,SFolder);
    save(strcat(SFolder,'\Syncronization.mat'),'SYNC')
    xlswrite(strcat(SFolder,'\Syncronization.xls'),SYNC)
    IN=find(SYNC(:,1)==Downscaling(:,1));
    Downscaling(:,z)=SYNC(IN:end,2);    
end
end