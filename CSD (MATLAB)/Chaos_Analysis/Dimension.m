%% PHASE SPACE RECONSTRUCTION: EMBEDDING DIMENSION
%
% Author: Santiago Duarte (Escuela Colombiana de Ingenieria Julio Garavito), santiagoduare09@gmail.com
% March 2017
%
% INPUTS:
%=============
% PRECIPITATION  : 
% DELAY          : 
%
% OUTPUT:
%=============
% DIM            : 

function [DIM, ENT]=Dimension(DATA, DELAY,FOLDER)

if nargin < 3
    FOLDER = '';
end
if exist('Precipitation.dat', 'file')~=2
    save Precipitation.dat DATA -ascii
end

%MAX DIMENSION
MO=1;
MMAX=floor(log10(length(DATA))*2);

%SUM OF CORRELATION
system([mfilename('fullpath'),'\..\Auxiliar\cygwin64\d2.exe -V0 ',pwd,'\Precipitation.dat -d',num2str(DELAY),' -M',num2str(MO),',',num2str(MMAX),' -o',pwd,'\cor']);
fdat = fopen('cor.c2');
Lineas=cell(2000,1);ln=1;
while ~feof(fdat)
    linea=fgetl(fdat);
    Lineas{ln,1}=linea;
    ln=ln+1;
end
fclose(fdat);
DSUM=zeros(100,2*MMAX);
for p=1:MMAX
    kl=0;b=0;
    c=['#dim= ' ,num2str(p)];
    for ii=1:length(Lineas)
        if strcmp(Lineas{ii,1},c)
            b=Lineas{ii,1};
        end
        if strcmp(b,c)
        kl=kl+1;
            if kl~=1
                DSUM(kl-1,2*p-1:2*p)=str2num(Lineas{ii,1});
            end
            if strcmp(Lineas{ii+1,1},'')
                b=0;
            end
        end
    end
end

%SLOPE OF CORRELATION
fdat2 = fopen('cor.d2','r+');
Lineas=cell(2000,1);ln=1;
while ~feof(fdat2)
    linea=fgetl(fdat2);
    Lineas{ln,1}=linea;
    ln=ln+1;
end
fclose(fdat2);
DSLOPE=zeros(100,2*MMAX);
for p=1:MMAX
    kl=0;b=0;
    c=['#dim= ' ,num2str(p)];
    for ii=1:length(Lineas)
        if strcmp(Lineas{ii,1},c)
            b=Lineas{ii,1};
        end
        if strcmp(b,c)
            kl=kl+1;
            if kl~=1
                DSLOPE(kl-1,2*p-1:2*p)=str2num(Lineas{ii,1});
            end
            if strcmp(Lineas{ii+1,1},'')
                b=0;
            end
        end
    end
end

%CORRELATION DIMENSION
DCOR=NaN(1, MMAX);
for p=1:MMAX
    dcor=NaN(1, 100);
    for x=1:100
        if DSLOPE(x,2*p)~=0
            dcor(x)=DSLOPE(x,2*p);
        end
    end
    DCOR(p)=nanmean(dcor,2);
end
for b=2:length(DCOR)-1
    if (abs(DCOR(b-1)- DCOR(b))/DCOR(b-1))<=0.01
        mdcor=DCOR(b); break
    else
        mdcor=MMAX; 
    end
end
DIM(1,1)=mdcor;

%CORRELATION ENTROPY
fdat3 = fopen('cor.h2','r+');
Lineas=cell(2000,1);ln=1;t=2;a=NaN(MMAX,1);
while ~feof(fdat3)
    linea=fgetl(fdat3);
    if strcmp(linea,['#dim= ' num2str(t)])
        a(t-1)=ln;
        t=t+1;
    end
    if strcmp(linea,'')
        b=ln;
    end
    Lineas{ln,1}=linea;
    ln=ln+1;
end
a(t-1)=b;COR=NaN(1,MMAX);
fclose(fdat3);
for n=1:MMAX
    COR(1,n)=dlmread('cor.h2','',[a(n)-4,1,a(n)-4,1]);
end
p=polyfit(1:MMAX,COR,1);DIM_LIMITE=100;
ENT=p(1)*DIM_LIMITE+p(2);

%FALSE NEAREST NEIGHBORS
MMAX2=2*MMAX+1;
system([mfilename('fullpath'),'\..\Auxiliar\cygwin64\false_nearest.exe -V0 ',pwd,'\Precipitation.dat  -M1,',num2str(MMAX2),' -d',num2str(DELAY),' -o',pwd,'\false_nearest.dat']);
FNN=load('false_nearest.dat') ;
kk=1;
if isempty(FNN)
    mfnn=MMAX2;
else
    while kk<=length(FNN(:,2))
        if FNN(kk,2)<=0.1
            mfnn=FNN(kk,1);
            break
        else
            mfnn=MMAX2;
            kk=kk+1;
        end
    end
end
DIM(2,1)=mfnn;

%CAO'S METHOD
if length(DATA)>(MMAX2*DELAY)
    [E1, E2] = cao_deneme_mod(DATA,DELAY,MMAX2);
    if isnan(E1)
        E1=ones(1,MMAX2-1);
    end
    mcao=1;
    while abs(E1(mcao)-E2(mcao))>0.01*E1(mcao)
        mcao=mcao+1;
        if mcao==MMAX2
            mcao=MMAX2-1;
            break
        end
    end
    mcao=mcao+1;
else
    mcao=mmax;
end
DIM(3,1)=mcao;
if min(DIM)==MMAX2
    DIM(1,1)=2;
end

%SAVE DIMENSION FIGURES (OPTIONAL)
if ~strcmp(FOLDER,'')   
    %PLOT CORRELATION DIMENSION
    figure('visible','off');   
    subplot(2,2,1)
    for p=1:MMAX
        plot(log(DSUM(:,2*p-1)),log(DSUM(:,2*p)));hold all
    end
    title('Sum of Correlation');xlabel('Log(r)');ylabel('Log(C(r))')  
    hold off
    subplot(2,2,2)
    for p=1:MMAX
        plot(log(DSLOPE(:,2*p-1)),DSLOPE(:,2*p));hold all
    end
    title('Slope of Correlation');xlabel('Log(r)');ylabel('Slope \nu ')
    hold off
    subplot(2,2,[3 4])
    plot(DCOR)
    title('Correlation Dimension');xlabel('Embedding dimension');ylabel('Correlation Exponent')
    saveas(gcf,strcat(FOLDER,'Correlation_Dimension'),'jpg'); hold off
    hold off
    
    %PLOT FALSE NEAREST NEIGHBORS
    subplot(1,1,1)
    plot(FNN(:,1),FNN(:,2));
    title('False Nearest Neighbors');xlabel('Embedding dimension');ylabel('% False Nearest Neighbors');
    saveas(gcf,strcat(FOLDER,'False_Nearest_Neighbors'),'jpg');
    
    %PLOT CAO'S METHOD
    plot(1:length(E1),E1,'k'); hold on
    plot(1:length(E2),E2); grid on
    title('Cao''s method'); legend('E1','E2',4);xlabel('Embedding dimension');ylabel('E1-E2');
    saveas(gcf,strcat(FOLDER,'Cao_Method'),'jpg');  
end

%DELETE CREATED FILES
delete('cor.c2','cor.d2','cor.h2','cor.stat','false_nearest.dat');
end