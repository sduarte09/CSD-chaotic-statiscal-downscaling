%% LYAPUNOV SPECTRUM
%
% Author: Santiago Duarte (Escuela Colombiana de Ingenieria Julio Garavito), santiagoduare09@gmail.com
% November 2016
%
% INPUTS:
%=============
% PRECIPITATION  :
% DIMENSION      :
% DELAY          : 
%
% OUTPUT:
%=============
% LAMBDA : 
% KY_DIM : 
% VERIFIED_DIMENSION : 
% TYPE_MOTION: 

function [LAMBDA,KY_DIM,VERIFIED_DIMENSION]=Lyapunov(DATA,DELAY, DIMENSION,FOLDER)

if nargin < 4
    FOLDER = '';
end
if exist('Precipitation.dat', 'file')~=2
    save Precipitation.dat DATA -ascii
end

%LYAPUNOV ESPECTRUM
verif='';comp=0;
while (strcmp(verif,'')==1 && DIMENSION<=24) || comp==0
    system([mfilename('fullpath'),'\..\Auxiliar\cygwin64\delay.exe -V0 ',pwd,'\Precipitation.dat -d',num2str(DELAY),' -m',num2str(DIMENSION),' -o',pwd,'\delvec.dat']);
    system([mfilename('fullpath'),'\..\Auxiliar\cygwin64\lyap_spec.exe -V0 ',pwd,'\delvec.dat -c1,',num2str(DIMENSION),' -m1,',num2str(DIMENSION),' -I -o',pwd,'\lyap_spec.dat']);
    fod=fopen('lyap_spec.dat','r+');
    verif=fscanf(fod,'%c',1);
    linea=fgetl(fod);
    comp=isempty(strfind(linea, '-1.#IND00e+000'));
    fclose(fod);
    DIMENSION =DIMENSION+1;
end
DIMENSION =DIMENSION-1;
system([mfilename('fullpath'),'\..\Auxiliar\cygwin64\lyap_spec.exe -V0 ',pwd,'\delvec.dat -c1,',num2str(DIMENSION),' -m1,',num2str(DIMENSION),' -I -o',pwd,'\lyap_spec.dat']);
fod = fopen('lyap_spec.dat','r+');
Lineas=cell(200,1);ln=1;nm=1;
while ~feof(fod)
    linea=fgetl(fod);
    Lineas{ln,1}=linea;
    ln=ln+1;
end
fclose(fod);
while length(Lineas{nm,1})~=50 || length(Lineas{nm+4,1})==41
    nm=nm+1;
    klp=0;
    if isempty(Lineas{nm,1})
        klp=99;
        break
    end
end
if length(linea)==33
    Lys=dlmread('lyap_spec.dat',' ',[0 0 nm-2 DIMENSION]);   
    LAMBDA=Lys(nm-1,2:end);
else
    klp=99; LAMBDA=99;Lys=0;
end
if klp~=99
    KY_DIM=dlmread('lyap_spec.dat',' ',[nm+3 2 nm+3 2]);
else
    KY_DIM=99;
end
VERIFIED_DIMENSION=DIMENSION;

%SAVE LYAPUNOV FIGURES (OPTIONAL)
if ~strcmp(FOLDER,'') &&  length(Lys(:,1))>=2
    %LYAPUNOV SPECTRUM
    figure('Visible','off');
    for s=2:DIMENSION+1
        plot(Lys(:,1),Lys(:,s)); hold all
    end
    leg=cell(DIMENSION,1);
    for k=1:DIMENSION
        leg{k,1}=num2str(k);
    end
    title('Lyapunov Spectrum');xlabel('Longitud(muestras de tiempo)');ylabel('Lyapunov Exponents');legend(leg);
    saveas(gcf,strcat(FOLDER,'Lyapunov_Spectrum'),'jpg'); hold off
end

%DELETE CREATED FILES
delete('lyap_spec.dat','delvec.dat','Precipitation.dat')
end