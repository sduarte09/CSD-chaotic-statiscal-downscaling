function [E1, E2] = cao_deneme_mod(x,tao,mmax)
%x : time series
%tao : time delay
%mmax : maximum embedding dimension
%reference:Cao, L. (1997), ``Practical method for determining the minimum
%embedding dimension of a scalar time series'', Physcai D, 110, pp. 43-50. 
%original author:"Merve Kizilkaya"
% Modified by Santiago Duarte and Nicolas Lopez
N=length(x);
a=NaN(length(x)-tao,1);
d=NaN(length(x)-tao,1);
E=NaN(length(mmax),1);
Ey=NaN(length(mmax),1);
for m=1:mmax
    M=N-m*tao;
    Y=psr_deneme(x,m,tao,M);
    for n=1:M
        distance=max(abs(bsxfun(@minus,Y,Y(n,:))),[],2);       
        distance2=distance;
        FI=find(distance==min(distance)); 
        if length(FI)==1
            distance2(FI)=max(distance);
            FG=find(distance2==min(distance2)); 
            FI(2)=FG(1);
        end
        if n==4647
            a=1;
        end
        nearpos=FI(2);
        neardis=abs(distance(nearpos)-distance(n));
        
        newpoint=[Y(n,:) x(n+m*tao)];
        newneig=[Y(nearpos,:) x(nearpos+m*tao)];
        R1=max(abs(newpoint-newneig),[],2);
        
        a(n)=R1/neardis;
        d(n)=abs(x(n+m*tao)-x(nearpos+m*tao));
        
    end
    E(m)=sum(a)/M;
    Ey(m)=sum(d)/M;
end
E1=E(2:end)./E(1:end-1);
E2=Ey(2:end)./Ey(1:end-1);

function Y=psr_deneme(x,m,tao,npoint)
%Phase space reconstruction
%x : time series
%m : embedding dimension
%tao : time delay
%npoint : total number of reconstructed vectors
%Y : M x m matrix
% original author:"Merve Kizilkaya"
% Modified by Santiago Duarte and Nicolas Lopez
N=length(x);
if nargin == 4
    M=npoint;
else
    M=N-(m-1)*tao;
end

Y2=zeros(m,M);

for i=1:m
    Y2(i,:)=x((1:M)+(i-1)*tao);
end
Y=Y2';

