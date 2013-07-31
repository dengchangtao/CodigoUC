close all
clear all
clc

%% Parámetros
L1=30;%cm debe ser >b/2
L2=15;%cm
L3=30;%cm
b=10;%cm
W=15;%cm
dx=0.25;%cm
dy=0.25;%cm,
pendiente=1/500;%sin(theta)=tan(theta)


%% los dx

x=0;%en cm
figure
N=1;

while x<=L1+L2+L3 %lineas verticales
    y1=0;%cm
    y2=W+b;%cm    
    X1(N)=x; %cm - el definitivo que se guarda    
    %dx=0.1*dx;%Aqui paso dx  de mm a cm
    line([x,x],[y1,y2], 'color','b', 'linewidth',0.1); hold on;
    x=x+dx ;
    N=N+1;    
end
% los dy
y=0;
M=1;
while y<=b+W %lineas horizontales
   Y1(M)=y;
   x1=0;
   x2=L1+L2+L3;         
   line([x1,x2],[y,y],'color','b');hold on;
   %dy=dy*0.1;
   y=y+dy;
   M=M+1;
end
axis equal

N=N-1; M=M-1;
X=zeros(M,N);
Y=zeros(M,N);
Z=zeros(M,N);
%% de cm a m
X1=X1/100;
Y1=Y1/100;
for i=1:M
   X(i,:)=X1; 
end
for i=1:N
   Y(:,i)=Y1'; 
end
X=[X;X1];
Y=[zeros(1,N); Y+0.01/100];
Y1=[0 Y1+0.01/100];

%% Canal base sin pared abajo y arriba
figure
Z=zeros(M,N);
M=M+1;
Z(1,:)=1;
Z(M,:)=1;
Z(find(Y1>=b/100),find(X1<L1/100))=1;%El canal no incluye el borde
Z(find(Y1>=b/100),find(X1>(L1+L2)/100))=1;
for i=1:size(Z,1) % la pendiente
   Z(i,:)=Z(i,:)+(X1(end)-X1)*pendiente;    
end
mesh(X,Y,Z);


XYZ=[reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)];
%save XYZsinparedes.dat -ASCII XYZ
%clear XYZ

%axis equal
%grid off
%% con pared abajo y arriba
figure
mesh(X,Y,Z);
axis equal
XYZ=[reshape(X',[],1) reshape(Y',[],1) reshape(Z',[],1)];
save XYZconparedes.dat -ASCII XYZ
%clear XYZ


%% SEÑALES
 close all
% 
%Oeste
indcanal=find(Y1<b/100 & Y1>0);
indnocanal=find(Y1>=b/100);
indnocanal=[indnocanal find(Y1==0)];
tf=300;
dt=1;
time=0:dt:tf;
Q=255/(100*100*100);%El caudal total m3/s
h0=0.01;%La altura



qx(indcanal)=Q/(length(indcanal)*h0*(dy/100)); %El caudal unitario por celda
qx(indnocanal)=0;
for i=1:length(time)
    Qx(:,i)=qx;
end
Qx=Qx';
Qx=reshape(Qx,[],1);


qy(indcanal)=0;
qy(indnocanal)=0;
for i=1:length(time)
    Qy(:,i)=qy;
end
Qy=Qy';
Qy=reshape(Qy,[],1)

eta(indcanal)=0;
eta(indnocanal)=0;
for i=1:length(time)
    Eta(:,i)=eta;
end
Eta=Eta';
Eta=reshape(Eta,[],1)

%ESTE
% 
% h(indcanal)=1/100;
% h(indnocanal)=0;
% 
% for i=1:length(time)    
%     if time(i)<=4
%        H(:,i)=1.2*h';    
%     else
%        H(:,i)=h';  
%     end
% end

hs2(1:length(time))=0.01;
%varia desde 0.012 hasta 0.01;
% for i=1:length(time)
%     if time(i)<=4        
%        hs2(i)=(0.012-0.01)/(0-4)*time(i)+0.012; 
%     else
%        hs2(i)=0.01;        
%     end
% end
%Parámetros de adimensionalizaciòn
L=1;
U=1;
H=1; 
T=L/U;


time=time/T;
hs2=hs2/H;

hs2=[time' hs2'];

save  hs2.dat   -ASCII hs2
save  Suh_xi0.dat -ASCII Qx
save  Svh_xi0.dat -ASCII Qy
save  Seta_xi0.dat -ASCII Eta
disp([ 'Nbx=' num2str(N)]);
disp(['Nby=' num2str(M)]);


