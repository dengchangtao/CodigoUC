%% Para visualizar
tic
clc
clear all
close all
% CASE INFORMATION
ruta='v0/results/'
param=load([ruta 'param.dat']);
time=load([ruta 'Time99.dat']);
number=param(1);
Nbx=param(2);
Nby=param(3);
n=param(6);
dit=param(7);
kappa=1e-6;

% Nsav=floor(n/dit); 
% Htot=zeros(Nbx,Nby); Utot=zeros(Nbx,Nby); Vtot=zeros(Nbx,Nby);
% Qxtot=zeros(Nbx,Nby); Qytot=zeros(Nbx,Nby);
%perp
mov_out=1;
figure('Position',  [  32         136        .7*1226         .7*729])

graf=1;
series=0;

it0=find(time<=0,1,'last');
itf=find(time<=2000,1,'last');
if mov_out==1
    Movie=avifile('2D.avi','FPS',5,'compression','none','quality',100);
end

for i=it0*dit:dit:itf*dit
    eval(['gunzip(''' ruta 'SOL2D.' int2str(i) '.dat.gz'')'])   % unzip file
    eval(['load ' ruta 'SOL2D.' int2str(i) '.dat'])             % load file
    system(['rm ' ruta 'SOL2D.' int2str(i) '.dat']);            % remove unzipped file
    S=reshape(SOL2D,Nbx,Nby,6);
    X=S(:,:,1); Y=S(:,:,2); Zf=S(:,:,3);
    H=S(:,:,4); U=S(:,:,5); V=S(:,:,6);
    H(H<=kappa)=nan;
    clear S
    Zf2=Zf;
    %Zf2(Zf<=-2)=nan;
    %Zf2(Zf>=2)=nan;
    switch graf
        case 1
            if i==it0*dit
                surf(X,Y,Zf,'edgecolor','none','FaceColor',[0.8 0.5 0]);hold on;
                p3=surf(X,Y,H+Zf,'edgecolor','none','facecolor',[0 1 1]);
                %p3=surf(X,Y,H+Zf,'edgealpha',0.31,'facecolor',[0 1 1]);
                t=title([num2str(time(i/dit+1)) '[s]' ],'fontsize',20);
                daspect([1 1 0.05])
                h=camlight;
                lighting phong
                lightangle(h,-125,45)   
                %lightangle(h,20,45)
                %view(125,45)
                %view(20,15)
                view(2)
                
                axis tight
                xlabel('x [m]','fontsize',14);
                ylabel('y [m]','fontsize',14);
                zlabel('z [m]','fontsize',14);
                set(gca,'fontsize',14);
                %colorbar
            else
                set(t,'string',[num2str(time(i/dit+1)) '[s]']);
                set(p3,'zdata',H+Zf);
                drawnow;
%                 for hor=1:5:100
%                     for vert=30:5:60 
%                 lightangle(h,-90+360*hor/100,90*vert/100)
%                 disp([ 'hor=' num2str(-90+360*hor/100)]);
%                 disp(['ver=' num2str(90*vert/100)]);
%                 drawnow
%                 
%                     end
%                 end
            end
        case 2
            rangox=1:1:size(X,1);
            rangoy=1:1:size(Y,2);
            if i==it0*dit
                subplot(121)
                p2=quiver(X(rangox,rangoy),Y(rangox,rangoy),U(rangox,rangoy),V(rangox,rangoy),1);
                axis equal
                %view(-90,90);
                view(2)
                axis tight
                t=title(num2str(time(i/dit+1)));
                subplot(122)
                surf(X,Y,Zf,'edgecolor','none','FaceColor',[0.8 0.3 0]);hold on;
                p3=surf(X,Y,H+Zf,'edgecolor','none');
                camlight
                lighting phong
                % shading interp
                %caxis([-2 5]);
                %view(-20,35);
                view(2)
                daspect([1 1 0.01])
                %axis equal;
                axis tight
                colorbar
            else
                set(p2,'udata',U(rangox,rangoy));
                set(p2,'vdata',V(rangox,rangoy));
                set(t,'string',num2str(time(i/dit+1)));

                set(p3,'zdata',H+Zf);
                drawnow;
            end
    end

    switch series
        case 1
            i/n
            if i==it0*dit
                ftcurvas=fopen([ruta 'tcurvas.dat'],'wt');
                curva0=load('curva1_0.dat');
                
                fhcurva0 = fopen([ruta 'hcurva0.dat'], 'wt');
                datoscurva=zeros(size(curva0,1),1);
                for ic=1:size(curva0,1)
                     datoscurva(ic,1)=H(curva0(ic,1),curva0(ic,2));
                end
                fprintf(fhcurva0, '%10.10f \n', datoscurva);
                fprintf(ftcurvas,'%10.10f \n', time(i/dit+1));

                curva500=load('curva1_500.dat');
                fhcurva500 = fopen([ruta 'hcurva500.dat'], 'wt');
                datoscurva=zeros(size(curva500,1),1);
                for ic=1:size(curva500,1)
                    datoscurva(ic,1)=H(curva500(ic,1),curva500(ic,2));
                end
                fprintf(fhcurva500, '%10.10f \n', datoscurva);
                %
                curva1000=load('curva1_1000.dat');
                fhcurva1000 = fopen([ruta 'hcurva1000.dat'], 'wt');
                datoscurva=zeros(size(curva1000,1),1);
                for ic=1:size(curva1000,1)
                    datoscurva(ic,1)=H(curva1000(ic,1),curva1000(ic,2));
                end
                fprintf(fhcurva1000, '%10.10f \n', datoscurva);

            else
                datoscurva=zeros(size(curva0,1),1);
                for ic=1:size(curva0,1)
                    datoscurva(ic,1)=H(curva0(ic,1),curva0(ic,2));
                end
                fprintf(fhcurva0, '%10.10f \n', datoscurva);
                fprintf(ftcurvas,'%10.10f \n', time(i/dit+1));

                datoscurva=zeros(size(curva500,1),1);
                for ic=1:size(curva500,1)
                    datoscurva(ic,1)=H(curva500(ic,1),curva500(ic,2));
                end
                fprintf(fhcurva500, '%10.10f \n', datoscurva);


                %curva1000=load('curva1000.dat');
                datoscurva=zeros(size(curva1000,1),1);
                for ic=1:size(curva1000,1)
                    datoscurva(ic,1)=H(curva1000(ic,1),curva1000(ic,2));
                end
                fprintf(fhcurva1000, '%10.10f \n', datoscurva);
                %
            end

    end
    %pause(0.5)    
    if mov_out==1 %HACE GRAFICOS PARA MOVIE
        frm=getframe(gcf);
        Movie=addframe(Movie,frm);
    end
end
fclose all

if mov_out==1
Movie=close(Movie);
close;
end

close all
