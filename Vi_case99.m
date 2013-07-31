%%% Visualization of the SV Model Results
% Shallow water equations solver
% Curvilinear coordinate system
% Approximate Riemann Solver VFRoe-ncv plus 2nd order Hydrostatic Reconstruction

% Maricarmen Guerra Paris
tic
clc
% close all
clear all

% CASE INFORMATION
cd v0/results/
param=load('param.dat');
% Case Number:
% number=99;
number=param(1);
% Name: Hydraulic Jump over bump
n=param(6);          % Number of result files
dit=param(7);          % Results every dit iteration
% Grid Info
% Rectangular and constant grid
Nbx=param(2);
Nby=param(3);
savefile='runxx_dx025_dy025_600s';

Nsav=floor(n/dit);
Htot=zeros(Nbx,Nby); Utot=zeros(Nbx,Nby); Vtot=zeros(Nbx,Nby);
Qxtot=zeros(Nbx,Nby); Qytot=zeros(Nbx,Nby);

mov_out=1;
figure('Position',  [  32         136        .7*1226         .7*729])

time=load('Time99.dat');
%
% if mov_out==1
% Movie=avifile('LHF_tmp','FPS',12,'compression','None','quality',100);
% end

Vol0=zeros(fix(n/dit)+1,1);

% j=1;
k=1;

for i=0:dit:n
    %  for i=9000:dit:n
    eval(['gunzip(''SOL2D.' int2str(i) '.dat.gz'')'])   % unzip file
    eval(['load SOL2D.' int2str(i) '.dat'])             % load file
    %     eval(['SOL2D=SOL2D_' int2str(i) ';'])
    %     eval(['clear SOL2D_' int2str(i) ';'])               % clear variable

    system(['rm SOL2D.' int2str(i) '.dat']);            % remove unzipped file
    %     S=SOL2D;
    S=reshape(SOL2D,Nbx,Nby,6);
    %     S=reshape(SOL2D,1+(Nbx-1)/4,1+(Nby-1)/4,6);
    X=S(:,:,1); Y=S(:,:,2); Zf=S(:,:,3);
    H=S(:,:,4); U=S(:,:,5); V=S(:,:,6);


    Theta=1/3;
    dx=X(2,1)-X(1,1);dy=Y(1,2)-Y(1,1);
    Gamma=Zf*0;
    b1=[ones(1,64) 1:-1/63:0]'; g1=[0:1/63:1 1:-1/63:0]'; r1=[ 0:1/63:1 ones(1,64)]';
    map=[r1 g1 b1];
    % %     % Transforming results into matrices
    Vol0(k)=sum(sum(H))*dx*dy;
    kappa=10e-5;
    iind=find(H<=kappa);
    H(iind)=0;
    U(iind)=0;
    V(iind)=0;

    Hzf=H+Zf;

    Htot=Htot+H/Nsav;
    Utot=Utot+U/Nsav;
    Vtot=Vtot+V/Nsav;
    Qxtot=Qxtot+H.*U/Nsav;
    Qytot=Qytot+H.*V/Nsav;

    if mov_out==1

        % subplot(1,2,1),
        % % % plot(X,Zf(:,16),'.-k'), hold on, plot(X(:,16),Zf(:,16)+H(:,16),'.-b')%, axis([19 26 0.3 1])
        % plot(Y,Zf(1,:),'.-k'), hold on, plot(Y(1,:),Zf(1,:)+H(1,:),'.-b')
        % plotyy(Y(1,:),Zf(1,:),Y(1,:),Zf(1,:)+H(1,:))

        % subplot(1,2,2),plot(X(:,16),U(:,16),'.-r')%, axis([19 26 -1 1])


        % % % Gamma(2:end-1,2:end-1)=1/dx*...
        % % %     (Theta*(V(3:end,3:end)-V(1:end-2,3:end))+...
        % % %     (1-2*Theta)*(V(3:end,2:end-1)-V(1:end-2,2:end-1))+...
        % % %     Theta*(V(3:end,1:end-2)-V(1:end-2,1:end-2)))-...
        % % %     1/dy*...
        % % %     (Theta*(U(3:end,3:end)-U(3:end,1:end-2))+...
        % % %     (1-2*Theta)*(U(2:end-1,3:end)-U(2:end-1,1:end-2))+...
        % % %     Theta*(U(1:end-2,3:end)-U(1:end-2,1:end-2)));
        % % %
        % % % colormap(map) , pcolor(Y,flipud(X),Gamma), caxis([-0.05 0.05]), shading interp, hold on%,'facealpha',0.5)% to get blue and red for vorticity
        % % % contour(Y,flipud(X),Zf,[0:0.05:2],'k','LineWidth',2),axis([0 30 5 30]), hold off
        % % %
        %     quiver(X,Y,U,V,1), axis([20 25 5 10])% , axis([0 30 0 30])
        %subplot(211)


        surf(X,Y,Hzf,'edgecolor','none','cdata',Hzf)%,'facealpha',0.1)
        %caxis([9 11]);
        hold on
        colormap jet
        colorbar
        caxis([-1 1])
        zlim([-1 1])
        surf(X,Y,Zf,'FaceColor',[0.8 0.5 0],'edgealpha',0.1),   hold off
        axis equal
        %subplot(212)
        %    plot(Hzf(:,1))
        %    ylim([-2 2])

        %       shading interp
        %        axis([4 23 0 30 0 1])
        %        axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y)) min(min(Zf)) max(max(Hzf))])
        %view(-70,50)   %view(-50,80)    %view(360,0)    %view(35,60)

        % % % % cd ../out_matlab
        % % % %     caxis([0 0.21])
        % % % %     %colormap(flipud(gray))
        % % % %     colorbar
        %axis equal

        t=time(k);
        k=k+1;

        title(['t = ',num2str(t,'%6.0f'),'[s]'],'FontSize',14)

        xlabel('X(m)'), ylabel('Y(m)'), zlabel('Z(m)')

        frm=getframe(gcf);
        %      frm=getframe;
        %Movie=addframe(Movie,frm);
        clf;

        % contour(X,Y,Zf,[0:0.05:2],'LineWidth',2), hold on, quiver(X,Y,Utot,Vtot,2), axis([5 25 0 30])
    else
        k=k+1;
    end
    if mod(i,10)==0
        clc, i, toc
    end
end

% if mov_out==1
% %     cd ..
% Movie=close(Movie);
% end
close;

% %% Vorticity
%
% Gamma(2:end-1,2:end-1)=1/dx*...
%     (Theta*(Vtot(3:end,3:end)-Vtot(1:end-2,3:end))+...
%     (1-2*Theta)*(Vtot(3:end,2:end-1)-Vtot(1:end-2,2:end-1))+...
%     Theta*(Vtot(3:end,1:end-2)-Vtot(1:end-2,1:end-2)))-...
%     1/dy*...
%     (Theta*(Utot(3:end,3:end)-Utot(3:end,1:end-2))+...
%     (1-2*Theta)*(Utot(2:end-1,3:end)-Utot(2:end-1,1:end-2))+...
%     Theta*(Utot(1:end-2,3:end)-Utot(1:end-2,1:end-2)));
% figure,
% colormap(map) , pcolor(X,Y,Gamma), caxis([-0.15 0.15]),shading interp, hold on%,'facealpha',0.5)% to get blue and red for vorticity
% colorbar, view(90,90)
% contour(X,Y,Zf,[0:0.05:2],'k','LineWidth',1),axis([5 25 0 30]), hold off
%          axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y)) min(min(Zf)) max(max(Hzf))])
%
% % contour(X,Y,Zf,[0:0.1:2],'LineWidth',2), hold on, quiver(X,Y,Utot,Vtot,5), axis([0 30 0 30])
% %surf(X,Y,Gamma,'facealpha',0.5)
%
% % % % colormap(map) , surf(X,Y,Gamma)%,'facealpha',0.5)% to get blue and red for vortivity
% % figure,
%
% % save tiempo.txt -ASCII time
% Utot0=Utot;Vtot0=Vtot;aa=find(gt(X,21));Utot(aa)=0;Vtot(aa)=0;
%
% figure, contour(X,Y,Zf,[0:0.05:2],'LineWidth',2), hold on
%     quiver(X,Y,Utot,Vtot,1), axis([5 25 0 30]), hold off, view(90,90)
%
% dx1=0.5;dy1=0.5;
% X1=X(1):dx1:X(end);X1=X(end):-dx1:X(1);Y1=Y(1):dy1:Y(end);[X2 Y2]=meshgrid(X1,Y1);
%
% Utot1=interp2(X',Y',Utot',X2,Y2,'cubic');    Vtot1=interp2(X',Y',Vtot',X2,Y2,'cubic');
%
% figure, contour(X,Y,Zf,[0:0.05:2],'LineWidth',2), hold on
% quiver(X2,Y2,10*Utot1,10*Vtot1,0,'b'), axis([5 25 0 30]), % *10 pour avoir vitesse en dm/s
% quiver(6.5,1, 0,1,0,'LineWidth',1);
% text(6, 1 , '0.1 m/s');
% hold off, view(90,90)
%
%
% % hold on
% % plot(10.28,15,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
% % plot(13.09*ones(18,1),[2.23 3.85 4.7 5.55 6.4 7.25 8.1 8.95 9.8 10.65 11.5 12.35 13.2 14.9 16.6 18.3 20.0 21.7],'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
% % plot(14.71*ones(18,1),[2.23 3.85 4.7 5.55 6.4 7.25 8.1 8.95 9.8 10.65 11.5 12.35 13.2 14.9 16.6 18.3 20.0 21.7],'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
%
% %     cd ../out_matlab
% %     if ne(exist([savefile '.mat']),2)
% %      save(savefile,'X', 'Y', 'Zf', 'Htot', 'Utot', 'Vtot', 'Gamma')
% %     cd ..
% %     else
% %      fprintf(1,'Warning, File already exist ! \n')
% %      cd ..
% %     end
%
% %     if mov_out==1
% %         cd results
% %         system('ffmpeg -y -i LHF_tmp.avi -sameq LHF01.avi');
% %         system('rm -v LHF_tmp.avi');
% %         cd ..
% %     end
%
% % fname = 'a.dat';
% % [fid, msg] = fopen(fname, 'r', 'ieee-be'); % Use 'ieee-le' for
% % Little-Endian
