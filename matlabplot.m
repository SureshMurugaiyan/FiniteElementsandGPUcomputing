clear all;
close all;
format long;

X = dlmread('xCellcenter.txt');
Y = dlmread('yCellcenter.txt');
U = dlmread('uxVelocity.txt');
V = dlmread('uyVelocity.txt');
[Nx,Ny] = size(X); % Grid size
Res = dlmread('ResidualPlotting.txt');

figure (1)
colormap(parula(25))
contourf(X,Y,U,25)
axis image

figure  (2)
colormap(parula(25))
contourf(X,Y,V,25)
axis image

figure  (3)
ax1 = subplot(1,2,1);
plot(X(Nx/2,:),V(Nx/2,:))
colormap(ax1,hot(8))
axis square
ax2 = subplot(1,2,2);
 plot(Y(:,Nx/2),U(:,Nx/2))
colormap(ax2,pink)
axis square



figure  (4)
ax1 = subplot(2,2,1);
plot(Res(:,1),Res(:,2))
colormap(ax1,hot(8))
axis square
ax2 = subplot(2,2,2);
plot(Res(:,1),Res(:,3))
colormap(ax2,pink)
axis square
ax3 = subplot(2,2,3);
plot(Res(:,1),Res(:,4))
colormap(ax2,pink)
axis square

figure  (5)
surface(X,Y,0*X)
axis image
title('GRID')
xlabel('X');
ylabel('Y');
%set(get(gca,'title'),'position',[0.5 0.8 0])
hold on;
plot(X(1,:),Y(1,:),'r','Linewidth',2);
plot(X(Ny,:),Y(Ny,:),'r','Linewidth',2);
plot(X(:,1),Y(:,1),'r','Linewidth',2);
plot(X(:,Nx),Y(:,Nx),'r','Linewidth',2);


% figure
% plot(X(Nx/2,:),V(Nx/2,:))
% 
% figure
% plot(Y(:,Nx/2),U(:,Nx/2))

%  figure
%  plot(Res(:,1),Res(:,2))
%  hold on
%  plot(Res(:,1),Res(:,3))
%  plot(Res(:,1),Res(:,4))