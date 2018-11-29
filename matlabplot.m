clear all;
close all;
Nx = 32;  % Grid size
Ny = 32;  % Grid size
X = dlmread('xCellcenter.txt');
Y = dlmread('yCellcenter.txt');
U = dlmread('uxVelocity.txt');
V = dlmread('uyVelocity.txt');

figure
colormap(parula(25))
% contourf(V,25)
contourf(X,Y,V,25)


% figure
% plot(X(Nx/2,:),V(Nx/2,:))
% 
% figure
% plot(Y(:,Nx/2),U(:,Nx/2))


figure
% ax1 = subplot(1,2,1);
% plot(X(Nx/2,:),V(Nx/2,:))
% colormap(ax1,hot(8))
% axis square
% ax2 = subplot(1,2,2);
%  plot(Y(:,Nx/2),U(:,Nx/2))
% colormap(ax2,pink)
% axis square

  plot(V(Ny/2,:))
% plot(U(:,Nx/2))