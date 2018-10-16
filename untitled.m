clear all;
close all;
X = dlmread('xCellcenter.txt');
Y = dlmread('yCellcenter.txt');
U = dlmread('uxVelocity.txt');
V = dlmread('uyVelocity.txt');

contourf(X,Y,V)

figure
plot(X(16,:),V(16,:))

figure
plot(Y(:,16),U(:,16))