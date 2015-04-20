%% This file visulaized the calculated transfer rates
clear all; clc; fig=0;
% close all;

%%

dir='C:\Users\amirhossein\Google Drive\Research\Exciton\Data\Environmental Effect\Resonance-Energy-Transfer-Rate\Transfer-Ex_A1(07,05)-to-Ex_A1(08,07)-Ckappa(1.0)\';
FileName=[dir,'transitionRates12.dat'];
kappa12=load(FileName);

FileName=[dir,'transitionRates21.dat'];
kappa21=load(FileName);

FileName=[dir,'theta.dat'];
theta=load(FileName);

FileName=[dir,'c2c.dat'];
c2c=load(FileName);

%%
% fig=fig+1; figure(fig); hold on; box on;
% surf(theta,c2c,kappa12,'EdgeColor','none');
% axis tight;
% 
% fig=fig+1; figure(fig); hold on; box on;
% surf(theta,c2c,kappa21,'EdgeColor','none');
% axis tight;

%%
nTheta = numel(theta);
fig=fig+1; figure(fig); box on;
plot(theta(2:nTheta),kappa12(2:nTheta),'-','LineWidth',3); hold on;
axis tight;

fig=fig+1; figure(fig); box on;
plot(theta(2:nTheta),kappa21(2:nTheta),'-','LineWidth',3); hold on;
axis tight;

%%
% fig=fig+1; figure(fig); hold on; box on;
% plot(c2c/1e-9,kappa12,'-k','LineWidth',3);
% 
% figure(fig); hold on; box on;
% plot(c2c/1e-9,kappa21,'-r','LineWidth',3);
% axis tight;