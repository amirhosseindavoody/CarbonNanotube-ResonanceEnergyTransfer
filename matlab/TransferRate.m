%% This file visulaized the calculated transfer rates
clear all; clc; fig=0;
% close all;

%%

dir='C:\Users\Amirhossein\Google Drive\Research\Exciton\Data\transition_rates_temperature_sweep\Transfer-(07,05)-iSub(1)-Length(00nm)-Center(00nm)-Ckappa(2.0)-to-(08,06)-iSub(1)-Length(00nm)-Center(00nm)-Ckappa(2.0)-C2C( 1.2nm)-Temperature(010K-500K)\';
FileName=[dir,'transition_rates.dat'];
raw_data=load(FileName);

raw_data = ctranspose(raw_data);

temperature = raw_data(1,:);
kappa_12_perp = raw_data(2,:);
kappa_21_perp = raw_data(3,:);
kappa_12_par = raw_data(4,:);
kappa_21_par = raw_data(5,:);

%%
fig=fig+1; figure(fig); box on;
plot(temperature,kappa_12_perp,'-*','LineWidth',3); hold on;
plot(temperature,kappa_21_perp,'-*','LineWidth',3); hold on;

fig=fig+1; figure(fig); box on;
plot(temperature,kappa_12_par,'-*','LineWidth',3); hold on;
plot(temperature,kappa_21_par,'-*','LineWidth',3); hold on;
axis tight;