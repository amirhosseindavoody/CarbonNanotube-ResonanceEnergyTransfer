%% This code visulizes the free electron transition energies calculated and saved by the frotran code
clear all; clc; fig=0;
close all;

%% Load data
dir='C:\Users\Amirhossein\Google Drive\Research\Exciton\Data\test-rates\Transfer-(07,05)-Ex0_Ep-iSub(1)-Length(00nm)-Center(00nm)-Ckappa(2.0)-to-(07,05)-Ex0_A2-iSub(1)-Length(00nm)-Center(00nm)-Ckappa(2.0)-C2C( 1.2nm- 1.2nm)-Theta(000-090)-partition(1)\';
FileName=[dir,'free_electron_transition_pp.dat'];
E_free_eh_pp=load(FileName);

FileName=[dir,'free_electron_transition_pm.dat'];
E_free_eh_pm=load(FileName);

FileName=[dir,'free_electron_transition_mp.dat'];
E_free_eh_mp=load(FileName);

FileName=[dir,'free_electron_transition_mm.dat'];
E_free_eh_mm=load(FileName);

nKcm=size(E_free_eh_pp,1);
Kcm_vec = linspace(1,nKcm,nKcm)-((nKcm-1)/2);
nkr=size(E_free_eh_pp,2);
kr_vec = linspace(1,nkr,nkr)-((nkr-1)/2);

%% plot E_free_eh_pp
fig=fig+1; figure(fig); box on;
plot(Kcm_vec,E_free_eh_pp,'-','LineWidth',3); hold on;
axis tight;

fig=fig+1; figure(fig); box on;
plot(kr_vec,E_free_eh_pp,'-','LineWidth',3); hold on;
axis tight;

%% plot E_free_eh_mm
fig=fig+1; figure(fig); box on;
plot(Kcm_vec,E_free_eh_mm,'-','LineWidth',3); hold on;
axis tight;

fig=fig+1; figure(fig); box on;
plot(kr_vec,E_free_eh_mm,'-','LineWidth',3); hold on;
axis tight;