%% This code visulizes the free electron transition energies calculated and saved by the frotran code
clear all; clc; fig=0;
close all;

%% Load data
dir='C:\Users\Amirhossein\Google Drive\Research\Exciton\Data\test-rates\Transfer-(07,05)-Ex0_A2-iSub(1)-Length(00nm)-Center(00nm)-Ckappa(2.0)-to-(08,07)-Ex0_A2-iSub(1)-Length(00nm)-Center(00nm)-Ckappa(2.0)-C2C( 1.2nm- 1.2nm)-Theta(000-000)-partition(1)\';
FileName=[dir,'free_electron_transition_pp.dat'];
E_free_eh_pp=load(FileName);

% FileName=[dir,'free_electron_transition_pm.dat'];
% E_free_eh_pm=load(FileName);
% 
% FileName=[dir,'free_electron_transition_mp.dat'];
% E_free_eh_mp=load(FileName);
% 
% FileName=[dir,'free_electron_transition_mm.dat'];
% E_free_eh_mm=load(FileName);

nkr=size(E_free_eh_pp,2);
kr_vec = linspace(1,nkr,nkr)-((nkr-1)/2);

FileName=[dir,'cnt1_kvec.dat'];
Kcm_vec=load(FileName);

FileName=[dir,'cnt1_Ex_t.dat'];
Ex_t=load(FileName);

FileName=[dir,'crossingPoints.dat'];
crossingPoints=load(FileName);

nCrossing = size(crossingPoints,1);

%% plot E_free_eh_pp
fig=fig+1; figure(fig); box on;
plot(Kcm_vec,E_free_eh_pp,'-','LineWidth',3); hold on;
axis tight;

plot(Kcm_vec,Ex_t,'-r','LineWidth',3); hold on;
axis tight;

for iC=1:nCrossing
    plot(Kcm_vec(crossingPoints(iC,3)),Ex_t(crossingPoints(iC,3),crossingPoints(iC,1)),'*k','LineWidth',3); hold on;
end;

% fig=fig+1; figure(fig); box on;
% plot(kr_vec,E_free_eh_pp,'-','LineWidth',3); hold on;
% axis tight;

return;

%% plot E_free_eh_mm
fig=fig+1; figure(fig); box on;
plot(Kcm_vec,E_free_eh_mm,'-','LineWidth',3); hold on;
axis tight;

fig=fig+1; figure(fig); box on;
plot(kr_vec,E_free_eh_mm,'-','LineWidth',3); hold on;
axis tight;