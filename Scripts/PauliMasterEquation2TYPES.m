%% This program solves the classical Pauli master equation in frequency and time domain
clc;
clear all;
close all;
fig = 0;
dir='C:\Users\amirhossein\Google Drive\Research\Exciton\Data\ForsterRate (07,05) to (08,06)\';

%% Load experimental data
load('C:\Users\amirhossein\Google Drive\Research\Exciton\Data\Experimental Data\Anisotropy(7,5)to(8,6)');

%% Calculate population percentage

nState = 10; %number of available states


theta = zeros(nState,1);
thetaMax = pi;
thetaMin = 0;
dTheta = (thetaMax-thetaMin)/nState;
for i = 1 : nState
   theta(i) = (i-1)*dTheta;
end

tmp = 0;
P0 = zeros(2*nState,1);
for i = 1:nState
   P0(i,1) = cos(theta(i))^2; 
   tmp = tmp + cos(theta(i))^2;
end
P0 = P0/tmp;

% tmp = 0;
% P0 = zeros(2*nState,1);
% for i = 1:nState
%    P0(i,1) = 0.6; 
%    tmp = tmp + P0(i,1);
% end
% for i = nState+1:2*nState
%    P0(i,1) = 0.4; 
%    tmp = tmp + P0(i,1);
% end
% P0 = P0/tmp;

FileName=[dir,'kappaMatrix.dat'];
kappa=load(FileName);

tranMatrix = zeros(2*nState,2*nState);
for i = 1 : 2*nState
    for j = 1 : 2*nState
        if (i ~= j )
            tranMatrix(i,i) = tranMatrix(i,i) + kappa (i,j);
            tranMatrix(i,j) = -kappa(j,i);
        end;
    end;
end;

[Ea,omega] = eig(tranMatrix);
C = Ea\P0;

Tmax = 100/sum(diag(omega));
nT = 50000;

t = linspace(0,Tmax,nT);
P = zeros(2*nState,nT);


for i = 1 : 2*nState
    for j = 1 : 2*nState
        P(i,:) = P(i,:) + Ea(i,j) * C(j) * exp(-omega(j,j)*t);
    end;
end;

P1 = sum(P(1:nState,:),1);
P2 = sum(P(nState+1:2*nState,:),1);
Ptotal = sum(P(1:2*nState,:),1);
fig = fig+1; figure(fig); hold on;
plot(t/1e-12,real(P1./Ptotal),'LineWidth',3);
plot(t/1e-12,real(P2./Ptotal),'LineWidth',3);
ylim([0,1]);
xlim([0,4]);
% axis tight;
% ylim([0,1]);
box on;

%% Calculate and plot anisotropy
num = zeros(1,nT);
denom = zeros(1,nT);
for i = 1 : nState
    num = num + P(i,:)*(cos(theta(i))^2-sin(theta(i))^2);
    denom = denom + P(i,:);
end;
anisotropy1 = num./denom;

num = zeros(1,nT);
denom = zeros(1,nT);
for i = 1 : nState
    num = num + P(nState+i,:)*(cos(theta(i))^2-sin(theta(i))^2);
    denom = denom + P(nState+i,:);
end;
anisotropy2 = num./denom;

fig = fig+1; figure(fig); hold on;
plot(t(2:nT)/1e-12,abs(anisotropy1(2:nT)),'-b','LineWidth',3);
plot(t(2:nT)/1e-12,abs(anisotropy2(2:nT)),'-r','LineWidth',3);
plot(Time2,Anisotropy1,'-b','LineWidth',3);
plot(Time2,Anisotropy2,'-r','LineWidth',3);
ylim([0,0.5]);
xlim([0,4]);
% axis tight;
box on;
