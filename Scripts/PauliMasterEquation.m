%% This program solves the classical Pauli master equation in frequency and time domain
clc;
clear all;
close all;
fig = 0;
dir='C:\Users\amirhossein\Google Drive\Research\Exciton\test Data\ForsterRate (08,06) to (08,06)\';

%% Upload the transfer rates table
% FileName=[dir,'transitionRates12.dat'];
% kappa12=load(FileName);
% 
% FileName=[dir,'transitionRates21.dat'];
% kappa21=load(FileName);
% 
% FileName=[dir,'theta.dat'];
% theta=load(FileName);
% 
% FileName=[dir,'c2c.dat'];
% c2c=load(FileName);
% 
% nTheta = numel(theta);
% dTheta = (theta(nTheta)-theta(1))/(nTheta-1);

%% Calculate population percentage

nState = 4; %number of available states
P0 = zeros(nState,1);
P0(1,1)=1;

angle = zeros(nState,1);
for i = 1 : nState
   angle(i) = 0+(i-1)*(3*pi/4/(nState-1));
end

tmp = 0;
for i = 1:nState
   P0(i,1) = cos(angle(i))^2; 
   tmp = tmp + cos(angle(i))^2;
end
P0 = P0/tmp;


kappa = zeros(nState,nState);
FileName=[dir,'kappaMatrix.dat'];
kappa=load(FileName);

tranMatrix = zeros(nState,nState);
for i = 1 : nState
    for j = 1 : nState
        if (i ~= j )
            tranMatrix(i,i) = tranMatrix(i,i) + kappa (i,j);
            tranMatrix(i,j) = -kappa(j,i);
        end;
    end;
end;

[Ea,omega] = eig(tranMatrix);
C = Ea\P0;

Tmax = 100/sum(diag(omega));
nT = 5000;

t = linspace(0,Tmax,nT);
P = zeros(nState,nT);


for i = 1 : nState
    for j = 1 : nState
        P(i,:) = P(i,:) + Ea(i,j) * C(j) * exp(-omega(j,j)*t);
    end;
end;

fig = fig+1; figure(fig); hold on;
plot(t/1e-12,real(P),'LineWidth',3);
ylim([0,1]);
xlim([0,4]);
% axis tight;
box on;

%% Calculate and plot anisotropy
num = zeros(1,nT);
denom = zeros(1,nT);

for i = 1 : nState
    num = num + P(i,:)*(cos(angle(i))^2-sin(angle(i))^2);
    denom = denom + P(i,:);
end;

anisotropy = num./denom;

fig = fig+1; figure(fig); hold on;
plot(t/1e-12,real(anisotropy),'LineWidth',3);
ylim([0,1]);
xlim([0,4]);
% axis tight;
box on;