%% This file visualizes the results of the fortran program for CNT bethe salpeter equation
clear all; clc; fig=10;
close all;
dir='C:\Users\amirhossein\Google Drive\Research\Exciton\test Data\ForsterFiles\';
eV=1.6e-19;

%% plot CNT unit cell
% FileName=[dir,'posA.dat'];
% posA=load(FileName);
% FileName=[dir,'posB.dat'];
% posB=load(FileName);
% 
% fig=fig+1; figure(fig); hold on; box on;
% plot(posA(:,1),posA(:,2),'b.','MarkerSize',20);
% plot(posB(:,1),posB(:,2),'r.','MarkerSize',20);
% 
% plot(posA(1,1),posA(1,2),'g.','MarkerSize',20);
% plot(posB(1,1),posB(1,2),'k.','MarkerSize',20);
% 
% axis equal; axis tight;
% return;

%% plot exciton energy dispersion of the first cnt
FileName=[dir,'cnt1_kvec.dat'];
k_vec1=load(FileName);
nKcm1 = numel(k_vec1);
FileName=[dir,'cnt1_Ex0_A2.dat'];
cnt1_Ex0_A2=load(FileName);

fig=fig+1; figure(fig); hold on; box on;
plot(k_vec1/1e9,cnt1_Ex0_A2/eV,'-b','LineWidth',3);
% return;
%% plot exciton energy dispersion of the second cnt
FileName=[dir,'cnt2_kvec.dat'];
k_vec2=load(FileName);
nKcm2 = numel(k_vec2);
FileName=[dir,'cnt2_Ex0_A2.dat'];
cnt2_Ex0_A2=load(FileName);

figure(fig); hold on; box on;
plot(k_vec2/1e9,cnt2_Ex0_A2/eV,'-r','LineWidth',3);

% axis tight;
% return;

%% plot crossing points
FileName=[dir,'crossingPoints.dat'];
crossingPoints=load(FileName);

nC = size(crossingPoints,1);
tmp = (nKcm1-1)/2+1;

figure(fig); hold on; box on;
for i = 1:nC
    plot(k_vec1(tmp+crossingPoints(i,3))/1e9,cnt1_Ex0_A2(tmp+crossingPoints(i,3),crossingPoints(i,1))/eV,'*-k','LineWidth',3);
end

axis tight;
% return;

%% plot points with equal energies
% FileName=[dir,'sameEnergy.dat'];
% sameEnergy=load(FileName);
% 
% nS = size(sameEnergy,1);
% tmp = (nKcm1-1)/2+1;
% 
% fig = fig+1; figure(fig); hold on; box on;
% for i = 1:nS
%     plot(k_vec1,cnt1_Ex0_A2/eV,'-b','LineWidth',3); hold on;
%     plot(k_vec2,cnt2_Ex0_A2/eV,'-r','LineWidth',4);
%     plot(k_vec1(tmp+sameEnergy(i,3)),cnt1_Ex0_A2(tmp+sameEnergy(i,3),sameEnergy(i,1))/eV,'*-g','LineWidth',3);
%     plot(k_vec2(tmp+sameEnergy(i,4)),cnt2_Ex0_A2(tmp+sameEnergy(i,4),sameEnergy(i,2))/eV,'*-k','LineWidth',4);
%     hold off;
%     pause;
% end
% 
% axis tight;
% return;

%% plot exciton density of states of the first cnt
FileName=[dir,'cnt1_kvec.dat'];
k_vec1=load(FileName);
nKcm1 = numel(k_vec1);
FileName=[dir,'cnt1_DOS.dat'];
cnt1_DOS=load(FileName);

fig=fig+1; figure(fig); hold on; box on;
plot(k_vec1,cnt1_DOS*eV,'-','LineWidth',3);
for i = 1:nC
    plot(k_vec1(tmp+crossingPoints(i,3)),cnt1_DOS(tmp+crossingPoints(i,3),crossingPoints(i,1))*eV,'*-k','LineWidth',4);
end
axis tight;
% for iX = 1:size(cnt1_DOS,2)
%     plot(k_vec1,cnt1_DOS(:,iX),'-','LineWidth',3);
%     pause;
% end;

% return;

%% plot exciton density of states of the first cnt
FileName=[dir,'cnt2_kvec.dat'];
k_vec2=load(FileName);
nKcm2 = numel(k_vec2);
FileName=[dir,'cnt2_DOS.dat'];
cnt2_DOS=load(FileName);

fig=fig+1; figure(fig); hold on; box on;
plot(k_vec2,cnt2_DOS*eV,'-','LineWidth',3);
for i = 1:nC
    plot(k_vec2(tmp+crossingPoints(i,3)),cnt2_DOS(tmp+crossingPoints(i,3),crossingPoints(i,2))*eV,'*-k','LineWidth',4);
end
axis tight;
% for iX = 1:size(cnt2_DOS,2)
%     plot(k_vec2,cnt2_DOS(:,iX),'-','LineWidth',3);
%     pause;
% end;

return;

%% plot exciton energy Ex0_A2
FileName=[dir,'Ex0_A2.dat'];
Ex0_A2=load(FileName);
[nKcm,nX]=size(Ex0_A2);
Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);

fig=fig+1; figure(fig); hold on; box on;
for i=1:nX
    plot(Kcm_vec,Ex0_A2(:,i)/eV,'*-r','LineWidth',3);
end;
axis tight;

%% plot exciton energy Ex_A1
FileName=[dir,'Ex1_A2.dat'];
Ex1_A2=load(FileName);
[nKcm,nX]=size(Ex1_A2);
Kcm_vec=dk*(-(nKcm-1)/2:+(nKcm-1)/2);

fig=fig+1; figure(fig); hold on; box on;
for i=1:nX
    plot(Kcm_vec,Ex1_A2(:,i)/eV,'-r','LineWidth',3);
end;
axis tight;

return;
%% plot exciton wavefunction in k-space
FileName=[dir,'Psi_A1.dat'];
tmp=load(FileName);
[nKcm,ncol]=size(tmp);
nk=ncol/nX/2;
Psi_A1=zeros(nKcm,nX,nk);

for i=1:nX
    for j=1:nk
        Psi_A1(:,i,j)=tmp(:,(i-1)*2*nk+2*j-1)+1i*tmp(:,(i-1)*2*nk+2*j);
    end;
end;
clear tmp;

FileName=[dir,'Psi0_A2.dat'];
tmp=load(FileName);
[nKcm,ncol]=size(tmp);
nk=ncol/nX/2;
Psi0_A2=zeros(nKcm,nX,nk);

for i=1:nX
    for j=1:nk
        Psi0_A2(:,i,j)=tmp(:,(i-1)*2*nk+2*j-1)+1i*tmp(:,(i-1)*2*nk+2*j);
    end;
end;
clear tmp;

FileName=[dir,'Psi1_A2.dat'];
tmp=load(FileName);
[nKcm,ncol]=size(tmp);
nk=ncol/nX/2;
kr_vec=dk*(-(nk-1)/2:+(nk-1)/2);
Psi1_A2=zeros(nKcm,nX,nk);

for i=1:nX
    for j=1:nk
        Psi1_A2(:,i,j)=tmp(:,(i-1)*2*nk+2*j-1)+1i*tmp(:,(i-1)*2*nk+2*j);
    end;
end;
clear tmp;

iKcm=(nKcm+1)/2;
fig=fig+1; figure(fig); box on; hold on;
for iX=1
    tmp(1,:)=Psi_A1(iKcm,iX,:);
    plot(kr_vec,abs(tmp),'-b','LineWidth',3);
    tmp(1,:)=Psi0_A2(iKcm,iX,:);
    plot(kr_vec,abs(tmp),'-r','LineWidth',3);
    tmp(1,:)=Psi1_A2(iKcm,iX,:);
    plot(kr_vec,abs(tmp),'-k','LineWidth',3);
end;
axis tight;
clear tmp;