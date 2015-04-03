%% This code plots the matrix element dependence on angle, distance, and k-vector mismatch
clear all; clc;
% close all;
fig = 0;

%% Dependence on angle

% nTheta = 2000;
% thetaMin = 0;
% thetaMax = pi;
% theta = linspace(thetaMin,thetaMax,nTheta);
% 
% K1 = 10;
% K2 = 11;
% 
% radius1 = 0.1;
% radius2 = 0.1;
% Distance = 0.3;
% 
% nPhi = 30;
% phiMin = 0;
% phiMax = 2*pi;
% phi = linspace(phiMin,phiMax,nPhi);
% dPhi = phi(2)-phi(1);
% 
% Jphi=zeros(1,nTheta);
% 
% for iPhi1 = 1:nPhi
%     for iPhi2 = 1:nPhi
%         tmp = sqrt(K2^2+K1^2-2*K1*K2*cos(theta));
%         Jphi = Jphi + dPhi^2 * exp(2*1i*(K1*(radius2*cos(phi(iPhi2))-radius1*cos(phi(iPhi1))*cos(theta))+ ...
%                                     K2*(radius1*cos(phi(iPhi1))-radius2*cos(phi(iPhi2))*cos(theta)))./sin(theta)) ...
%                              .* exp(-2*(Distance+radius2*sin(phi(iPhi2))-radius1*sin(phi(iPhi1)))./sin(theta).*tmp)./tmp;
%         
%     end;
% end;
% 
% fig = fig+1; figure(fig); hold on; box on;
% plot(theta,abs(Jphi)/max(abs(Jphi)),'-b','LineWidth',3);
% axis tight;
% 
% % fig = fig+1; figure(fig); hold on; box on;
% % plot(theta,abs(Jphi),'-k','LineWidth',3);
% % axis tight;

%% Dependence on distance

theta = 0.7;

K1 = 10;
K2 = 10.5;

radius1 = 0.1;
radius2 = 0.1;

nDistance = 200;
distanceMin=0.3;
distanceMax=1.0;
distance = linspace(distanceMin,distanceMax,nDistance);

nPhi = 30;
phiMin = 0;
phiMax = 2*pi;
phi = linspace(phiMin,phiMax,nPhi);
dPhi = phi(2)-phi(1);

Jphi=zeros(1,nDistance);

for iPhi1 = 1:nPhi
    for iPhi2 = 1:nPhi
        tmp = sqrt(K2^2+K1^2-2*K1*K2*cos(theta));
        Jphi = Jphi + dPhi^2 * exp(2*1i*(K1*(radius2*cos(phi(iPhi2))-radius1*cos(phi(iPhi1))*cos(theta))+ ...
                                    K2*(radius1*cos(phi(iPhi1))-radius2*cos(phi(iPhi2))*cos(theta)))/sin(theta)) ...
                             * exp(-2*(distance+radius2*sin(phi(iPhi2))-radius1*sin(phi(iPhi1)))/sin(theta)*tmp)/tmp;
        
    end;
end;

% fig = fig+1; figure(fig); hold on; box on;
% plot(distance,abs(Jphi)/max(abs(Jphi)),'-b','LineWidth',3);
% axis tight;

fig = fig+1; figure(fig); hold on; box on;
plot(distance,abs(Jphi),'-','LineWidth',3);
axis tight;