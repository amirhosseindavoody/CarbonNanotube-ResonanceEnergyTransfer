%% Integrals
clear all; clc; fig=0;
% close all;

%%
theta = 0;
sinTheta = sin(theta);
cosTheta = cos(theta);

x = 0;
xp = 1;

c2c = 1.2;
radius = 0.4;
radiusp = 0.4;

nPhi = 1000;
phi = linspace(0,2*pi,nPhi);
sinPhi = sin(phi);
cosPhi = cos(phi);
dPhi = phi(2)-phi(1);

iPhip = 1;

nmu = 1000;
mu_max = 73;
mu_min = -mu_max;
mu = linspace(mu_min,mu_max,nmu);
F = zeros(nmu,1);

for imu = 1:nmu
    imu
    arg1 = (xp*cosTheta -  radiusp*cosPhi(iPhip)*sinTheta - x).^2 + (xp*sinTheta + radiusp*cosPhi(iPhip)*cosTheta - radius*cosPhi).^2 + (c2c + radiusp*sinPhi(iPhip) - radius*sinPhi).^2;
    arg2 = exp(1i*mu(imu)*phi) .* sqrt(arg1);
    F(imu,1) = sum(arg2)*dPhi;
end;
    

fig=fig+1; figure(fig); box on; hold on;
% plot(mu,F,'LineWidth',3);
plot(mu,abs(F),'LineWidth',3);