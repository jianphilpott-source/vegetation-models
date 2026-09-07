clear all;
clc;
clear;

k = linspace(0,100,1000);

%Parameters
A = 7;
m = 2;
D_b = 0.001;
D_w = 0.1;
nu = 100;

b_star =  (A+sqrt(A^2 - 4*m^2)) / (2*m);
w_star =  (A -sqrt(A^2 - 4*m^2)) / 2;

growthRate = zeros(size(k));

for i = 1:length(k)

    M_k = [-1-b_star^2 - D_w * k(i)^2 + 1i*nu*k(i), -2*w_star*b_star;
            b_star^2, 2*w_star*b_star-m-D_b*k(i)^2];

    lambda = eig(M_k);

    growthRate(i) = max(real(lambda));
end

% Most unstable wavenumber
[maxGrowth,kIndex] = max(growthRate);
k_max = k(kIndex);

% Corresponding wavelength
wavelength = 2*pi/k_max;

figure;
plot(k,growthRate,'LineWidth',2);
hold on;

plot(k_max,maxGrowth,'ko','MarkerFaceColor','k');

yline(0,'k--');

xlabel('$k$','Interpreter','latex');
ylabel('$\mathrm{Re}(\lambda_{\max})$','Interpreter','latex');
title('Growth Rate Against Wavenumber','Interpreter','latex');

disp("k_max = " + k_max)
disp("wavelength = " + wavelength)
grid on;
box on;