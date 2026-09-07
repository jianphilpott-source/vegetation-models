clear;
clc;
clear all;

%Parameter Input, Mortality and Rainfall Range
m = 1;
A = linspace(0,5,1000);

%initialise equilibria vectors
b0 = zeros(size(A));
b_up = zeros(size(A));
b_down = zeros(size(A));

idx = A >= 2*m;

b_up(idx) = (A(idx)+sqrt(A(idx).^2 - 4*m^2)) ./ (2*m);
b_down(idx) = (A(idx)-sqrt(A(idx).^2 - 4*m^2)) ./ (2*m);

% Calculate stability conditions

trace_up = m - 1 - b_up(idx).^2;
det_up   = m .* (b_up(idx).^2 - 1);

trace_down = m - 1 - b_down(idx).^2;
det_down   = m .* (b_down(idx).^2 - 1);

stable_up = trace_up < 0 & det_up > 0;
stable_down = trace_down < 0 & det_down > 0;
%Plotting Equilibria
figure;
hold on;

%BareSoil is stable Everywhere
plot(A, b0, 'b-', 'LineWidth', 2);

% Upper vegetated equilibrium

A_veg = A(idx);
b_veg_up = b_up(idx);

plot(A_veg(stable_up), b_veg_up(stable_up), ...
    'g-', 'LineWidth', 2);
plot(A_veg(~stable_up), b_veg_up(~stable_up), ...
    'g--', 'LineWidth', 2);

% Lower vegetated equilibrium

b_veg_down = b_down(idx);
plot(A_veg(stable_down), b_veg_down(stable_down), ...
    'g-', 'LineWidth', 2);
plot(A_veg(~stable_down), b_veg_down(~stable_down), ...
    'g--', 'LineWidth', 2);

xlabel('Rainfall, $A$', 'Interpreter', 'latex');
ylabel('Equilibrium Biomass, $b_*$', 'Interpreter', 'latex');

%Saddle-Node Bifurcation
plot(2*m, 1, 'ko', 'MarkerFaceColor', 'r');

legend('Bare soil (stable)', ...
'Vegetated (stable)', ...
'Vegetated (unstable)', ...
'Saddle-node', ...
'Location','best', 'Interpreter', 'latex');

title('Bifurcation Diagram of Homogeneous Klausmeier Model', 'Interpreter', 'latex');

grid on;
box on;

hold off;

