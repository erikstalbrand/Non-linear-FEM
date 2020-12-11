%% Exercise 3.1 St. Venant-Kirchhoff Model

clear;
close all;
clc;

v = 0.3;
E = 1;

mpara = [E, v];
flag = [1 2 3];
% defgrad = F;
tau12_plot = zeros(100, 1);
tau22_plot = zeros(100, 1);
sigma12_plot = zeros(100, 1);
sigma22_plot = zeros(100, 1);
k_plot = zeros(100, 1);
n = 1;

for k = 0:0.01:1

defgrad = [1 k 0 1]';

stress = stressMater2D1(2, mpara, defgrad);         % flag = 2 for kirchoff, tau
tau12_plot(n) = stress(3);          % 3 is (1, 2)
tau22_plot(n) = stress(2);          % 3 is (1, 2)

stress = stressMater2D1(3, mpara, defgrad);         % flag = 3 for cauchy, sigma
sigma12_plot(n) = stress(3);
sigma22_plot(n) = stress(2);

k_plot(n) = k;

n = n + 1;
end

figure(1)
p1 = plot(k_plot, tau12_plot, '-ob', 'Linewidth', 0.5);
hold on 
p2 = plot(k_plot, sigma12_plot, '-r', 'Linewidth', 1);
legend([p1 p2], '\tau_{12}', '\sigma_{12}')
xlabel('\kappa')
ylabel('\tau/E or \sigma/E')
title('St. Venant-Kirchoff Model')
grid on

figure(2)
p3 = plot(k_plot, tau22_plot, '-ob', 'Linewidth', 0.5);
hold on 
p4 = plot(k_plot, sigma22_plot, '-r', 'Linewidth', 1);
legend([p3 p4], '\tau_{22}', '\sigma_{22}')
xlabel('\kappa')
ylabel('\tau/E or \sigma/E')
title('St. Venant-Kirchoff Model')
grid on

%% Neo-Hookean Model

tau12_plot = zeros(100, 1);
tau12_plot = zeros(100, 1);
sigma12_plot = zeros(100, 1);
sigma22_plot = zeros(100, 1);
k_plot = zeros(100, 1);
n = 1;

for k = 0:0.01:1

defgrad = [1 k 0 1]';

stress = stressMater2D2(2, mpara, defgrad);         % flag = 2 for kkirchoff, tau
tau12_plot(n) = stress(3);          % 3 is (1, 2)
tau22_plot(n) = stress(2);          % 3 is (1, 2)

stress = stressMater2D2(3, mpara, defgrad);         % flag = 3 for cauchy, sigma
sigma12_plot(n) = stress(3);
sigma22_plot(n) = stress(2);

k_plot(n) = k;

n = n + 1;
end

figure(3)
p1 = plot(k_plot, tau12_plot, '-ob', 'Linewidth', 0.5);
hold on 
p2 = plot(k_plot, sigma12_plot, '-r', 'Linewidth', 1);
legend([p1 p2], '\tau_{12}', '\sigma_{12}')
xlabel('\kappa')
ylabel('\tau/E or \sigma/E')
title('Neo-Hookean Model')
grid on

figure(4)
p3 = plot(k_plot, tau22_plot, '-ob', 'Linewidth', 0.5);
hold on 
p4 = plot(k_plot, sigma22_plot, '-r', 'Linewidth', 1);
legend([p3 p4], '\tau_{22}', '\sigma_{22}')
xlabel('\kappa')
ylabel('\tau/E or \sigma/E')
title('Neo-Hookean Model')
grid on
