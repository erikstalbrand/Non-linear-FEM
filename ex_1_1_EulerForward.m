%% Exercise 1.1 a) Euler forward method
clc;
clear;
close all;

load_steps = 40;
P = 0;                                  % Initial load
y0 = sind(60);                          % Initial angle
u = 0;                                  % Initial displacement
%du = 0;                                % Initial displacement
kbar = 1;                               % Initial given value
force_curve = zeros(load_steps, 1);     % Normalised force
ubar = zeros(load_steps, 1);            % Normalised diplacement

dP = 0.15;

for n = 1:load_steps

    P = P + dP;
    
    Kt = 2*(1 + (y0^2-1)/((1-2*u*y0 + u^2)^(3/2))) + kbar;
    
    du = dP/Kt;
    
    u = u + du;
    
    force_curve(n+1) = P;         % Numerical force results from Euler Forward
    ubar(n+1) = u;                % Normalised displacement
end

Pbar = [];
u = 0:0.05:2;
for i = 1:length(u)  % True equilibrium path
    
    lambda = sqrt(1 - 2*u(i)*y0 + u(i)^2);
    p = 2*(1/lambda - 1)*(y0-u(i)) + kbar*u(i);

    Pbar = [Pbar; p];
end

figure(1)
p1 = plot(ubar, force_curve, 'LineWidth', 2);
hold on 
p2 = plot(u, Pbar, '^', 'LineWidth', 1.5);
hold off
xlabel('Displacement')
ylabel('Force')
legend([p1 p2], 'Euler Forward', 'True path')
xlim([0 2])
