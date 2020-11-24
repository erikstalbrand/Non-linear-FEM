%% Exercise 1.1 b) Newton-Raphson method

clc;
clear
close all;

load_steps = 300;
y0 = sind(60);      % Initial angle
u = 0;              % Initial displacement
du = 0;             % Initial displacement
kbar = 3;           % Initial given value
force_curve = zeros(load_steps, 1);   % Normalised force
ubar = zeros(load_steps, 1);          % Normalised diplacement
TOL = 1E-10;
duTOL = 1E-10;
f = 0;
df = 0.02;
f_int = 0; %zeros(load_steps, 1);

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    f = f + df;
    iter = 1;
    res = f_int - f;
    
    while abs(res) > TOL && abs(du) > duTOL || iter == 1
        disp(['Residual: ', num2str(abs(res))]);

        Kt = 2*(1 + (y0^2-1)/((1-2*u*y0 + u^2)^(3/2))) + kbar;
        
        du = -res/Kt;
        u = u + du;
        
        lambda = sqrt(1 - 2*u*y0 + u^2);
        f_int = 2*(1/lambda - 1)*(y0-u) + kbar*u;
        
        % New states
        res = f_int - f;
        iter = 2;

    end
    
    force_curve(n+1) = f;     % Numerical force results from Euler Forward
    ubar(n+1) = u;            % Normalised displacement
end


Pbar = [];
p = 0;
u = 0:0.05:2
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
legend([p1 p2], 'Newon-Raphson', 'True path')
xlim([0 2])


