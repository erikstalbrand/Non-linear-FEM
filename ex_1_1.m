%% Exercise 1.1 a) Euler forward method, k=1

clc;
clear;
close all;

%------------------------------------
kbar = 1;                               % CHANGE FOR DIFFERENT SPRING STIFFNESS
%------------------------------------

load_steps = 40;
P = 0;                                  % Initial load
y0 = sind(60);                          % Initial angle
u = 0;                                  % Initial displacement
%du = 0;                                % Initial displacement
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

if kbar == 1
    title('Euler forward, kbar = 1')
else
    title('Euler forward, kbar = 3')
end

%% %% Exercise 1.1 a) Euler forward method, k=3

clc;
clear;
%------------------------------------
kbar = 3;                               % CHANGE FOR DIFFERENT SPRING STIFFNESS
%------------------------------------

load_steps = 40;
P = 0;                                  % Initial load
y0 = sind(60);                          % Initial angle
u = 0;                                  % Initial displacement
%du = 0;                                % Initial displacement
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

figure(2)
p1 = plot(ubar, force_curve, 'LineWidth', 2);
hold on 
p2 = plot(u, Pbar, '^', 'LineWidth', 1.5);
hold off
xlabel('Displacement')
ylabel('Force')
legend([p1 p2], 'Euler Forward', 'True path')
xlim([0 2])

if kbar == 1
    title('Euler forward, kbar = 1')
else
    title('Euler forward, kbar = 3')
end
%% b) Newton-Raphson, k = 1
clc;
clear;

%------------------------------------
kbar = 1;                               % CHANGE FOR DIFFERENT SPRING STIFFNESS
%------------------------------------

load_steps = 300;
y0 = sind(60);      % Initial angle
u = 0;              % Initial displacement
du = 0;             % Initial displacement
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


figure(3)
p1 = plot(ubar, force_curve, 'LineWidth', 2);
hold on 
p2 = plot(u, Pbar, '^', 'LineWidth', 1.5);
hold off
xlabel('Displacement')
ylabel('Force')
legend([p1 p2], 'Newon-Raphson', 'True path')
xlim([0 2])
if kbar == 1
    title('Newton-Raphson, kbar = 1')
else
    title('Newton-Raphson, kbar = 3')
end

%% b) Newton-Raphson, k = 3
clc;
clear;

%------------------------------------
kbar = 3;                               % CHANGE FOR DIFFERENT SPRING STIFFNESS
%------------------------------------

load_steps = 300;
y0 = sind(60);      % Initial angle
u = 0;              % Initial displacement
du = 0;             % Initial displacement
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


figure(4)
p1 = plot(ubar, force_curve, 'LineWidth', 2);
hold on 
p2 = plot(u, Pbar, '^', 'LineWidth', 1.5);
hold off
xlabel('Displacement')
ylabel('Force')
legend([p1 p2], 'Newon-Raphson', 'True path')
xlim([0 2])
if kbar == 1
    title('Newton-Raphson, kbar = 1')
else
    title('Newton-Raphson, kbar = 3')
end
%% c) Modified Newton-Raphson method, k = 1

clc;
clear
%------------------------------------
kbar = 1;                               % CHANGE FOR DIFFERENT SPRING STIFFNESS
%------------------------------------
load_steps = 200;
y0 = sind(60);                          % Initial angle
u = 0;                                  % Initial displacement
du = 0;                                 % Initial displacement
force_curve = zeros(load_steps, 1);     % Normalised force
ubar = zeros(load_steps, 1);            % Normalised diplacement
TOL = 1E-10;
duTOL = 1E-10;
f = 0;
df = 0.05;
f_int = 0;

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    f = f + df;
    iter = 1;
    res = f_int - f;
    
    Kt = 2*(1 + (y0^2-1)/((1-2*u*y0 + u^2)^(3/2))) + kbar;  % Initial tangent stiffness
    
    while abs(res) > TOL && abs(du) > duTOL || iter == 1
        disp(['Residual: ', num2str(abs(res))]);

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
u = 0:0.05:2;

for i = 1:length(u)        % True equilibrium path
    
    lambda = sqrt(1 - 2*u(i)*y0 + u(i)^2);
    p = 2*(1/lambda - 1)*(y0-u(i)) + kbar*u(i);
    
    Pbar = [Pbar; p];

end


figure(5)
p1 = plot(ubar, force_curve, 'LineWidth', 2);
hold on 
p2 = plot(u, Pbar, '^', 'LineWidth', 1.5);
hold off
xlabel('Displacement')
ylabel('Force')
legend([p1 p2], 'Modified Newton-Raphson', 'True path')
xlim([0 2])
if kbar == 1
    title('Modified Newton-Raphson, kbar = 1')
else
    title('Modified Newton-Raphson, kbar = 3')
end

%% c) Modified Newton-Raphson method, k = 3

clc;
clear
%------------------------------------
kbar = 3;                               % CHANGE FOR DIFFERENT SPRING STIFFNESS
%------------------------------------
load_steps = 200;
y0 = sind(60);                          % Initial angle
u = 0;                                  % Initial displacement
du = 0;                                 % Initial displacement
force_curve = zeros(load_steps, 1);     % Normalised force
ubar = zeros(load_steps, 1);            % Normalised diplacement
TOL = 1E-10;
duTOL = 1E-10;
f = 0;
df = 0.05;
f_int = 0;

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    f = f + df;
    iter = 1;
    res = f_int - f;
    
    Kt = 2*(1 + (y0^2-1)/((1-2*u*y0 + u^2)^(3/2))) + kbar;  % Initial tangent stiffness
    
    while abs(res) > TOL && abs(du) > duTOL || iter == 1
        disp(['Residual: ', num2str(abs(res))]);

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
u = 0:0.05:2;

for i = 1:length(u)        % True equilibrium path
    
    lambda = sqrt(1 - 2*u(i)*y0 + u(i)^2);
    p = 2*(1/lambda - 1)*(y0-u(i)) + kbar*u(i);
    
    Pbar = [Pbar; p];

end


figure(6)
p1 = plot(ubar, force_curve, 'LineWidth', 2);
hold on 
p2 = plot(u, Pbar, '^', 'LineWidth', 1.5);
hold off
xlabel('Displacement')
ylabel('Force')
legend([p1 p2], 'Modified Newton-Raphson', 'True path')
xlim([0 2])
if kbar == 1
    title('Modified Newton-Raphson, kbar = 1')
else
    title('Modified Newton-Raphson, kbar = 3')
end
