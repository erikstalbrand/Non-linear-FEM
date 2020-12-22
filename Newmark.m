%% Exercise 7.1 Newmark algorithm, df = 200, dt = 0.0001

clear;
clc;
close all;

load('geometry_E71.mat')

load_steps = 100;
TOL = 1E-3;

rho = 7800e-9;
beta = 1/4;
gamma = 2*beta;

flag = 1;                         % Cauchy stress
Kt = zeros(nrdof, nrdof);
f_int = zeros(nrdof, 1);
f = zeros(nrdof, 1);

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);

d1 = 0;
M = zeros(nrdof, nrdof);
C = d1*M;
t = 1;
E = 210E3;
vp = 0.3;
mpara = [E, vp];

dt = 0.0001;   % Time-step
tmax = dt*load_steps;    % Max time

a = zeros(nrdof, 1);
v = zeros(nrdof, 1);
acc = zeros(nrdof, 1);
load = zeros(nrdof,1);
res_eff = zeros(nrdof, 1);

load(find(P)) = -2000000;    % Loading rate, dt*load = df
count = 0;

nelm = nrelem;

% Mass matrix
for el = 1:nelm
    ec = [ex(el, :); ey(el, :)];
    Me = cont2D3m(ec, t, rho);
    indx = edof(el, 2:end);
    M(indx, indx) = M(indx, indx) + Me;
end

% Force controlled Newton-Raphson
for n = 0:dt:tmax
    count = count + 1;
    disp(['LOAD STEP:', num2str(count)]);
    % New load level    
    
    f = n*load;       % f_n+1
%     force_plot(count) = -f(96);
    
    % Predictor step with same acceleration, 7.19
    a = a + dt*v + (dt^2)*(1/2)*acc;
    v = v + dt*acc;

    res_eff = f_int - f + M*acc;
    iter = 1;
    
    while norm(res_eff) > TOL || iter == 1
%         disp(['Residual: ', num2str(norm(res))])
%         disp(['Effective Residual: ', num2str(norm(res_eff))])
        
        % Reset stiffness matrix and force vector
        Kt = zeros(nrdof, nrdof);
        f_int = zeros(nrdof, 1);
        
        for el = 1:nelm
            ed = extract(edof(el, :), a);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            Ke = cont2D3te(ec, t, D, ed, stress);
            
            indx = edof(el, 2:end);
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        K_eff = Kt + (1/(beta*dt^2))*M;            % 7.16

        da = solveq(K_eff, -res_eff, bc);
        
        %  Displacement, velocity and acceleration vectors
        a = a + da;
        v = v + ((gamma)/(beta*dt))*da;
        acc = acc + (1/(beta*dt^2))*da;

        ed = extract(edof, a);
        % Stresses and strains
        for el = 1:nelm
%             ed = extract(edof(el, :), a);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed(el, :));
            stress = stressMater2D2(flag, mpara, defgrad);
            fe = cont2D3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
        end
        
        % Out of balance forces

        res_eff = f_int + M*acc - f;
        res_eff(bc(:, 1)) = 0;
        disp(['Effective Residual: ', num2str(norm(res_eff))])
        
        iter = 2;
    end
    
    force_plot(count) = -f(96);
    disp_plot(count) = -a(96);

    ed = extract(edof, a);
    figure(1)
    drawnow update
    clf(figure(1));
%     eldraw2(ex, ey, [1 2 1], edof(:,1))
    eldraw2(ex,ey,[1 2 1]);
    eldisp2(ex, ey, ed, [1 4 1], 10);
    title('df: 200N, dt: 0.0001 s');
end
figure(2)
plot(disp_plot, force_plot)
title('Force-displacement curve for dof=96')
xlabel('Displacement[mm]');
ylabel('Force[N]');
legend('df: 200N    dt: 0.0001 s')

%% Exercise 7.1 Newmark algorithm, df = 2000, dt = 1e-5

load('geometry_E71.mat')

load_steps = 100;
TOL = 1E-3;

rho = 7800e-9;
beta = 1/4;
gamma = 2*beta;

flag = 1;                         % Cauchy stress
Kt = zeros(nrdof, nrdof);
f_int = zeros(nrdof, 1);
f = zeros(nrdof, 1);

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);

d1 = 0;
M = zeros(nrdof, nrdof);
C = d1*M;
t = 1;
E = 210E3;
vp = 0.3;
mpara = [E, vp];

dt = 1e-5;   % Time-step
tmax = dt*load_steps;    % Max time

a = zeros(nrdof, 1);
v = zeros(nrdof, 1);
acc = zeros(nrdof, 1);
load = zeros(nrdof,1);
res_eff = zeros(nrdof, 1);

load(find(P)) = -200000000;    % Loading rate, dt*load = df
count = 0;

nelm = nrelem;

% Mass matrix
for el = 1:nelm
    ec = [ex(el, :); ey(el, :)];
    Me = cont2D3m(ec, t, rho);
    indx = edof(el, 2:end);
    M(indx, indx) = M(indx, indx) + Me;
end

% Force controlled Newton-Raphson
for n = 0:dt:tmax
    count = count + 1;
    disp(['LOAD STEP:', num2str(count)]);
    % New load level    
    
    f = n*load;       % f_n+1
%     force_plot(count) = -f(96);
    
    % Predictor step with same acceleration, 7.19
    a = a + dt*v + (dt^2)*(1/2)*acc;
    v = v + dt*acc;

    res_eff = f_int - f + M*acc;
    iter = 1;
    
    while norm(res_eff) > TOL || iter == 1
%         disp(['Residual: ', num2str(norm(res))])
%         disp(['Effective Residual: ', num2str(norm(res_eff))])
        
        % Reset stiffness matrix and force vector
        Kt = zeros(nrdof, nrdof);
        f_int = zeros(nrdof, 1);
        
        for el = 1:nelm
            ed = extract(edof(el, :), a);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            Ke = cont2D3te(ec, t, D, ed, stress);
            
            indx = edof(el, 2:end);
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        K_eff = Kt + (1/(beta*dt^2))*M;            % 7.16

        da = solveq(K_eff, -res_eff, bc);
        
        %  Displacement, velocity and acceleration vectors
        a = a + da;
        v = v + ((gamma)/(beta*dt))*da;
        acc = acc + (1/(beta*dt^2))*da;

        ed = extract(edof, a);
        % Stresses and strains
        for el = 1:nelm
%             ed = extract(edof(el, :), a);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed(el, :));
            stress = stressMater2D2(flag, mpara, defgrad);
            fe = cont2D3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
        end
        
        % Out of balance forces

        res_eff = f_int + M*acc - f;
        res_eff(bc(:, 1)) = 0;
        disp(['Effective Residual: ', num2str(norm(res_eff))])
        
        iter = 2;
    end
    
    force_plot(count) = -f(96);
    disp_plot(count) = -a(96);

    ed = extract(edof, a);
    figure(3)
    drawnow update
    clf(figure(3));
%     eldraw2(ex, ey, [1 2 1], edof(:,1))
    eldraw2(ex,ey,[1 2 1]);
    eldisp2(ex, ey, ed, [1 4 1], 10);
    title('df: 2000N, dt: 1e-05 s');
    
end
figure(4)
plot(disp_plot, force_plot)
title('Force-displacement curve for dof=96')
xlabel('Displacement[mm]');
ylabel('Force[N]');
legend('df: 2000N   dt: 1e-05 s')

