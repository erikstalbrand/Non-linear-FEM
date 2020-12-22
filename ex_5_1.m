%% Exercise 5.1

clear;
clc;
close all;

load('geometry_E51.mat');

% material parameters
E = 210E3;
v = 0.3;
mpara = [E v];
t=1;

ndof = max(max(edof));
nelm = length(edof);
for i = 1:nelm
ec = [ex(i, :); ey(i, :)];       % element coordinates
end

load_steps = 100;
TOL = 1E-6;

flag = 1;               % St. Venant-Kirchhoff assumed
a = zeros(ndof, 1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
df = 0.01;
f = zeros(ndof, 1);

disp_nodes = bc(27:end, 1);       % Non-zero boundary condition on nodes 27-31, nodes with displacement
bc1 = bc;
scale = 1/10;
bc1(:, 2) = bc1(:, 2)*scale;

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);

% Displacement controlled Newton-Raphson
for n = 1:load_steps
disp(['Load Step: ', num2str(n)])

bc = bc1;
a0 = a;

res = f_int - f;

iter = 1;

    while max(abs(res)) > TOL || iter == 1
        disp(['Residual: ', num2str(norm(res))])
        
        Kt = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        D = cell(4,1);
        stress = cell(4,1);
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D4ts(ec, ed);
            for j = 1:4
                D{j,1} = dMater2D2(flag, mpara, defgrad{j,1});
                stress{j,1} = stressMater2D2(flag, mpara, defgrad{j,1});
            end
            Ke = cont2D4te(ec, t, D, ed, stress);
        
            indx = edof(el, 2:end);
        
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        
        da = solveq(Kt, -res, bc);
        a0 = a0 + da;
        for el = 1:nelm
            
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D4ts(ec, ed);
            for j= 1:4
                D{j, 1} = dMater2D2(flag, mpara, defgrad{j, 1});
                stress{j, 1} = stressMater2D2(flag, mpara, defgrad{j, 1});
            end
            fe = cont2D4tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe;
            
        end
        
        res = f_int - f;
        res(bc(:, 1)) = 0;
        
        bc(:, 2) = 0;       % Reset boundary conditions
        
        iter = 2;
    end

    a = a0;
    
    force_plot(n+1) = f_int(33);
    disp_plot(n+1) = a(33);
    
end

ed = extract(edof, a);
figure(1)
eldraw2(ex, ey, [1 2 1])
hold on
eldisp2(ex, ey, ed, [1,4,1], 1);

figure(2)
plot(disp_plot, force_plot);
xlim([0 100]);
xlabel('Displacement [mm]')
ylabel('Internal force [N]')
title('Force-displacement curve for dof=33');