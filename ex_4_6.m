%% Exercise 4.6

clear;
clc;
close all;

load('ring_mesh.mat');

load_steps = 100;
TOL = 1E-3;

ndof = max(max(edof));
nelm=max(edof(:,1));
E = 210E3;
v = 0.3;
mpara = [E, v];
t = 1;

flag = 1;               % St. Venant-Kirchhoff assumed
a = zeros(ndof, 1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
f = zeros(ndof, 1);

disp_nodes = bc([1, 3], 1);       % Non-zero boundary condition on nodes 1 and 3, nodes with displacement
bc1 = bc;
scale = 1/5;
bc1(:, 2) = bc1(:, 2)*scale;

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);
res = zeros(ndof, 1);
% Displacement controlled Newton-Raphson
for n = 1:load_steps
disp(['Load Step: ', num2str(n)])

bc = bc1;
a0 = a;

% res = f_int - f;
res = res*0;
iter = 1;

    while max(abs(res)) > TOL || iter == 1
        disp(['Residual: ', num2str(norm(res))])
        
        Kt = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = contAxi3ts(ec, ed);
            D = dMaterAxi2(flag, mpara, defgrad);
            stress = stressMaterAxi2(flag, mpara, defgrad);
            Ke = contAxi3te(ec, t, D, ed, stress);
        
            indx = edof(el, 2:end);
        
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        da = solveq(Kt, -res, bc);
        a0 = a0 + da;
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = contAxi3ts(ec, ed);
            D = dMaterAxi2(flag, mpara, defgrad);
            stress = stressMaterAxi2(flag, mpara, defgrad);
            fe = contAxi3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
            
        end
        
        res = f_int - f;
        res(bc(:, 1)) = 0;
        
        bc(:, 2) = 0;       % Reset boundary conditions
        
        iter = 2;
    end

    a = a0;
    
    force_plot(n) = min(f_int(disp_nodes(:, 1)));
    disp_plot(n) = max(a(disp_nodes(:, 1)));
    
end

ed = extract(edof, a);
figure(1)
eldraw2(ex, ey, [1 2 1])
hold on
eldisp2(ex, ey, ed, [1,4,1], 1);
title('Original(green) and deformed mesh(red)');
xlabel('[mm]')
ylabel('[mm]')

figure(2)
plot(disp_plot, force_plot);
xlim([0 100]);
xlabel('Displacement [mm]')
ylabel('Force [N]')
title('Force-displacement curve for dof=45');
xlim([0 20]);
