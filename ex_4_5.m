%% Exercise 4.5, St. Venant-Kirchhoff model
clear;
clc;
close all;

load('data.mat');

load_steps = 100;
TOL = 1E-6;

flag = 3;               % St. Venant-Kirchhoff assumed
a = zeros(ndof, 1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
f = zeros(ndof, 1);

disp_nodes = bc(27:end, 1);       % Non-zero boundary condition on nodes 27-31, nodes with displacement
bc1 = bc;
scale = 1/10;
bc1(:, 2) = bc1(:, 2)*scale;

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);
% defgrad = zeros(4, 1);

defgrad_old = zeros(nelm, 4);
defgrad_old(:, [1, 4]) = 1;
da = zeros(ndof, 1);
res = zeros(ndof, 1);

exx = ex;
eyy = ey;
% Displacement controlled Newton-Raphson
for n = 1:load_steps
disp(['Load Step: ', num2str(n)])

bc = bc1;
a0 = a;

res = f_int - f;

iter = 1;

    while norm(res) > TOL || iter == 1
        disp(['Residual: ', num2str(norm(res))])
        
        Kt = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el = 1:nelm
            edinc = extract(edof(el, :), da);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3us(ec, edinc, defgrad_old(el, :));
            D = dMater2D1(2, mpara, defgrad);
            stress = stressMater2D1(flag, mpara, defgrad);
            Ke = cont2D3ue(ec, t, D, stress);
        
            indx = edof(el, 2:end);
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        da = solveq(Kt, -res, bc);
        a0 = a0 + da;
        ed = extract(edof, da);
        ex = ex + ed(:, 1:2:end);
        ey = ey + ed(:, 2:2:end);
        
        for el = 1:nelm
            edinc = extract(edof(el, :), da);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3us(ec, edinc, defgrad_old(el, :));
            stress = stressMater2D1(flag, mpara, defgrad);
            fe = cont2D3uf(ec, t, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
            
            defgrad_old(el, :) = defgrad';              % Sparar defgrad för att använda i nästa loop(updated Lagrangian)
        end
        
        res = f_int - f;
        res(bc(:, 1)) = 0;  % Res is zero where boundary conditions are applied
        
        bc(:, 2) = 0;       % Reset boundary conditions to not apply displacement each iteration
        
        iter = 2;
    end
    % Accept quantities
    a = a0;
    
    % Vectors for plotting
    force_plot(n+1) = f_int(369);
    disp_plot(n+1) = a(369);
end

ed = extract(edof, a);
figure(1)
eldraw2(exx, eyy, [1 2 1])
hold on
eldisp2(exx, eyy, ed, [1,4,1], 1);
title('Original(green), deformed linear mesh(red) and deformed non-linear mesh(blue)');
xlabel('[mm]')
ylabel('[mm]')

figure(2)
p1 = plot(disp_plot, force_plot);
hold on
xlim([0 100]);
title('Force-displacement curve for dof=369');
xlabel('Displacement [mm]')
ylabel('Internal force [N]')

%% Exercise 4.5, Neo-Hookean model
% clear;
% clc;
% close all;

load('data.mat');

load_steps = 100;
TOL = 1E-3;

flag = 3;               % St. Venant-Kirchhoff assumed
a = zeros(ndof, 1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
f = zeros(ndof, 1);

disp_nodes = bc(27:end, 1);       % Non-zero boundary condition on nodes 27-31, nodes with displacement
bc1 = bc;
scale = 1/10;
bc1(:, 2) = bc1(:, 2)*scale;

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);
% defgrad = zeros(4, 1);

defgrad_old = zeros(nelm, 4);
defgrad_old(:, [1, 4]) = 1;
da = zeros(ndof, 1);
res = zeros(ndof, 1);

% Displacement controlled Newton-Raphson
for n = 1:load_steps
disp(['Load Step: ', num2str(n)])

bc = bc1;
a0 = a;

res = f_int - f;

iter = 1;

    while norm(res) > TOL || iter == 1
        disp(['Residual: ', num2str(norm(res))])
        
        Kt = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el = 1:nelm
            edinc = extract(edof(el, :), da);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3us(ec, edinc, defgrad_old(el, :));
            D = dMater2D2(2, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            Ke = cont2D3ue(ec, t, D, stress);
        
            indx = edof(el, 2:end);
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        da = solveq(Kt, -res, bc);
        a0 = a0 + da;
        ed = extract(edof, da);
        ex = ex + ed(:, 1:2:end);
        ey = ey + ed(:, 2:2:end);        
        
        for el = 1:nelm
            edinc = extract(edof(el, :), da);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3us(ec, edinc, defgrad_old(el, :));
            stress = stressMater2D2(flag, mpara, defgrad);
            fe = cont2D3uf(ec, t, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
            
            defgrad_old(el, :) = defgrad';              % Sparar defgrad för att använda i nästa loop(updated Lagrangian)
        end
        
        res = f_int - f;
        res(bc(:, 1)) = 0;  % Res is zero where boundary conditions are applied
        
        bc(:, 2) = 0;       % Reset boundary conditions to not apply displacement each iteration
        
        iter = 2;
    end
    % Accept quantities
    a = a0;
    
    % Vectors for plotting
    force_plot(n+1) = f_int(369);
    disp_plot(n+1) = a(369);
end

ed = extract(edof, a);
figure(1)
eldraw2(exx, eyy, [1 2 1])
hold on
eldisp2(exx, eyy, ed, [1,5,1], 1);

figure(2)
p2 = plot(disp_plot, force_plot);
legend([p1 p2], 'Linear', 'Non-linear');
