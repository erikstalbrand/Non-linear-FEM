%% Exercise 5.2 6-node isoparametric elements
clear;
clc;
close all;

load('geom_E52.mat');

load_steps = 10;
TOL = 1E-3;

% material parameters
E = 210E3;
v = 0.3;
mpara = [E v];
t = 1;

ndof = max(max(edof));
nelm = length(edof);
for i = 1:nelm
ec = [ex(i, :); ey(i, :)];       % element coordinates
end
exx = ex;
eyy = ey;

flag = 3;               % Cauchy stress
a = zeros(ndof, 1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
f = zeros(ndof, 1);
bc1 = bc;
da = zeros(ndof, 1);

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);

D = cell(3, 1);
stress = cell(3, 1);
defgrad = cell(3, 1);

c = cell(3, 1);
defgrad_old = cell(3, nelm);
defgrado = cell(3, 1);
for el = 1:nelm
    defgrad_old{1, el} = [1 0 0 1]';
    defgrad_old{2, el} = [1 0 0 1]';
    defgrad_old{3, el} = [1 0 0 1]';
end

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
            
            for j = 1:3
                c{j, 1} = defgrad_old{j, el};
            end
            defgrad = cont2D6us(ec, edinc, c);
            for j = 1:3
                D{j,1} = dMater2D2(2, mpara, defgrad{j,1});
                stress{j,1} = stressMater2D2(flag, mpara, defgrad{j,1});
            end
            
            Ke = cont2D6ue(ec, t, D, stress);
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
            
            for j = 1:3
                c{j, 1} = defgrad_old{j, el};
            end
            defgrad = cont2D6us(ec, edinc, c);
            for j = 1:3
                D{j,1} = dMater2D2(2, mpara, defgrad{j,1});
                stress{j,1} = stressMater2D2(flag, mpara, defgrad{j,1});
            end
            fe = cont2D6uf(ec, t, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe;
            for i = 1:3
            defgrad_old{i, el} = defgrad{i, 1}';              % Sparar defgrad för att använda i nästa loop(updated Lagrangian)
            end
        end
        
        res = f_int - f;
        res(bc(:, 1)) = 0;
        
        bc(:, 2) = 0;       % Reset boundary conditions
        
        iter = 2;
    end

    a = a0;
    
    force_plot(n+1) = f_int(313);
    disp_plot(n+1) = a(313);
    
end

ed = extract(edof, a);
figure(1)
eldraw2(exx, eyy, [1 2 1])
hold on
eldisp2(exx, eyy, ed, [1,4,1], 1);
title('Original(green) and deformed mesh(red)');
xlabel('[mm]')
ylabel('[mm]')

figure(2)
plot(disp_plot, force_plot);
xlim([0 100]);
title('Force-displacement curve for dof=313');
xlabel('Displacement [mm]')
ylabel('Internal force [N]')
