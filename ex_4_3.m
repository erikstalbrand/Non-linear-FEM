%% Exercise 4.3
clear;
close all;
clc;
% Coordinates
coord=[ 0    0  0   %1 node number
        0   10  0   %2
       10   10  0   %3
       10    0  0]; %4

% Degrees of freedom for each node
dof=[1  2       %1 node number
     3  4       %2
     5  6       %3
     7  8];     %4

% Topology matrix
edof=[1  dof(1,:) dof(2,:) dof(3,:)
      2  dof(1,:) dof(3,:) dof(4,:)];
 
% Extract element coordinates
[Ex,Ey,Ez]=coordxtr(edof,coord,dof,3);

% figure(1);
% eldraw3(Ex,Ey,Ez,[1 5 1]);

bc = [1 0;
      2 0;
      3 0;
      5 -1;
      7 -1;
      8 0];
  
disp_nodes = [5 7]';

% global vectors sizes
ndof=max(max(dof));
nelm=max(edof(:,1));

for i = 1:nelm
ec = [Ex(i, :); Ey(i, :)];       % element coordinates
end

%% Material 1
E = 210E3;
v = 0.3;
mpara = [E v];
t = 0.001;

load_steps = 150;
TOL = 1E-5;

flag = 1;               % St. Venant-Kirchhoff assumed
a = zeros(ndof, 1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
f = 0;

bc1 = bc;
scale = 1/80;
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
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [Ex(el, :); Ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D1(flag, mpara, defgrad);
            stress = stressMater2D1(flag, mpara, defgrad);
            Ke = cont2D3te(ec, t, D, ed, stress);
        
            indx = edof(el, 2:end);
        
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        

        da = solveq(Kt, -res, bc);
        a0 = a0 + da;
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [Ex(el, :); Ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D1(flag, mpara, defgrad);
            stress = stressMater2D1(flag, mpara, defgrad);
            fe = cont2D3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
            
        end
        
        bc(:,2) = 0;
        
        res = f_int - f;
        res(bc(:, 1)) = 0;
                
        iter = 2;
    end

    a = a0;

    force_plot(n) = f_int(7);
    disp_plot(n) = a(7);
end

ed = extract(edof, a);
figure(1)
eldraw2(Ex, Ey, [1 2 1])
hold on
eldisp2(Ex, Ey, ed, [1,4,1], 1);

figure(2)
hold on
p1 = plot(-disp_plot, force_plot);

%% Exercise 4.3 material 2
% Coordinates
coord=[ 0    0  0   %1 node number
        0   10  0   %2
       10   10  0   %3
       10    0  0]; %4

% Degrees of freedom for each node
dof=[1  2       %1 node number
     3  4       %2
     5  6       %3
     7  8];     %4

% Topology matrix
edof=[1  dof(1,:) dof(2,:) dof(3,:)
      2  dof(1,:) dof(3,:) dof(4,:)];
 
% Extract element coordinates
[Ex,Ey,Ez]=coordxtr(edof,coord,dof,3);

% figure(1);
% eldraw3(Ex,Ey,Ez,[1 5 1]);

bc = [1 0;
      2 0;
      3 0;
      5 -1;
      7 -1;
      8 0];
  
disp_nodes = [5 7]';

% global vectors sizes
ndof=max(max(dof));
nelm=max(edof(:,1));

for i = 1:nelm
ec = [Ex(i, :); Ey(i, :)];       % element coordinates
end

E = 210E3;
v = 0.3;
mpara = [E v];
t = 0.001;

load_steps = 150;
TOL = 1E-5;

flag = 1;               % St. Venant-Kirchhoff assumed
a = zeros(ndof, 1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
f = 0;

bc1 = bc;
scale = 1/80;
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
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [Ex(el, :); Ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            Ke = cont2D3te(ec, t, D, ed, stress);
        
            indx = edof(el, 2:end);
        
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        

        da = solveq(Kt, -res, bc);
        a0 = a0 + da;
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [Ex(el, :); Ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            fe = cont2D3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
            
        end
        
        bc(:,2) = 0;
        
        res = f_int - f;
        res(bc(:, 1)) = 0;
                
        iter = 2;
    end

    a = a0;

    force_plot(n) = f_int(7);
    disp_plot(n) = a(7);
end

ed = extract(edof, a);
figure(1)
eldraw2(Ex, Ey, [1 2 1])
hold on
eldisp2(Ex, Ey, ed, [1,5,1], 1);
title('Undeformed configuration(green), Deformed linear material(red), Deformed non-linear material(blue)')
xlabel('[mm]')
ylabel('[mm]')

figure(2)
p2 = plot(-disp_plot, force_plot);
xlabel('Displacement [mm]');
ylabel('Force [N]')
title('Force-displacement for dof=7')
legend([p1 p2], 'Linear', 'Non-linear')