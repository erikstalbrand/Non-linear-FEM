%% Exercise 5.2


clear;
clc;
close all;

load('geometry_E51.mat');

%Generate data for total dof and number of elements
ndof = max(max(edof));
nelm = length(edof);
for i = 1:nelm
ec = [ex(i, :); ey(i, :)];       % element coordinates
end

%Initial values
defgrad_old_global = cell(3,nelm);
for el = 1:nelm
    defgrad_old_global{1,el} = [1 0 0 1];
    defgrad_old_global{2,el} = [1 0 0 1];
    defgrad_old_global{3,el} = [1 0 0 1];
end
load_steps = 100;
TOL = 1E-3;

% material parameters
E = 210E3;
v = 0.3;
mpara = [E v];
t=1;



%Useful data
flagstress = 3;     %Cauchystress
flagmater = 2;      %Pushforward

%Useful vectors
a = zeros(ndof, 1);
da = zeros(ndof,1);
Kt = zeros(ndof, ndof);
f_int = zeros(ndof, 1);
f = 0;
S = 0;

%Displacement control
disp_nodes = bc(27:end, 1);       % Non-zero boundary condition on nodes 27-31, nodes with displacement
bc1 = bc;
scale = 10;
bc1(:, 2) = bc1(:, 2)/scale;

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
        %res = f-fint;
        disp(['Residual: ', num2str(norm(res))])

        %Reset K and fint when iterating
        Kt = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        D = cell(4,1);
        stress = cell(4,1);
        defgrad_old = cell(4,1);
        for el = 1:nelm
            edinc = extract(edof(el,:),da);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D6us(ec, edinc,defgrad_old_global{:,el});
            if iter==1
                for j = 1:3
                    D{j,1} = dMater2D1(flagmater, mpara, defgrad{j,1});
                    stress{j,1} = stressMater2D1(flagstress, mpara, defgrad{j,1});
                end
            end
                
            Ke = cont2D6ue(ec, t, D, ed, stress);
        
            %Assemble K
            indx = edof(el, 2:end);
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        %Calculate and update displacement
        da = solveq(Kt, -res, bc);
        a0 = a0 + da;
        
        %Update coordinates
        ed = extract(edof, da);
        ex = ex + ed(:,1:2:end);
        ey = ey + ed(:,1:2:end);
            
        for el = 1:nelm
            ec = [ex(el, :); ey(el, :)];
            edinc = extract(edof(el,:),da);
            
            defgrad = cont2D6us(ec, edinc,defgrad_old_global{el,:});
            for j = 1:3
                D{j,1} = dMater2D1(flagmater, mpara, defgrad{j,1});
                stress{j,1} = stressMater2D1(flagstress, mpara, defgrad{j,1});
            end
            defgrad_old_global(el,:) = defgrad;
            fe = cont2D6uf(ec, t, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe;
            
        end
        %Only apply displacement in the first iteration
        bc(:,2) = 0;
        
        %Out of balance forces
        res = f_int - f;
        res(bc(:, 1)) = 0;
                
        iter = 2;
    end

    a = a0;
    
    
    force_plot(n) = min(f_int(disp_nodes(:, 1)));
    disp_plot(n) = min(a(disp_nodes(:, 1)));
    
    
end

ed = extract(edof, a);
figure(1)
eldraw2(ex, ey, [1 2 1])
hold on
eldisp2(ex, ey, ed, [1,4,1], 1);

figure(2)
plot(disp_plot, force_plot);




