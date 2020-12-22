%% Exercise 4.4 Chrisfield

clear;
clc;
close all;

load('arch_mesh.mat');

load_steps = 360;
TOL = 1E-6;
cTOL = 1E-6;

%--------------------------------------------
flag = 1;
a = zeros(nrdof, 1);
Kt = zeros(nrdof, nrdof);
f_int = zeros(nrdof, 1);
f = zeros(nrdof, 1);

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);
%--------------------------------------------
v = 0.3;
E = 210E3;
mpara = [E, v];
t = 1;

da_r = zeros(nrdof, 1);
da_tot = zeros(nrdof, 1);
da_p = zeros(nrdof, 1);
lambda = 0;
dlambda = 0;
dlambda_tot = 0;
psi = 0;
l = 5;
nelm = nrelem;
c = 0;

res = zeros(nrdof, 1);

for n = 1:load_steps
    disp(['Load Step: ', num2str(n)])
        
    % Iteration quantities
    lambda0 = lambda;
    a0 = a;
    
    % Reset variables
    da_tot = zeros(nrdof, 1);
    dlambda_tot = 0;
    iter = 1;
    res = f_int - lambda0*P;
    while max(abs(res)) > TOL || abs(c) > cTOL || iter == 1
        disp(['Residual: ', num2str(norm(res))])
        
        Kt = zeros(nrdof, nrdof);
        f_int = zeros(nrdof, 1);
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            Ke = cont2D3te(ec, t, D, ed, stress);
            
            indx = edof(el, 2:end);
            
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        % Pseudo displacement increments
        da_r = solveq(Kt, -res, bc);
        da_p = solveq(Kt, P, bc);
        da_tot = a0 - a;
        dlambda_tot = lambda0 - lambda;
        
        % Load parameter
        if iter == 1            % The first iteration we need to calculate dlambda differently
            if n > 1
                da_n = a - a_old;     % Increment between the last two converged equilibrium points
                s = sign(da_n'*da_p);
                dlambda = s*(l/sqrt(da_p'*da_p + psi*(P')*P));
            else
                dlambda = l/sqrt(da_p'*da_p + psi*(P')*P);     % s = 1 in first iteration
            end
        else
            % Load parameter
            a1 = da_p'*da_p + psi*(P')*P;
            a2 = 2*da_p'*(da_tot+da_r) + 2*psi*dlambda_tot*(P')*P;
            a3 = (da_tot+da_r)'*(da_tot+da_r)+psi*dlambda_tot^2*(P')*P-l^2;
            a4 = da_tot'*(da_tot+da_r);
            a5 = da_tot'*da_p;
            
            eq = [a1 a2 a3];
            dlambdas = roots(eq);
            
            if isreal(dlambdas) == 1
                if a4 + a5*dlambdas(1) > a4 + a5*dlambdas(2)
                    dlambda = dlambdas(1);
                else
                    dlambda = dlambdas(2);
                end
            else
                break;          % If complex dlambda, reduce l and restart
            end
        end

        % Update lambda
        lambda0 = lambda0 + dlambda;
        
        % Displacement increment
        da = da_r + dlambda*da_p;       
        a0 = a0 + da;
        
        % Stresses and strains
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            fe = cont2D3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*P;
        c = da_tot'*da_tot + psi*dlambda_tot^2*(P')*P - l^2;
        
        res(bc(:,1)) = 0;
        iter = 2;
    end
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    force_plot(n) = -P(96)*lambda;
    disp_plot(n) = a(96);
    
end

z = zeros(load_steps, 1);
figure(1)
p1 = plot(disp_plot, force_plot);
hold on
grid on
xlabel('Displacement[mm]')
ylabel('Force[N]')
title('Force-displacement curve for dof=96')

figure(2)
ed = extract(edof, a);
eldraw2(ex, ey, [1 4 1]);
hold on
grid on
eldisp2(ex, ey, ed, [1 2 1], 1);

figure(6)
plot3(z, disp_plot, force_plot)
hold on
grid on
xlabel('x-Displacement[mm]')
ylabel('y-Displacement[mm]')
zlabel('Force[N]')

%% Exercise 4.4 Chrisfield

% clear;
% clc;
% close all;
% load('arch_mesh.mat');

load_steps = 250;
TOL = 1E-6;
cTOL = 1E-6;

%--------------------------------------------
flag = 1;
a = zeros(nrdof, 1);
Kt = zeros(nrdof, nrdof);
f_int = zeros(nrdof, 1);
f = zeros(nrdof, 1);

force_plot = zeros(load_steps, 1);
disp_plot_vert = zeros(load_steps, 1);
disp_plot_hor = zeros(load_steps, 1);
%--------------------------------------------
v = 0.3;
E = 210E3;
mpara = [E, v];
t = 1;

da_r = zeros(nrdof, 1);
da_tot = zeros(nrdof, 1);
da_p = zeros(nrdof, 1);
lambda = 0;
dlambda = 0;
dlambda_tot = 0;
psi = 0;
l = 5;
nelm = nrelem;
c = 0;
P(95) = -0.1;

res = zeros(nrdof, 1);

for n = 1:load_steps
    disp(['Load Step: ', num2str(n)])
        
    % Iteration quantities
    lambda0 = lambda;
    a0 = a;
    
    % Reset variables
    da_tot = zeros(nrdof, 1);
    dlambda_tot = 0;
    iter = 1;
    res = f_int - lambda0*P;
    while max(abs(res)) > TOL || abs(c) > cTOL || iter == 1
        disp(['Residual: ', num2str(norm(res))])
        
        Kt = zeros(nrdof, nrdof);
        f_int = zeros(nrdof, 1);
        
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            Ke = cont2D3te(ec, t, D, ed, stress);
            
            indx = edof(el, 2:end);
            
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        
        
        % Pseudo displacement increments
        da_r = solveq(Kt, -res, bc);
        da_p = solveq(Kt, P, bc);
        da_tot = a0 - a;
        dlambda_tot = lambda0 - lambda;
        
        % Load parameter
        if iter == 1            % The first iteration we need to calculate dlambda differently
            if n > 1
                da_n = a - a_old;     % Increment between the last two converged equilibrium points
                s = sign(da_n'*da_p);
                dlambda = s*(l/sqrt(da_p'*da_p + psi*(P')*P));
            else
                dlambda = l/sqrt(da_p'*da_p + psi*(P')*P);     % s = 1 in first iteration
            end
        else
            % Load parameter
            a1 = da_p'*da_p + psi*(P')*P;
            a2 = 2*da_p'*(da_tot+da_r) + 2*psi*dlambda_tot*(P')*P;
            a3 = (da_tot+da_r)'*(da_tot+da_r)+psi*dlambda_tot^2*(P')*P-l^2;
            a4 = da_tot'*(da_tot+da_r);
            a5 = da_tot'*da_p;
            
            eq = [a1 a2 a3];
            dlambdas = roots(eq);
            
            if isreal(dlambdas) == 1
                if a4 + a5*dlambdas(1) > a4 + a5*dlambdas(2)
                    dlambda = dlambdas(1);
                else
                    dlambda = dlambdas(2);
                end
            else
                break;          % If complex dlambda, reduce l and restart
            end
        end

        % Update lambda
        lambda0 = lambda0 + dlambda;
        
        % Displacement increment
        da = da_r + dlambda*da_p;       
        a0 = a0 + da;
        
        % Stresses and strains
        for el = 1:nelm
            ed = extract(edof(el, :), a0);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            fe = cont2D3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*P;
        c = da_tot'*da_tot + psi*dlambda_tot^2*(P')*P - l^2;
        
        res(bc(:,1)) = 0;
        iter = 2;
    end
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    force_plot(n) = -P(96)*lambda;
    disp_plot_vert(n) = a(96);
    disp_plot_hor(n) = a(95);
end

figure(1)
p2 = plot(disp_plot_vert, force_plot);
hold on
grid on
legend([p1 p2], 'Ideal system', 'Perturbed system');

figure(4)
plot(disp_plot_hor, force_plot);
hold on
grid on
xlabel('Displacement[mm]')
ylabel('Force[N]')
title('Force-displacement curve for dof=95')

figure(5)
ed = extract(edof, a);
eldraw2(ex, ey, [1 4 1]);
hold on
grid on
eldisp2(ex, ey, ed, [1 2 1], 1);

figure(6)
plot3(disp_plot_hor, disp_plot_vert, force_plot)
grid on

