%% Project 1, Exercise 2, crifn
clear;
close all;
clc;
load('geom_2020.mat')
%--------------------------------------------------------------%
path_steps = 100;
da_r = zeros(ndof, 1);
da_tot = zeros(ndof, 1);
da_p = zeros(ndof, 1);
lambda = 0;
dlambda = 0;
dlambda_tot = 0;
a = zeros(ndof, 1);             % Global displacement vector
l = 0.1;                        % Sphere radius
TOL = 1E-5;                     % Convergence criteria
psi = 1;

K = zeros(ndof, ndof);          % Global stiffness matrix
f = zeros(ndof, 1);             % Global force vector
f_int = zeros(ndof, 1);         % Internal force vector
lambdaplot2 = zeros(path_steps, 1);

uplot2 = zeros(path_steps, 1);
vplot2 = zeros(path_steps, 1);
wplot2 = zeros(path_steps, 1);

draw = 1;
for n = 1:path_steps
    disp(['Path Step: ', num2str(abs(n))]);
    
    % Iteration quantities
    iter = 1;
    a0 = a;
    lambda0 = lambda;
    dlambda_tot = 0;
    
    res = zeros(ndof, 1);
    da_tot = zeros(ndof, 1);
    
    while norm(res) > TOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);
        
        % Reset quantities
        K = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el = 1:nelm
            ec = [Ex(el, 1), Ex(el, 2);             % Element coordinates
                Ey(el, 1), Ey(el, 2);
                Ez(el, 1), Ez(el, 2)];
            ed = extract(edof(el, :), a0);          % Element displacements
            
            lambdaNL = stretch1D(ec, ed);           % Stretch
            sg = stress1D(ep, lambdaNL);
            es = ep(2)*sg;              % Normal force es
            Et = dmat1D(ep, lambdaNL);              % Non-linear element stiffness
            
            Ke = bar3geNL(ec, ed, ep, es, Et);                        % Element stiffness matrix
            
            indx_el = Edof(el, 2:end);                          % Dof for element
            % Assemble stiffness matrix
            K(indx_el, indx_el) = K(indx_el, indx_el) + Ke;
        end
        
        % Pseudo displacement increments
        da_r = solveq(K, -res, bc);
        da_p = solveq(K, P, bc);
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
        
        lambda0 = lambda0 + dlambda;
        da = da_r + dlambda*da_p;       % Displacement increment
        a0 = a0 + da;
        
        % STRESSES AND STRAINS
        for el = 1:nelm
            ed = extract(edof(el,:), a0);       % Element displacement
            ec = [Ex(el, 1), Ex(el, 2);         % Element coordinates
                  Ey(el, 1), Ey(el, 2);
                  Ez(el, 1), Ez(el, 2)];
            
            lambdaNL = stretch1D(ec, ed);       % Stretch
            sg = stress1D(ep, lambdaNL);
            es = ep(2)*sg;                      % Normal force es
            
            indx_el = edof(el, 2:end);          % Dof for element
            
            % Internal element force vectors
            fe = bar3gf(ec, es, ed);
            
            % Assemble force vector
            f_int(indx_el) = f_int(indx_el) + fe;
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*P;
        
        res(bc(:,1)) = 0;
        iter = 2;
    end
    % Accept quantities
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    lambdaplot2(n) = lambda;
    uplot2(n) = a(plotdof_u);
    vplot2(n) = a(plotdof_v);
    wplot2(n) = a(plotdof_w);
    
    coordUpd = coord;
    %Draws the structure every 20th step
    if(mod(n, 20) == 0 && draw == 1 || n ==1)
        for k=1:3
            coordUpd(:,k) = coordUpd(:,k)+a(k:3:(end+k-3));
            pause(0.01)
        end
        [Exx,Eyy,Ezz]=coordxtr(edof,coordUpd,dof,2);
        clf;
        
        figure(1)
        eldraw3(Exx,Eyy,Ezz,[1 4 1]);
        drawnow;
        title('Deformed configuration task 2')
        hold on
    end
end

load('plotvariables_ex1.mat');

figure(2)
p1 = plot(uplot2, lambdaplot2, 'r');
hold on
p2 = plot(uplot, lambdaplot, 'b');
xlabel('Displacement u')
ylabel('Load parameter, \lambda')
legend([p1 p2], 'Non-linear elastic material model', 'Linear elastic material model')
title('Equilibrium paths')

figure(3)
p3 = plot(wplot2, lambdaplot2, 'r');
hold on
p4 = plot(wplot, lambdaplot, 'b');
xlabel('Displacement w')
ylabel('Load parameter, \lambda')
legend([p3 p4], 'Non-linear elastic material model', 'Linear elastic material model')
title('Equilibrium paths')

figure(4)
p5 = plot(wplot2, vplot2, 'r');
hold on
p6 = plot(wplot, vplot, 'b');
xlabel('Displacement w')
ylabel('Displacement v')
legend([p5 p6], 'Non-linear elastic material model', 'Linear elastic material model')
title('Equilibrium paths')

figure(5)
ed2 = extract(edof, a0);
eldraw3(Ex, Ey, Ez, [1 4 1]);
hold on
eldisp3(Ex, Ey, Ez, ed2, [1 2 1], 1);
legend('Reference configuration')
title('Constrained path-following for non-linear material')

figure(6)
ed3 = extract(edof1, a1);
%eldraw3(Ex, Ey, Ez, [1 4 1]);
eldisp3(Ex, Ey, Ez, ed2, [1 2 1], 1);
hold on
eldisp3(Ex1, Ey1, Ez1, ed3, [1 4 1], 1);
legend('Non-linear elastic material model')
title('Deformed configuration for the linear and non-linear elastic material model')

Ex2 = Ex;
Ey2 = Ey;
Ez2 = Ez;
a_2 = a0;
edof2 = edof;

save('plotvariables_ex2.mat', 'uplot2', 'vplot2', 'wplot2', 'lambdaplot2', 'Ex2', 'Ey2', 'Ez2', 'a_2', 'edof2');