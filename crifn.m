%% Project 1, Exercise 2, crifn
clear;
close all;
clc;

% Coordinates
x1=1.697;
coord=[-x1 -1 0  %1 node number
    0  -1 0  %2
    x1 -1 0  %3
    -x1  1 0  %4
    0   1 0  %5
    x1  1 0  %6
    -x1  0 1  %7
    0   0 1  %8
    x1  0 1];%9

% Degrees of freedom for each node
dof=[1 2 3      %1 node number
    4 5 6      %2
    7 8 9      %3
    10 11 12   %4
    13 14 15   %5
    16 17 18   %6
    19 20 21   %7
    22 23 24   %8
    25 26 27]; %9

% Topology matrix
edof=[1  dof(1,:) dof(7,:)
    2  dof(2,:) dof(7,:)
    3  dof(2,:) dof(8,:)
    4  dof(2,:) dof(9,:)
    5  dof(3,:) dof(9,:)
    6  dof(4,:) dof(7,:)
    7  dof(5,:) dof(7,:)
    8  dof(5,:) dof(8,:)
    9  dof(5,:) dof(9,:)
    10 dof(6,:) dof(9,:)
    11 dof(7,:) dof(8,:)
    12 dof(8,:) dof(9,:)];

% Extract element coordinates
[Ex,Ey,Ez]=coordxtr(edof,coord,dof,2);

enod = [1 2;
        2 3;
        3 4;
        4 5;
        5 6;
        6 7;
        7 8;
        8 9;
        9 10;
        10 11;
        11 12;
        12 13];

% figure(5);
% eldraw3(Ex,Ey,Ez,[1 5 1]);
% hold on

% Essential boundary conditions
bc=[dof(1,:)' [0 0 0]'
    dof(2,:)' [0 0 0]'
    dof(3,:)' [0 0 0]'
    dof(4,:)' [0 0 0]'
    dof(5,:)' [0 0 0]'
    dof(6,:)' [0 0 0]'];

% pause
% global vectors sizes
ndof=max(max(dof));
nelm=max(edof(:,1));

% External load increment
P = zeros(ndof,1);
P([dof(7,3) dof(8,3) dof(9,3)]) = -0.03*[1.5 1 1.5];

% Degrees of freedom to plot (u,v, and w)
plotdof_w = dof(8,3);  % Center top node (displacement w)
plotdof_v = dof(9,3);  % Edge node (displacement v)
plotdof_u = dof(9,1);  % Edge node (displacement u)

Edof=edof;

% Truss properties
ep=[1 1];
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

K = zeros(ndof, ndof);
f = zeros(ndof, 1);
f_int = zeros(ndof, 1);
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

        K = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el = 1:nelm
            %ec = [coord(el, 1:2:end); coord(el, 2:2:end)]';     % Element coordinates, 3x2 matrix
            ec = [Ex(el, 1), Ex(el, 2);
                  Ey(el, 1), Ey(el, 2);
                  Ez(el, 1), Ez(el, 2)];
            ed = extract(edof(el, :), a0);
            
            lambdaNL = stretch1D(ec, ed);
            sg = stress1D(ep, lambdaNL);
            [es, eg] = bar3gsNL(ec, ed, ep, sg);                      % Normal force es and Green's strain eg
            Et = dmat1D(ep, lambdaNL);
            
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
        if iter == 1            % The first iteration we need to calculate dlambda differently due to denominator in c becoming zero
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
                break;          % KANSKE REPETERA KODEN SÅ ATT L MINSKAS DIREKT ISTÄLLET FÖR TOTAL OMSTART?
            end
            
        end
        
        lambda0 = lambda0 + dlambda;
        
        da = da_r + dlambda*da_p;       % Displacement increment
        a0 = a0 + da;
        
        
        % STRESSES AND STRAINS???
        
        
        for el = 1:nelm
            %ec = [coord(enod(el, 1), :); coord(enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(edof(el,:), a0);                                 % Element displacement
            
            lambdaNL = stretch1D(ec, ed);
            sg = stress1D(ep, lambdaNL);
            [es, eg] = bar3gsNL(ec, ed, ep, sg);                      % Normal force es and Green's strain eg
            %Et = dmat1D(ep, lambdaNL);
            
            indx_el = edof(el, 2:end);                                    % Dof for element
            
            % Internal element force vectors
            fe = bar3gf(ec, es, ed);
            
            % Assemble force vector
            f_int(indx_el) = f_int(indx_el) + fe;
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*P;
        %c = da_tot'*da_tot + psi*dlambda_tot^2*(P')*P - l^2;
        
        res(bc(:,1)) = 0;
        % New states
        iter = 2;
    end
    % Accept quantities
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    lambdaplot2(n) = lambda;
    %     uplot(n) = a(24);                 % a(24) is the displacement in u-dir on center top node
    %     vplot(n) = a(25);                 % a(25) is the displacement in v-dir on edge node
    %     wplot(n) = a(27);                 % a(27) is the displacement in w-dir on edge node
    uplot2(n) = a(plotdof_u);
    vplot2(n) = a(plotdof_v);
    wplot2(n) = a(plotdof_w);
    %     plotdof_w = dof(8,3);             % Center top node (displacement w)
    %     plotdof_v = dof(9,3);             % Edge node (dispalcement v)
    %     plotdof_u = dof(9,1);             % Edge node (dispalcement u)
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
p1 = plot(uplot2, lambdaplot2);
hold on
p2 = plot(uplot, lambdaplot);
xlabel('Displacement u')
ylabel('Load parameter, lambda')
legend([p1 p2], 'Task 2', 'Task 1')

figure(3)
p3 = plot(wplot2, lambdaplot2);
hold on
p4 = plot(wplot, lambdaplot);
xlabel('Displacement w')
ylabel('Load parameter, lambda')
legend([p3 p4], 'Task 2', 'Task 1')

figure(4)
p5 = plot(wplot2, vplot2);
hold on
p6 = plot(wplot, vplot);
xlabel('Displacement w')
ylabel('Displacement v')
legend([p5 p6], 'Task 2', 'Task 1')

figure(5)
ed2 = extract(edof, a0);
eldraw3(Ex, Ey, Ez, [1 4 1]);
hold on
eldisp3(Ex, Ey, Ez, ed2, [1 3 0], 1);
legend('Reference configuration')
title('Constrained path-following for non-linear material')

figure(6)
ed3 = extract(edof1, a1);
%eldraw3(Ex, Ey, Ez, [1 4 1]);
eldisp3(Ex, Ey, Ez, ed2, [1 4 1], 1);
hold on
eldisp3(Ex1, Ey1, Ez1, ed3, [1 3 0], 1);
legend('Task 2')
title('Deformed configuration for task 1 and task 2')

Ex2 = Ex;
Ey2 = Ey;
Ez2 = Ez;
a2 = a0;
edof2 = edof;

save('plotvariables_ex2.mat', 'uplot2', 'vplot2', 'wplot2', 'lambdaplot2', 'Ex2', 'Ey2', 'Ez2', 'a2', 'edof2');
