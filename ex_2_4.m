%% Exercise 2.4 a) Fully linearized, ideal system

clear;
clc;

path_steps = 300;

fi = 65;
x0A = -cosd(fi);        %*
x0B = 0;
x0C = cosd(fi);         %*
y0B = sind(fi);
z0B = 0;
alpha1 = 0.55*y0B;
alpha2 = 0;
alpha = [alpha1 alpha2];
TOL = 1E-5;
cTOL = 1E-5;

P = zeros(1, path_steps);
ybar = zeros(1, path_steps);

x = zeros(3, 1);
y = zeros(3, 1);
z = zeros(3, 1);
x(1) = x0A;
%x(2) = x0B;
y(2) = y0B;
%z(2) = z0B;
x(3) = x0C;
coord = [x y z];

ndof = 11;                      % No zero entries

f = zeros(ndof, 1);
f_int = zeros(ndof, 1);
da = zeros(ndof, 1);
a = zeros(ndof, 1);
K = zeros(ndof, ndof);

Edof=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9];             % Dof for each element
Edof_spring = [1 6 10];
%2 5 11];                       % Dof for each spring
Enod = [1 2;
    2 3];
nelm = 2;
ep = [1 1];
bc = [1 0;
    2 0;
    3 0;
    7 0;
    8 0;
    9 0;
    10 0;
    11 0];

da_r = zeros(ndof, 1);
da_tot = zeros(ndof, 1);        % scalar instead of vector?
da_p = zeros(ndof, 1);
lambda = 0;
dlambda = 0;
dlambda_tot = 0;
c = 0;
psi = 1;
p = zeros(ndof, 1);
p(5) = -0.1;
res = zeros(ndof, 1);
l = 0.018;                       % radius of sphere

for n = 1:path_steps
    disp(['Path Step: ', num2str(abs(n))]);
    iter = 1;
    
    a0 = a;
    lambda0 = lambda;
    
    dlambda_tot = 0;
    res = zeros(ndof, 1);
    da_tot = zeros(ndof, 1);
    
    while norm(res) > TOL && abs(c) > cTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);
        
        K = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el=1:nelm
            % Calculate stiffness matrix
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            Ke = bar3ge(ec, ed, ep, es);                                  % Element stiffness matrix
            Ks = bar1e(ep, alpha(el));                                    % Spring stiffness matrix
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            % Assemble stiffness matrix
            K(indx_el, indx_el) = K(indx_el, indx_el) + Ke;
            K(indx_s, indx_s) = K(indx_s, indx_s) + Ks;
        end
        
        % Pseudo displacement increments
        da_r = solveq(K, -res, bc);
        da_p = solveq(K, p, bc);
        da_tot = a0 - a;                    % eq. 2.127, current da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, current dlambda_tot
        
        % Load parameter
        if iter == 1            % The first iteration we need to calculate dlambda differently due to denominator in c becoming zero
            if n > 1
                da_n = a - a_old;     % ???
                s = sign(da_n'*da_p);
                dlambda = s*(l/sqrt(da_p'*da_p + psi*(p')*p));
            else
                dlambda = l/sqrt(da_p'*da_p + psi*(p')*p);     % s = 1 i första iterationen
            end
        else
            % Load parameter
            dlambda = -(c + 2*da_tot'*da_r)/(2*da_tot'*da_p + 2*psi*(p')*p*dlambda_tot);
        end
        
        lambda0 = lambda0 + dlambda;
        
        % Displacement vector
        da = da_r + dlambda*da_p;
        a0 = a0 + da;               % New displacement vector
        
        % Stresses and strains
        for el=1:nelm
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            
            % Internal element force vectors
            fe = bar3gf(ec, es, ed);
            fe_s = bar1f(alpha(el), eds);
            
            % Assemble force vector
            f_int(indx_el) = f_int(indx_el) + fe;
            f_int(indx_s) = f_int(indx_s) + fe_s;
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*p;
        c = da_tot'*da_tot + psi*dlambda_tot^2*(p')*p - l^2;
        %c = -(2*da_tot'*da_r + 2*da_tot'*da_p*dlambda + 2*psi*(p')*p*dlambda_tot*dlambda);
        
        res(bc(:,1)) = 0;
        % New states
        iter = 2;
    end
    
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    P(n) = -lambda*p(5);
    ybar(n) = y0B + a(5);               % a(5) is the displacement in y-dir on node 2
end

figure(1)
p1 = plot(ybar, P, 'r');              % Plotting force-displacement for
hold on
xlabel('y')
ylabel('P')
title('Fully linearized force-displacement curve')
xlim([-2 1]);
ylim([-0.5 3]);
%legend('Ideal system')


%% 2.4 a) Fully linearized, tilted system

path_steps = 300;

fi = 65;
x0A = -cosd(fi);        %
x0C = cosd(fi);         %*
y0B = sind(fi);
z0B = 0.1;
alpha1 = 0.55*y0B;
alpha2 = 0;
alpha = [alpha1 alpha2];
TOL = 1E-5;
cTOL = 1E-5;

P = zeros(1, path_steps);
ybar = zeros(1, path_steps);

x = zeros(3, 1);
y = zeros(3, 1);
z = zeros(3, 1);
x(1) = x0A;
y(2) = y0B;
z(2) = z0B;
x(3) = x0C;
coord = [x y z];

ndof = 11;                  % No zero entries

f = zeros(ndof, 1);
f_int = zeros(ndof, 1);
da = zeros(ndof, 1);
a = zeros(ndof, 1);
K = zeros(ndof, ndof);

Edof=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9];          % Dof for each element
Edof_spring = [1 6 10];
%2 5 11];        % Dof for each spring
Enod = [1 2;
    2 3];
nelm = 2;
ep = [1 1];
bc = [1 0;
    2 0;
    3 0;
    7 0;
    8 0;
    9 0;
    10 0;
    11 0];

da_r = zeros(ndof, 1);
da_tot = zeros(ndof, 1);        % scalar instead of vector?
da_p = zeros(ndof, 1);
lambda = 0;
dlambda = 0;
dlambda_tot = 0;
c = 0;
psi = 1;
p = zeros(ndof, 1);
p(5) = -0.1;
res = zeros(ndof, 1);
l = 0.018; % radius on sphere

for n = 1:path_steps
    disp(['Path Step: ', num2str(abs(n))]);
    iter = 1;
    
    a0 = a;
    lambda0 = lambda;
    
    dlambda_tot = 0;
    res = zeros(ndof, 1);
    da_tot = zeros(ndof, 1);
    
    while norm(res) > TOL && abs(c) > cTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);
        
        K = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el=1:nelm
            % Calculate stiffness matrix
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            Ke = bar3ge(ec, ed, ep, es);                                  % Element stiffness matrix
            Ks = bar1e(ep, alpha(el));                                    % Spring stiffness matrix
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            % Assemble stiffness matrix
            K(indx_el, indx_el) = K(indx_el, indx_el) + Ke;
            K(indx_s, indx_s) = K(indx_s, indx_s) + Ks;
        end
        
        % Pseudo displacement increments
        da_r = solveq(K, -res, bc);
        da_p = solveq(K, p, bc);
        da_tot = a0 - a;                    % eq. 2.127, current da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, current dlambda_tot
        
        % Load parameter
        if iter == 1            % The first iteration we need to calculate dlambda differently due to denominator in c becoming zero
            if n > 1
                da_n = a - a_old;     % *
                s = sign(da_n'*da_p);
                dlambda = s*(l/sqrt(da_p'*da_p + psi*(p')*p));
            else
                dlambda = l/sqrt(da_p'*da_p + psi*(p')*p);     % s = 1 i första iterationen
            end
        else
            % Load parameter
            dlambda = -(c + 2*da_tot'*da_r)/(2*da_tot'*da_p + 2*psi*(p')*p*dlambda_tot);
        end
        
        lambda0 = lambda0 + dlambda;
        
        % Displacement vector
        da = da_r + dlambda*da_p;
        a0 = a0 + da;               % New displacement vector
        
        % Stresses and strains
        for el=1:nelm
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            
            % Internal element force vectors
            fe = bar3gf(ec, es, ed);
            fe_s = bar1f(alpha(el), eds);
            
            % Assemble force vector
            f_int(indx_el) = f_int(indx_el) + fe;
            f_int(indx_s) = f_int(indx_s) + fe_s;
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*p;
        c = da_tot'*da_tot + psi*dlambda_tot^2*(p')*p - l^2;
        %c = -(2*da_tot'*da_r + 2*da_tot'*da_p*dlambda + 2*psi*(p')*p*dlambda_tot*dlambda);
        
        res(bc(:,1)) = 0;
        % New states
        iter = 2;
    end
    
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    P(n) = -lambda*p(5);
    ybar(n) = y0B + a(5);               % a(5) is the displacement in y-dir on node 2
end

p2 = plot(ybar, P, 'b');              % Plotting force-displacement for
hold off
% xlabel('y')
% ylabel('P')
legend([p1 p2], 'Ideal system', 'Tilted system')

%% Exercise 2.4 b) Fulfilling constraint equation, ideal system

clear;
clc;

path_steps = 300;

fi = 65;
x0A = -cosd(fi);
x0B = 0;
x0C = cosd(fi);         %*
y0B = sind(fi);
z0B = 0;
alpha1 = 0.55*y0B;
alpha2 = 0;
alpha = [alpha1 alpha2];
TOL = 1E-5;
cTOL = 1E-5;

P = zeros(1, path_steps);
ybar = zeros(1, path_steps);

x = zeros(3, 1);
y = zeros(3, 1);
z = zeros(3, 1);
x(1) = x0A;
%x(2) = x0B;
y(2) = y0B;
%z(2) = z0B;
x(3) = x0C;
coord = [x y z];

ndof = 11;                  % No zero entries

f = zeros(ndof, 1);
f_int = zeros(ndof, 1);
da = zeros(ndof, 1);
a = zeros(ndof, 1);
K = zeros(ndof, ndof);

Edof=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9];          % Dof for each element
Edof_spring = [1 6 10];
%2 5 11];        % Dof for each spring
Enod = [1 2;
        2 3];
nelm = 2;
ep = [1 1];
bc = [1 0;
    2 0;
    3 0;
    7 0;
    8 0;
    9 0;
    10 0;
    11 0];

da_r = zeros(ndof, 1);
da_tot = zeros(ndof, 1);
da_p = zeros(ndof, 1);
lambda = 0;
dlambda = 0;
dlambda_tot = 0;
c = 0;
psi = 1;
p = zeros(ndof, 1);
p(5) = -0.3;
res = zeros(ndof, 1);
l = 0.018; % radius on sphere

for n = 1:path_steps
    disp(['Path Step: ', num2str(abs(n))]);
    iter = 1;
    
    a0 = a;
    lambda0 = lambda;
    
    dlambda_tot = 0;
    res = zeros(ndof, 1);
    da_tot = zeros(ndof, 1);
    
    while norm(res) > TOL && abs(c) > cTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);
        
        K = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el=1:nelm
            % Calculate stiffness matrix
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates for each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            Ke = bar3ge(ec, ed, ep, es);                                  % Element stiffness matrix
            Ks = bar1e(ep, alpha(el));                                    % Spring stiffness matrix
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            % Assemble stiffness matrix
            K(indx_el, indx_el) = K(indx_el, indx_el) + Ke;
            K(indx_s, indx_s) = K(indx_s, indx_s) + Ks;
        end
        
        % Pseudo displacement increments
        da_r = solveq(K, -res, bc);
        da_p = solveq(K, p, bc);
        da_tot = a0 - a;                    % eq. 2.127, current da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, current dlambda_tot
        
        % Load parameter
        if iter == 1            % The first iteration we need to calculate dlambda differently due to denominator in c becoming zero
            if n > 1
                da_n = a - a_old;     % ???
                s = sign(da_n'*da_p);
                dlambda = s*(l/sqrt(da_p'*da_p + psi*(p')*p));
            else
                dlambda = l/sqrt(da_p'*da_p + psi*(p')*p);     % s = 1 i första iterationen
            end
        else
            % Load parameter
            a1 = da_p'*da_p + psi*(p')*p;
            a2 = 2*da_p'*(da_tot+da_r) + 2*psi*dlambda_tot*(p')*p;
            a3 = (da_tot+da_r)'*(da_tot+da_r)+psi*dlambda_tot^2*(p')*p-l^2;
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
                break;
            end
            
        end
        
        lambda0 = lambda0 + dlambda;
        
        % Displacement vector
        da = da_r + dlambda*da_p;
        a0 = a0 + da;               % New displacement vector
        
        % Stresses and strains
        for el=1:nelm
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            
            % Internal element force vectors
            fe = bar3gf(ec, es, ed);
            fe_s = bar1f(alpha(el), eds);
            
            % Assemble force vector
            f_int(indx_el) = f_int(indx_el) + fe;
            f_int(indx_s) = f_int(indx_s) + fe_s;
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*p;
        c = da_tot'*da_tot + psi*dlambda_tot^2*(p')*p - l^2;
        %c = -(2*da_tot'*da_r + 2*da_tot'*da_p*dlambda + 2*psi*(p')*p*dlambda_tot*dlambda);
        
        res(bc(:,1)) = 0;
        % New states
        iter = 2;
    end
    
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    P(n) = -lambda*p(5);
    ybar(n) = y0B + a(5);               % a(5) is the displacement in y-dir on node 2
end

figure(2)
p1 = plot(ybar, P, 'r');              % Plotting force-displacement for
hold on
xlabel('y')
ylabel('P')
xlim([-2 1]);
ylim([-0.5 3]);
title('Constraint equation, force-displacement curve')
%legend('Ideal system')

%% Exercise 2.4 b) Fulfilling constraint equation, tilted system

path_steps = 300;

fi = 65;
x0A = -cosd(fi);        % - ???
x0B = 0;
x0C = cosd(fi);         %???
y0B = sind(fi);
z0B = 0.1;
alpha1 = 0.55*y0B;
alpha2 = 0;
alpha = [alpha1 alpha2];
TOL = 1E-5;
cTOL = 1E-5;

P = zeros(1, path_steps);
ybar = zeros(1, path_steps);

x = zeros(3, 1);
y = zeros(3, 1);
z = zeros(3, 1);
x(1) = x0A;
%x(2) = x0B;
y(2) = y0B;
z(2) = z0B;
x(3) = x0C;
coord = [x y z];

ndof = 11;                  % No zero entries

f = zeros(ndof, 1);
f_int = zeros(ndof, 1);
da = zeros(ndof, 1);
a = zeros(ndof, 1);
K = zeros(ndof, ndof);

Edof=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9];          % Dof for each element
Edof_spring = [1 6 10];
%2 5 11];        % Dof for each spring
Enod = [1 2;
    2 3];
nelm = 2;
ep = [1 1];
bc = [1 0;
    2 0;
    3 0;
    7 0;
    8 0;
    9 0;
    10 0;
    11 0];

da_r = zeros(ndof, 1);
da_tot = zeros(ndof, 1);        % scalar instead of vector?
da_p = zeros(ndof, 1);
lambda = 0;
dlambda = 0;
dlambda_tot = 0;
c = 0;
psi = 1;
p = zeros(ndof, 1);
p(5) = -0.3;
res = zeros(ndof, 1);
l = 0.018; % radius on sphere

for n = 1:path_steps
    disp(['Path Step: ', num2str(abs(n))]);
    iter = 1;
    
    a0 = a;
    lambda0 = lambda;
    
    dlambda_tot = 0;
    res = zeros(ndof, 1);
    da_tot = zeros(ndof, 1);
    
    while norm(res) > TOL && abs(c) > cTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);
        
        K = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el=1:nelm
            % Calculate stiffness matrix
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            Ke = bar3ge(ec, ed, ep, es);                                  % Element stiffness matrix
            Ks = bar1e(ep, alpha(el));                                    % Spring stiffness matrix
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            % Assemble stiffness matrix
            K(indx_el, indx_el) = K(indx_el, indx_el) + Ke;
            K(indx_s, indx_s) = K(indx_s, indx_s) + Ks;
        end
        
        % Pseudo displacement increments
        da_r = solveq(K, -res, bc);
        da_p = solveq(K, p, bc);
        da_tot = a0 - a;                    % eq. 2.127, current da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, current dlambda_tot
        
        % Load parameter
        if iter == 1            % The first iteration we need to calculate dlambda differently due to denominator in c becoming zero
            if n > 1
                da_n = a - a_old;     % ???
                s = sign(da_n'*da_p);
                dlambda = s*(l/sqrt(da_p'*da_p + psi*(p')*p));
            else
                dlambda = l/sqrt(da_p'*da_p + psi*(p')*p);     % s = 1 i första iterationen
            end
        else
            % Load parameter
            a1 = da_p'*da_p + psi*(p')*p;
            a2 = 2*da_p'*(da_tot+da_r) + 2*psi*dlambda_tot*(p')*p;
            a3 = (da_tot+da_r)'*(da_tot+da_r)+psi*dlambda_tot^2*(p')*p-l^2;
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
                break;
            end
            
        end
        
        lambda0 = lambda0 + dlambda;
        
        % Displacement vector
        da = da_r + dlambda*da_p;
        a0 = a0 + da;               % New displacement vector
        
        % Stresses and strains
        for el=1:nelm
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(1,:), a0);                          % Spring displacement
            [es, eg] = bar3gs(ec, ed, ep);                                % Normal force es and Green's strain eg
            
            indx_el = Edof(el, 2:end);                                    % Dof for element
            indx_s = Edof_spring(1, 2:end);                               % Dof for spring
            
            % Internal element force vectors
            fe = bar3gf(ec, es, ed);
            fe_s = bar1f(alpha(el), eds);
            
            % Assemble force vector
            f_int(indx_el) = f_int(indx_el) + fe;
            f_int(indx_s) = f_int(indx_s) + fe_s;
        end
        
        da_tot = a0 - a;                    % eq. 2.127, new da_tot
        dlambda_tot = lambda0 - lambda;     % eq. 2.127, new dlambda_tot
        
        % Out of balance forces
        res = f_int - lambda0*p;
        c = da_tot'*da_tot + psi*dlambda_tot^2*(p')*p - l^2;
        %c = -(2*da_tot'*da_r + 2*da_tot'*da_p*dlambda + 2*psi*(p')*p*dlambda_tot*dlambda);
        
        res(bc(:,1)) = 0;
        % New states
        iter = 2;
    end
    
    a_old = a;
    a = a0;
    lambda = lambda0;
    
    P(n) = -lambda*p(5);
    ybar(n) = y0B + a(5);               % a(5) is the displacement in y-dir on node 2
end

% figure(2)
p2 = plot(ybar, P, 'b');              % Plotting force-displacement for
hold on
% xlabel('y')
% ylabel('P')
legend([p1 p2], 'Ideal system', 'Tilted system')