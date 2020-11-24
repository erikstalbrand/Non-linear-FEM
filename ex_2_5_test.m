%% Exercise 2.5, Fully linearized
clear;
clc;

path_steps = 300;

fi = 65;
x0A = -cosd(fi);        % *
x0B = 0;
x0C = cosd(fi);         % *
y0B = sind(fi);
z0B = 0;
alpha1 = 0.55*y0B;
alpha2 = 1.5*y0B;
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

p = zeros(ndof, 1);
df = zeros(ndof, 1);
p(5) = -0.01;
f_int = zeros(ndof, 1);
da = zeros(ndof, 1);
a = zeros(ndof, 1);
K = zeros(ndof, ndof);

Edof=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9];          % Dof for each element
Edof_spring = [1 6 10;
    2 5 11];        % Dof for each spring
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
res = zeros(ndof, 1);
l = 0.02; % radius on sphere

for n = 1:path_steps
    disp(['Path Step: ', num2str(abs(n))]);
    iter = 1;
    
    a0 = a;
    lambda0 = lambda;
    
    dlambda_tot = 0;
    res = zeros(ndof, 1);
    da_tot = zeros(ndof, 1);
    
    %f = f + df;
    %res = f_int - lambda0*f;
    %res(:, 1) = 0;
    
    while norm(res) > TOL && abs(c) > cTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);
        
        K = zeros(ndof, ndof);
        f_int = zeros(ndof, 1);
        
        for el=1:nelm
            % Calculate stiffness matrix
            ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';        % Coordinates fo each element
            ed = extract(Edof(el,:), a0);                                 % Element displacement
            eds = extract(Edof_spring(el,:), a0);                         % Spring displacement
            
            lambdaNL = stretch1D(ec, ed);          % stretch
            sg = stress1D(ep, lambdaNL);           % stress
            Et = dmat1D(ep, lambdaNL);             % Material tangent stiffness
            
            [es, eg] = bar3gsNL(ec, ed, ep, sg);                          % Normal force es and Green's strain eg
            
            Ke = bar3geNL(ec, ed, ep, es, Et);                            % Element stiffness matrix
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
            
            lambdaNL = stretch1D(ec, ed);          % stretch
            sg = stress1D(ep, lambdaNL);           % stress
            Et = dmat1D(ep, lambdaNL);             % Material tangent stiffness
            [es, eg] = bar3gsNL(ec, ed, ep, sg);                                % Normal force es and Green's strain eg
            
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
    
    %P(n) = abs(f(5));
    P(n) = -lambda*p(5);
    ybar(n) = y0B + a(5);               % a(5) is the displacement in y-dir on node 2
end

figure(1)
p1 = plot(ybar, P, 'b', 'Linewidth', 1.5);              % Plotting force-displacement for
hold on
xlabel('y')
ylabel('P')
title('Force-displacement curve')
%legend('Ideal system')


%% Tilted system

load_steps = 300;

fi = 65;
x0A = -cosd(fi);        % - ???
x0B = 0;
x0C = cosd(fi);         %???
y0B = sind(fi);
z0B = 0.1;
alpha1 = 0.55*y0B;
alpha2 = 1.5*y0B;
alpha = [alpha1 alpha2];    %Spring stiffness
TOL = 1E-5;
daTOL = 1E-5;

P = zeros(load_steps, 1);
ybar = zeros(load_steps, 1);

%----------------------------------------------------%

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
df = zeros(ndof, 1);
df(5) = -0.01;

f_int = zeros(ndof, 1);
da = zeros(ndof, 1);
a = zeros(ndof, 1);
K = zeros(ndof, ndof);
res = zeros(ndof, 1);

Edof=[1 1 2 3 4 5 6; 
       2 4 5 6 7 8 9];          % Dof for each element
Edof_spring = [1 6 10; 
                2 5 11];        % Dof for each spring
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

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    f = f + df;
    iter = 1;
    res = f_int - f;
    res(:, 1) = 0;
    
    
    while norm(res) > TOL && norm(da) > daTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);

        K = K*0;
        f_int = f_int*0;
        
        for el = 1:nelm
           ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';        % Coordinates fo each element
           ed = extract(Edof(el,:), a);                                 % Element displacement
           eds = extract(Edof_spring(el,:), a);                         % Spring displacement
           
           lambda = stretch1D(ec, ed);          % stretch
           sg = stress1D(ep, lambda);           % stress
           Et = dmat1D(ep, lambda);             % Material tangent stiffness
           
           [es, eg] = bar3gsNL(ec, ed, ep, sg);                               % Normal force es and Green's strain eg
           
           Ke = bar3geNL(ec, ed, ep, es, Et);                                 % Element stiffness matrix
           Ks = bar1e(ep, alpha(el));                                         % Spring stiffness matrix
           
           indx_el = Edof(el, 2:end);                   % Dof for element
           indx_s = Edof_spring(el, 2:end);             % Dof for spring
           
           % Assemble stiffness matrix
           K(indx_el, indx_el) = K(indx_el, indx_el) + Ke;
           K(indx_s, indx_s) = K(indx_s, indx_s) + Ks;
           
           % Internal element force vectors
           fe = bar3gf(ec, es, ed);
           fe_s = bar1f(alpha(el), eds);
           
           % Assemble force vector
           f_int(indx_el) = f_int(indx_el) + fe;
           f_int(indx_s) = f_int(indx_s) + fe_s;

%            f_int = insert(Edof(el, 2:end), f_int, fe);
%            f_int = insert(Edof(el, 2:end), f_int, fe_s);
        end
        
        % New states
        res = f_int - f;
        res(bc(:,1)) = 0;
        da = solveq(K, -res, bc);
        a = a + da;
        iter = 2;
    end
    
% Numerical force results from Newton Raphson
P(n) = abs(f(5));                                     % f(5) is the force applied in y-dir on node 2
ybar(n) = y0B + a(5);
end

% figure(1)
p2 = plot(ybar, P, 'r', 'Linewidth', 1.5);              % Plotting force-displacement for
hold on
% xlabel('y')
% ylabel('P')
legend([p1 p2], 'Ideal system', 'Tilted system')