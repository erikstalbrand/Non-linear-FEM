%% Exercise 2.2 a)

clc;
clear all;

load_steps = 300;

fi = 65;
x0A = -cosd(fi);        % - ???
x0B = 0;
x0C = cosd(fi);         %???
y0B = sind(fi);
z0B = 0;
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
%z(2) = z0B;
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
    
    while norm(res) > TOL || norm(da) > daTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);

        K = K*0;
        f_int = f_int*0;
        
        for el = 1:nelm
           ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
           ed = extract(Edof(el,:), a);                                 % Element displacement
           eds = extract(Edof_spring(el,:), a);                  % Spring displacement
           [es, eg] = bar3gs(ec, ed, ep);                               % Normal force es and Green's strain eg
           
           Ke = bar3ge(ec, ed, ep, es);                                 % Element stiffness matrix
           Ks = bar1e(ep, alpha(el));                                   % Spring stiffness matrix
           
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
P(n) = abs(f(5));                            % f(5) is the force applied in y-dir on node 2
ybar(n) = y0B + a(5);
end


figure(1)
p1 = plot(ybar, P, 'g', 'Linewidth', 1);              % Plotting force-displacement for
hold on
xlabel('y')
ylabel('P')
%legend('Ideal system')

%% b)

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
    
    while norm(res) > TOL || norm(da) > daTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);

        K = K*0;
        f_int = f_int*0;
        
        for el = 1:nelm
           ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
           ed = extract(Edof(el,:), a);                                 % Element displacement
           eds = extract(Edof_spring(el,:), a);                  % Spring displacement
           [es, eg] = bar3gs(ec, ed, ep);                               % Normal force es and Green's strain eg
           
           Ke = bar3ge(ec, ed, ep, es);                                 % Element stiffness matrix
           Ks = bar1e(ep, alpha(el));                                   % Spring stiffness matrix
           
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
P(n) = abs(f(5));                            % f(5) is the force applied in y-dir on node 2
ybar(n) = y0B + a(5);
end

%figure(2)
p2 = plot(ybar, P, 'r', 'Linewidth', 1);              % Plotting force-displacement for
hold on
xlabel('y')
ylabel('P')
%legend('Tilted frame structure')

%% c)

load_steps = 300;

fi = 65;
x0A = -cosd(fi);        % - ???
x0B = 0.1;
x0C = cosd(fi);         %???
y0B = sind(fi);
z0B = 0;
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
x(2) = x0B;
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
    
    while norm(res) > TOL || norm(da) > daTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);

        K = K*0;
        f_int = f_int*0;
        
        for el = 1:nelm
           ec = [coord(Enod(el, 1), :); coord(Enod(el, 2), :)]';         % Coordinates fo each element
           ed = extract(Edof(el,:), a);                                 % Element displacement
           eds = extract(Edof_spring(el,:), a);                  % Spring displacement
           [es, eg] = bar3gs(ec, ed, ep);                               % Normal force es and Green's strain eg
           
           Ke = bar3ge(ec, ed, ep, es);                                 % Element stiffness matrix
           Ks = bar1e(ep, alpha(el));                                   % Spring stiffness matrix
           
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
P(n) = abs(f(5));                            % f(5) is the force applied in y-dir on node 2
ybar(n) = y0B + a(5);
end


%figure(3)
p3 = plot(ybar, P, 'b', 'Linewidth', 1);              % Plotting force-displacement for
hold on
xlabel('y')
ylabel('P')
%legend('Different bar lengths')
legend([p1 p2 p3], 'Ideal system', 'Tilted frame structure', 'Different bar lengths')
