%% Exercise 2.1 a) load controlled Newton-Raphson scheme for ideal system
clear;
clc;

P = 0;
load_steps = 50;
df = [0; -0.1; 0];
f_int = zeros(3, 1);
f = [0; P; 0];
%f = zeros(3, 1);
x = 0;
y = 0;
z = 0;
fi = 65;
x0A = cosd(fi);
x0B = 0;
y0B = sind(fi);
z0B = 0;
alpha1 = 0.55*y0B;
alpha2 = 1.5*y0B;
TOL = 1E-4;
duTOL = 1E-4;
da = zeros(3, 1);
a = zeros(3, 1);

P = zeros(3, load_steps);
ybar = zeros(3, load_steps);

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    f = f + df;
    iter = 1;
    res = f_int - f;
    
    while abs(norm(res)) > TOL && abs(norm(da)) > duTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);

        k11 = 3*x^2 + y^2 + z^2 + 2*x0A^2 - y0B^2 - z0B^2;
        k22 = x^2 + 3*y^2 + z^2 + alpha2 - y0B^2 - z0B^2;
        k33 = x^2 + y^2 + 3*z^2 + alpha1 - y0B^2 - z0B^2;
        
        Kt = [k11 2*x*y 2*x*z;
              2*y*x k22 2*y*z;
              2*z*x 2*z*y k33];
        
        da = -Kt\res;        
        a = a + da;
        
        x = x0B + a(1);
        y = y0B + a(2);
        z = z0B + a(3);
        
        f_int = [x*(x^2+y^2+z^2 + 2*x0A^2 - y0B^2 - z0B^2);
                 y*(x^2+y^2+z^2 - y0B^2 - z0B^2) + alpha2*(y-y0B);
                 z*(x^2+y^2+z^2 - y0B^2 - z0B^2) + alpha1*(z-z0B)];
        
        % New states
        res = f_int - f;
        iter = 2;

    end
    
    P(:, n) = f;                   % Numerical force results from Newton Raphson
    ybar(:, n) = [x; 
                  y; 
                  z];             % Normalised displacement
end

figure(1)
p1 = plot(ybar(2, :), -P(2, :), 'g');              % Plotting force-displacement for
hold on
xlabel('y')
ylabel('P')
% legend('Ideal system')

%% Exercise 2.1 b) load controlled Newton-Raphson scheme for tilted system

P = 0;
load_steps = 50;
df = [0; -0.1; 0];
f_int = zeros(3, 1);
f = [0; P; 0];
%f = zeros(3, 1);
x = 0;
y = 0;
z = 0;
fi = 65;
x0A = cosd(fi);
x0B = 0;
y0B = sind(fi);
z0B = 0.1;
alpha1 = 0.55*y0B;
alpha2 = 1.5*y0B;
TOL = 1E-10;
duTOL = 1E-10;
da = 0;
a = 0;

P = zeros(3, load_steps);
ybar = zeros(3, load_steps);

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    f = f + df;
    iter = 1;
    res = f_int - f;
    
    while abs(norm(res)) > TOL && abs(norm(da)) > duTOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);

        k11 = 3*x^2 + y^2 + z^2 + 2*x0A^2 - y0B^2 - z0B^2;
        k22 = x^2 + 3*y^2 + z^2 + alpha2 - y0B^2 - z0B^2;
        k33 = x^2 + y^2 + 3*z^2 + alpha1 - y0B^2 - z0B^2;
        
        Kt = [k11 2*x*y 2*x*z;
              2*y*x k22 2*y*z;
              2*z*x 2*z*y k33];
        
        da = -Kt\res;           % ???
        %du = -res./Kt;
        
        a = a + da;
        
        x = x0B + a(1);
        y = y0B + a(2);
        z = z0B + a(3);
        
        f_int = [x*(x^2+y^2+z^2 + 2*x0A^2 - y0B^2 - z0B^2);
                 y*(x^2+y^2+z^2- y0B^2 - z0B^2) + alpha2*(y-y0B);
                 z*(x^2+y^2+z^2- y0B^2 - z0B^2) + alpha1*(z-z0B)];
        
        % New states
        res = f_int - f;
        iter = 2;

    end
    
    P(:, n) = f;                  % Numerical force results from Newton Raphson
    ybar(:, n) = [x; 
                  y; 
                  z];             % Normalised displacement
end

p2 = plot(ybar(2, :), -P(2, :), 'b');              % Plotting force-displacement for
%legend('Tilted system')
legend([p1 p2], 'Ideal system', 'Tilted system')
hold off
