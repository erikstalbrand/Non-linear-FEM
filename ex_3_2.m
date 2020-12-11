%% Exercise 3.2, St. Venant-Kirchhoff
clear;
clc;
close all;

load_steps = 110;
res = zeros(3, 1);
dS = 0.001;
S = zeros(3, 1);
M = zeros(3, 1);
C = zeros(3);
dC = zeros(3);

v = 0.3;
E = 1;
mpara = [E, v];
e1 = 0;
e2 = 0;
TOL = 1E-5;

S_plot = zeros(load_steps, 1);
e1_plot = zeros(load_steps, 1);
e2_plot = zeros(load_steps, 1);
sigma_plot = zeros(load_steps, 1);
defgrad = [1+e1 0 0 1+e2]';

F_init = [defgrad(1) defgrad(2) 0;
          defgrad(3) defgrad(4) 0;
             0          0       1];
    
C = F_init'*F_init;


for n = 1:load_steps
    
    S(1, 1) = S(1, 1) + dS;     % Spänningen vi vill ha
    
    M = stressMater2D1(1, mpara, defgrad);      % Spänningen vi har
    res = S - M;                                % Felet som uppstår
    %res = res*0;
    
    iter = 1;
    while norm(res) > TOL || iter == 1
        
        D = dMater2D1(1, mpara, defgrad);         % Tangential material stiffness
        
        Kt = -D/2;
        dC = -Kt\res;
        C(1, 1) = C(1, 1) + dC(1);         % Dependence on deformation instead of strain when we have large deformation
        C(2, 2) = C(2, 2) + dC(2);
        C(3, 3) = C(3, 3) + dC(3);
        
        F = sqrt(C);        % Deformation gradient tensor
        defgrad = [F(1, 1) F(1, 2) F(2, 1) F(2, 2)]';
        
        M = stressMater2D1(1, mpara, defgrad);      % Ny spänning
        
        %---------------------------------
        sigma = stressMater2D1(2, mpara, defgrad);   % FEL I LÖSNIINGEN? DE HAR tau OCH INTE sigma, dvs. flag == 2 istället för flag == 3
        %---------------------------------
        res = S - M;
        iter = 2;
    end
    e1 = defgrad(1) - 1;
    e2 = defgrad(4) - 1;
    
    S_plot(n) = S(1, 1);
    e1_plot(n) = e1;
    e2_plot(n) = e2;
    sigma_plot(n) = sigma(1);
    
end

figure(1)
p1 = plot(e1_plot, S_plot);
xlabel('\epsilon_1');
ylabel('S_{11} and \sigma_{11}')
title('St. Venant-Kirchhoff Model')
hold on
grid on
p2 = plot(e1_plot, sigma_plot);
legend([p1 p2], 'S_{11}', '\sigma_{11}')

figure(2)
p3 = plot(e2_plot, S_plot);
xlabel('\epsilon_2');
ylabel('S_{11} and \sigma_{11}')
title('St. Venant-Kirchhoff Model')
hold on
grid on
p4 = plot(e2_plot, sigma_plot);
legend([p3 p4], 'S_{11}', '\sigma_{11}')

%% Exercise 3.2, Neo-Hookean Model

load_steps = 1000;
res = zeros(3, 1);
dS = 0.0001;
S = zeros(3, 1);
M = zeros(3, 1);
C = zeros(3);
dC = zeros(3);

v = 0.3;
E = 1;
mpara = [E, v];
e1 = 0;
e2 = 0;
TOL = 1E-5;

S_plot = zeros(load_steps, 1);
e1_plot = zeros(load_steps, 1);
e2_plot = zeros(load_steps, 1);
sigma_plot = zeros(load_steps, 1);
defgrad = [1+e1 0 0 1+e2]';
    F_init = [defgrad(1) defgrad(2) 0;
        defgrad(3) defgrad(4) 0;
        0          0       1];
    
    
    C = F_init'*F_init;


for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    S(1) = S(1) + dS;                                   % Spänningen vi vill ha
    
    M = stressMater2D2(1, mpara, defgrad);              % Spänningen vi har
    res = S - M;                                        % Felet som uppstår
        
    iter = 1;
    
    while norm(res) > TOL || iter == 1
        disp(['Residual: ', num2str(abs(norm(res)))]);
        
        D = dMater2D2(1, mpara, defgrad);               % Tangential material stiffness
        
        Kt = -D/2;
        dC = -Kt\res;
        
        C(1, 1) = C(1, 1) + dC(1);                      % Dependence on deformation instead of strain when we have large deformation
        C(2, 2) = C(2, 2) + dC(2);
        C(3, 3) = C(3, 3) + dC(3);
        
        F = sqrt(C);                                    % Deformation gradient tensor, diagonal matrix
        defgrad = [F(1, 1) F(1, 2) F(2, 1) F(2, 2)]';
        
        M = stressMater2D2(1, mpara, defgrad);          % Ny spänning
        
        %---------------------------------
        sigma = stressMater2D2(2,mpara,defgrad);        % FEL I LÖSNIINGEN? DE HAR tau OCH INTE sigma, dvs. flag == 2 istället för flag == 3
        
        %---------------------------------
        res = S - M;
        iter = 2;
    end
    e1 = defgrad(1) - 1;
    e2 = defgrad(4) - 1;
    
    S_plot(n) = S(1);
    e1_plot(n) = e1;
    e2_plot(n) = e2;
    sigma_plot(n) = sigma(1);
    
end

figure(3)
p5 = plot(e1_plot, S_plot);
xlabel('\epsilon_1');
ylabel('S_{11} and \sigma_{11}')
title('Neo-Hookean Model')
hold on
grid on
p6 = plot(e1_plot, sigma_plot);
legend([p5 p6], 'S_{11}', '\sigma_{11}')


figure(4)
p7 = plot(e2_plot, S_plot);
xlabel('\epsilon_2');
ylabel('S_{11} and \sigma_{11}')
title('Neo-Hookean Model')
hold on
grid on
p8 = plot(e2_plot, sigma_plot);
legend([p7 p8], 'S_{11}', '\sigma_{11}')

