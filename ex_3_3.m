%% Exercise 3.3, St. Venant-Kirchhoff Model

clear;
clc;
close all;

load_steps = 100;
res = zeros(3, 1);
dtau = [0;
        0;
        0.001];
tau = zeros(3, 1);
N = zeros(3, 1);

% F = zeros(3);
dF = zeros(3);
Kt = zeros(3);

v = 0.3;
E = 1;
mpara = [E, v];

F12 = 0;
F21 = 0;
TOL = 1E-4;

tau_plot = zeros(load_steps, 1);
F_plot = zeros(load_steps, 1);
sigma_plot = zeros(load_steps, 1);

defgrad = [1, F12, F21, 1]';

F = [defgrad(1) defgrad(2) 0;
          defgrad(3) defgrad(4) 0;
             0          0       1];

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    tau = tau + dtau;                                           % Spänningen vi vill ha
    
    N = stressMater2D1(2, mpara, defgrad);                      % Spänningen vi har, samma som om vi får ut S och multiplicerar med F,F'...
             
    res = tau - N;                                              % Felet som uppstår

    iter = 1;

    while max(abs(res)) > TOL || iter == 1
        
        Kt = zeros(3);
        %------------------------------------------------
        defgrad = [F(1, 1)+1E-6, F(1, 2), F(1, 2), F(2, 2)]';           % Symmetric
        N = stressMater2D1(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        resi = tau - N;
        Kt(:, 1) = (res(:)-resi(:))/(1E-6);
        %------------
        defgrad = [F(1, 1), F(1, 2), F(1, 2), F(2, 2)+1E-6]';           % Symmetric
        N = stressMater2D1(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        resi = tau - N;
        Kt(:, 2) = (res(:)-resi(:))/(1E-6);
        %------------
        defgrad = [F(1, 1), F(1, 2)+1E-6, F(1, 2)+1E-6, F(2, 2)]';           % Symmetric
        N = stressMater2D1(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        resi = tau - N;
        Kt(:, 3) = (res(:)-resi(:))/(1E-6);
        %------------------------------------------------
        
        dF = Kt\res;
        F(1, 1) = F(1, 1) + dF(1);
        F(2, 2) = F(2, 2) + dF(2);
        F(1, 2) = F(1, 2) + dF(3);
        
        defgrad = [F(1, 1) F(1, 2) F(1, 2) F(2, 2)]';
        
        N = stressMater2D1(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        
        res = tau - N;
        disp(['Residual: ', num2str(abs(norm(res)))]);
        iter = 2;
    end
    sigma = stressMater2D1(3,mpara,defgrad);
    F12 = defgrad(2);
    
    tau_plot(n+1) = N(3);
    F_plot(n+1) = F12;
    sigma_plot(n+1) = sigma(3);
end

figure(1)
p1 = plot(F_plot, tau_plot);
xlabel('F_{12}');
ylabel('\tau_{12}/E and \sigma_{12}/E')
title('St. Venant-Kirchhoff Model')
hold on
grid on
p2 = plot(F_plot, sigma_plot);
legend([p1 p2], '\tau_{12}', '\sigma_{12}')

%% Exercise 3.3, Neo-Hookean Model

load_steps = 100;
res = zeros(3, 1);
dtau = [0;
        0;
        0.001];
tau = zeros(3, 1);
N = zeros(3, 1);

dF = zeros(3);
Kt = zeros(3);

v = 0.3;
E = 1;
mpara = [E, v];

F12 = 0;
F21 = 0;
TOL = 1E-4;

tau_plot = zeros(load_steps, 1);
F_plot = zeros(load_steps, 1);
sigma_plot = zeros(load_steps, 1);

defgrad = [1, F12, F21, 1]';

F = [defgrad(1) defgrad(2) 0;
          defgrad(3) defgrad(4) 0;
             0          0       1];

for n = 1:load_steps
    disp(['Load Step: ', num2str(abs(n))]);
    
    tau = tau + dtau;                                           % Spänningen vi vill ha
    
    N = stressMater2D2(2, mpara, defgrad);                      % Spänningen vi har, samma som om vi får ut S och multiplicerar med F,F'...
             
    res = tau - N;                                              % Felet som uppstår

    iter = 1;
    
    while max(abs(res)) > TOL || iter == 1
        
        Kt = zeros(3);
        %------------------------------------------------
        defgrad = [F(1, 1)+1E-6, F(1, 2), F(1, 2), F(2, 2)]';           % Symmetric
        N = stressMater2D2(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        resi = tau - N;
        Kt(:, 1) = (res(:)-resi(:))/(1E-6);
        %------------
        defgrad = [F(1, 1), F(1, 2), F(1, 2), F(2, 2)+1E-6]';           % Symmetric
        N = stressMater2D2(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        resi = tau - N;
        Kt(:, 2) = (res(:)-resi(:))/(1E-6);
        %------------
        defgrad = [F(1, 1), F(1, 2)+1E-6, F(1, 2)+1E-6, F(2, 2)]';           % Symmetric
        N = stressMater2D2(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        resi = tau - N;
        Kt(:, 3) = (res(:)-resi(:))/(1E-6);
        %------------------------------------------------
        
        dF = Kt\res;
        F(1, 1) = F(1, 1) + dF(1);
        F(2, 2) = F(2, 2) + dF(2);
        F(1, 2) = F(1, 2) + dF(3);
        
        defgrad = [F(1, 1) F(1, 2) F(1, 2) F(2, 2)]';
        
        N = stressMater2D2(2, mpara, defgrad);                  % Ny spänning, flag == 2 för tau
        
        res = tau - N;
        disp(['Residual: ', num2str(abs(norm(res)))]);
        iter = 2;
    end
    sigma = stressMater2D2(3,mpara,defgrad);
    F12 = defgrad(2);
    
    tau_plot(n+1) = N(3);
    F_plot(n+1) = F12;
    sigma_plot(n+1) = sigma(3);
end

figure(2)
p3 = plot(F_plot, tau_plot);
xlabel('F_{12}');
ylabel('\tau_{12}/E and \sigma_{12}/E')
title('Neo-Hookean Model')
hold on
grid on
p4 = plot(F_plot, sigma_plot);
legend([p3 p4], '\tau_{12}', '\sigma_{12}')

