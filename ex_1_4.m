%% Exercise 1.4
clc;
clear
close all;

path_steps = 50;
y0 = sind(60);      % Initial angle
u = 0;              % Initial displacement
du = 0;             % Initial displacement
kbar = 0;           % No spring
force_curve = zeros(path_steps, 1);   % Normalised force
ubar = zeros(path_steps, 1);          % Normalised diplacement
TOL = 1E-10;
cTOL = 1E-10;

f_ext = 0;
%df = 0.1;
f_int = 0;
du_tot = 0;
dP_tot = 0;
P = 0;
p = 1;
fi = 1;
l = 0.1; % Radius of circle
dlambda = 0;
c = 0;
lambda = 0;

for n = 1:path_steps
    disp(['Path Step: ', num2str(abs(n))]);
    
    iter = 1;
    res = 0;
    c = 0;
    P0 = f_ext;
    u0 = u;
    
    while abs(res) > TOL || abs(c) > cTOL || iter == 1      % Equilibrium iteration
        disp(['Residual: ', num2str(abs(res))]);
        
        Kt = (2/5)*((-3*(u-y0)/((1-2*u*y0+u^2)^(5/2))-(2*u-2*y0))*(y0-u)-(1/(lambda^3)-lambda^2));
        
        % Pseudo displacement increment
        du_r = -res/Kt;
        du_p = p/Kt;
        
        if iter == 1            % The first iteration we need to calculate dlambda differently due to denominator in c will become zero
            if n > 1
                s = sign(du_tot*du_p+dP_tot);
                %s = sign(du*du_p + dlambda);       VARFÖR FUNKAR INTE
                %DENNA???
                dlambda = s*(l/sqrt(du_p^2+fi*p^2));
            else
                dlambda = l/sqrt(du_p^2+fi*p^2); % s = 1 i första iterationen
            end
        else
            % Load parameter
            dlambda = -(c+2*du_tot*du_r)/(2*du_tot*du_p+2*fi*dP_tot*p);
        end
        
        % Update load
        dP = dlambda*p;
        P = P + dP;
        dP_tot = P - P0;     % eq. 1.92???
        
        % Displacement vector
        du = du_r + dlambda*du_p;
        u = u + du;
        du_tot = u - u0;    % eq. 1.92???
        
        %du_tot = du_tot + du;          VARFÖR INTE DESSA???
        %dP_tot = dP_tot + dP;
        
        lambda = sqrt(1 - 2*u*y0 + u^2);
        f_int = (2/5)*(1/(lambda^3)-lambda^2)*(y0-u); % + kbar*u;
        f_ext = P;
        
        % Out of balance forces
        res = f_int - f_ext;
        c = du_tot^2 + fi*dP_tot^2 - l^2;
        
        iter = 2;
        
    end
    
    
    force_curve(n+1) = f_ext;     % Numerical force results from Euler Forward
    ubar(n+1) = u;            % Normalised displacement
    
end


Pbar = [];
p = 0;
u = 0:0.04:2;
for i = 1:length(u)  % True equilibrium path
    
    lambda = sqrt(1 - 2*u(i)*y0 + u(i)^2);
    p = (2/5)*(1/(lambda^3)-lambda^2)*(y0-u(i));
    
    Pbar = [Pbar; p];
    
end


figure(1)
p1 = plot(ubar, force_curve, 'LineWidth', 2);
hold on
p2 = plot(u, Pbar, '^', 'LineWidth', 1.5);
hold off
xlabel('Displacement')
ylabel('Force')
legend([p1 p2], 'Fully linearized', 'True path')
xlim([0 2])