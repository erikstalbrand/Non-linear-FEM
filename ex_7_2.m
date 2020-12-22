%% Exercise 7.2 Newmark algorithm

clear;
clc;
close all;

load('geom_E72.mat');

load_steps = 200;
TOL = 1E-3;
nmax = load_steps;

rho = 7800e-9;
beta = 1/4;
gamma = 2*beta;

nrelem = length(edof(:,1));
nrdof = max(max(edof));
flag = 1;                         % Cauchy stress
Kt = zeros(nrdof, nrdof);
f_int = zeros(nrdof, 1);
f = zeros(nrdof, 1);

force_plot = zeros(load_steps, 1);
disp_plot = zeros(load_steps, 1);
time_plot = zeros(load_steps, 1);

d1 = 0;
M = zeros(nrdof, nrdof);
C = d1*M;
t = 1;
E = 210E3;
vp = 0.3;
mpara = [E, vp];

dt = 1e-04;   % Time-step
tmax = dt*load_steps;    % Max time

a = zeros(nrdof, 1);
v = zeros(nrdof, 1);
acc = zeros(nrdof, 1);
load = zeros(nrdof,1);
res_eff = zeros(nrdof, 1);

load(find(P)) = 20000000;    % Loading rate, dt*load = df, df_tot = 10 kN
count = 0;

a_tot = zeros(nrdof, load_steps);

nelm = nrelem;

% Mass matrix
for el = 1:nelm
    ec = [ex(el, :); ey(el, :)];
    Me = cont2D3m(ec, t, rho);
    indx = edof(el, 2:end);
    M(indx, indx) = M(indx, indx) + Me;
end

% Force controlled Newton-Raphson
for n = 0:dt:tmax
    disp(['LOAD STEP:', num2str(count)]);
    
    % New load level    
    if count <= 5
    f = n*load;       % f_n+1, n=dt*50=0.005
    else if count == 10
        f = f*0;
        end
    end
    count = count + 1;
%     force_plot(count) = -f(96);
    
    % Predictor step with same acceleration, 7.19
    a = a + dt*v + (dt^2)*(1/2)*acc;
    v = v + dt*acc;

    res_eff = f_int - f + M*acc;
    iter = 1;
    
    while norm(res_eff) > TOL || iter == 1
%         disp(['Residual: ', num2str(norm(res))])
%         disp(['Effective Residual: ', num2str(norm(res_eff))])
        
        % Reset stiffness matrix and force vector
        Kt = zeros(nrdof, nrdof);
        f_int = zeros(nrdof, 1);
        
        for el = 1:nelm
            ed = extract(edof(el, :), a);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed);
            D = dMater2D2(flag, mpara, defgrad);
            stress = stressMater2D2(flag, mpara, defgrad);
            Ke = cont2D3te(ec, t, D, ed, stress);
            
            indx = edof(el, 2:end);
            Kt(indx, indx) = Kt(indx, indx) + Ke;
        end
        
        K_eff = Kt + (1/(beta*dt^2))*M;            % 7.16

        da = solveq(K_eff, -res_eff, bc);
        
        %  Displacement, velocity and acceleration vectors
        a = a + da;
        v = v + ((gamma)/(beta*dt))*da;
        acc = acc + (1/(beta*dt^2))*da;

        ed = extract(edof, a);
        % Stresses and strains
        for el = 1:nelm
%             ed = extract(edof(el, :), a);
            ec = [ex(el, :); ey(el, :)];
            defgrad = cont2D3ts(ec, ed(el, :));
            stress = stressMater2D2(flag, mpara, defgrad);
            fe = cont2D3tf(ec, t, ed, stress);
            
            indx = edof(el, 2:end);
            f_int(indx) = f_int(indx) + fe';
        end
        
        % Out of balance forces

        res_eff = f_int + M*acc - f;
        res_eff(bc(:, 1)) = 0;
        disp(['Effective Residual: ', num2str(norm(res_eff))])
        
        iter = 2;
    end
    
    a_tot(:, count) = a;
    force_plot(count) = f(10);
    disp_plot(count) = a(10);
    time_plot(count) = n;

    ed = extract(edof, a);
    figure(1)
    grid on
    drawnow update
    clf(figure(1));
%     eldraw2(ex, ey, [1 2 1], edof(:,1))
    eldraw2(ex,ey,[1 2 1]);
    eldisp2(ex, ey, ed, [1 4 1], 10);
    title('df: 200N, dt: 0.0001 s');
end

figure(3)
plot(time_plot, force_plot)

%% Post-processing
% Plots parameters

fontsize = 15;

wsize = [300, 300, 900, 600];

h1 = figure(4);

a = get(gca,'XTickLabel');

set(gca,'XTickLabel',a,'FontName','Times','fontsize',fontsize);

set(gcf, 'Position',  wsize);

mag = 10;

fr = 0;


for load=1:nmax % nmax is the number of load steps

    fr = fr+1;

    ed(:,:) = extract(edof, a_tot(:, load)); %(all_ed(load,:,:); % all_ed is the element displacements saved for each load step, (nmax, nelm, 6)

    clf;

    grid on;

%     title(['load step: ', num2str(load),'      magnification: ',num2str(mag),'     Load rate: ',num2str(max(abs(dfext))/dt)])

    eldisp2(ex, ey, ed, [1 4 1], mag);

    hold on

    edmirr = ed; % Plot mirrored circle

    edmirr(:,1:2:end) = -edmirr(:,1:2:end);

    eldisp2(-ex, ey, edmirr, [1 4 1], mag);

    xlim([-70 70])

    ylim([0 250])

    drawnow;

    F(fr) = getframe(gcf);

end

% Save the movie

v = VideoWriter('ring.avi');

open(v)

writeVideo(v,F)

close(v)
