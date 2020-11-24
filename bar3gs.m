function [es, ee] = bar3gs(ec, ed, ep)
x_0 = ec(:, 2)-ec(:, 1);
l0 = norm(x_0);

u = ec + [ed(1:3)', ed(4:6)'];

x = u(:, 2) - u(:, 1);
l = norm(x);

E = ep(1);
A = ep(2);

ee = (l^2-l0^2)/(2*l0^2);   % Green's strain
es = E*A*ee;                % Normal force
end

