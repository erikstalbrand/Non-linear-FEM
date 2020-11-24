function [fe] = bar3gf(ec, es, ed)

x_0 = ec(:, 2)-ec(:, 1);
l0 = norm(x_0);

dx = [(ec(1,2)+ed(4)) - (ec(1,1)+ed(1));
     (ec(2,2)+ed(5)) - (ec(2,1)+ed(2));
     (ec(3,2)+ed(6)) - (ec(3,1)+ed(3))];

fe = (es/l0)*[-dx; dx];     %eq. 2.15

% u = ec + [ed(1:3); ed(4:6)]';
% x = u(:, 2) - u(:, 1);
% 
% 
% fB = (es/l0)*x;          %???
% fA = -fB;
% 
% fe = vertcat(fA, fB);
end