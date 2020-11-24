function [Ke] = bar3geNL(ec, ed, ep, es, Et)
x_0 = ec(:, 2)-ec(:, 1);
l0 = norm(x_0);


% dx = [(ec(1,2)+ed(4)) - (ec(1,1)+ed(1));
%      (ec(2,2)+ed(5)) - (ec(2,1)+ed(2));
%      (ec(3,2)+ed(6)) - (ec(3,1)+ed(3))];
du = ed(4:6)' - ed(1:3)';
dx = x_0 + du;

%E = ep(1);
A = ep(2);

I = eye(3);

Ke = (Et*A/(l0^3))*[dx*dx' -dx*dx'; -dx*dx' dx*dx'] + es/l0*[I -I; -I I];

end

