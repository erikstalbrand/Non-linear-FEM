function [lambda] = stretch1D(ec, ed)
x_0 = ec(:, 2)-ec(:, 1);
l0 = norm(x_0);

% dx_0 = [ec(1,2) - ec(1,1);
%         ec(2,2) - ec(2,1);
%         ec(3,2) - ec(3,1)];
% l0 = sqrt((dx_0')*dx_0);

dx = [(ec(1,2)+ed(4)) - (ec(1,1)+ed(1));
     (ec(2,2)+ed(5)) - (ec(2,1)+ed(2));
     (ec(3,2)+ed(6)) - (ec(3,1)+ed(3))];
l = norm(dx);

%sqrt((dx')*dx);
 
lambda = l/l0;
end

