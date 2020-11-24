function [fe_s] = bar1f(alpha, eds)

dx = eds(2) - eds(1);

fe_s = alpha*[-dx; dx];     % eq. 2.15


% E = ep(1);
% A = ep(2);
% 
% l0 = 1;
% k = alpha*E*A/l0;
% 
% Ks = [k -k;
%       -k k];

end