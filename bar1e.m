function [Ks] = bar1e(ep, alpha)

E = ep(1);
A = ep(2);

l0 = 1;
k = alpha*E*A/l0;

Ks = [k -k;
      -k k];
end