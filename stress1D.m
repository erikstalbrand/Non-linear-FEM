function [sg] = stress1D(ep, lambda)

sg = (ep(1)/5)*(lambda^2-1/(lambda^3));
%sg = (E/5)*(lam^2-1/(lam^3));
end

