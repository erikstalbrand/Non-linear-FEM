function [Me] = cont2D3m(ec, t, rho)

C = [1         1       1;
    ec(1,1) ec(1,2) ec(1,3);
    ec(2,1) ec(2,2) ec(2,3)];

A0e = 1/2*det(C);

Me = rho*A0e*t/12*[2 0 1 0 1 0;
                   0 2 0 1 0 1;
                   1 0 2 0 1 0;
                   0 1 0 2 0 1;
                   1 0 1 0 2 0;
                   0 1 0 1 0 2];
end

