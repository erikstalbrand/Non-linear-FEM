function [defgrad] = cont2D3ts(ec, ed)

C = [1           1        1;
      ec(1, 1) ec(1, 2) ec(1, 3);
      ec(2, 1) ec(2, 2) ec(2, 3)];
  
A0e = 1/2*det(C);
% v0 = A0e*t;

dN1dx0 = 1/(2*A0e)*(ec(2, 2) - ec(2, 3));
dN2dx0 = 1/(2*A0e)*(ec(2, 3) - ec(2, 1));
dN3dx0 = 1/(2*A0e)*(ec(2, 1) - ec(2, 2));
dN1dy0 = 1/(2*A0e)*(ec(1, 3) - ec(1, 2));
dN2dy0 = 1/(2*A0e)*(ec(1, 1) - ec(1, 3));
dN3dy0 = 1/(2*A0e)*(ec(1, 2) - ec(1, 1));

H0e = [dN1dx0 0 dN2dx0 0 dN3dx0 0;      % 4x6
      dN1dy0 0 dN2dy0 0 dN3dy0 0;
      0 dN1dx0 0 dN2dx0 0 dN3dx0;
      0 dN1dy0 0 dN2dy0 0 dN3dy0];
 
ae = H0e*ed';       % 4x1

% F = [1+ae(1) ae(2)  0;
%      ae(3)  1+ae(4) 0
%      0      0       1];
 
 defgrad = [1+ae(1) ae(2) ae(3) 1+ae(4)]';

end

