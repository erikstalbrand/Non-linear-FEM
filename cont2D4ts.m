function [defgrad] = cont2D4ts(ec, ed)
x1 = ec(1,1);
x2 = ec(1,2);
x3 = ec(1,3);
x4 = ec(1,4);
y1 = ec(2,1);
y2 = ec(2,2);
y3 = ec(2,3);
y4 = ec(2,4);

gausscoord = 0.577350269189626;
gaussmat = [-gausscoord -gausscoord;  % xsi eta 1
             gausscoord -gausscoord;  % xsi eta 2
             gausscoord  gausscoord;  % xsi eta 3
            -gausscoord  gausscoord]; % xsi eta 4

defgrad = cell(4,1);
for n = 1:4
    
xsi = gaussmat(n, 1);
eta = gaussmat(n, 2);

% Components for Jacobian, 5.24
dx_dxsi = (1/4)*((eta-1)*x1 + (1-eta)*x2 + (eta+1)*x3 + (-1-eta)*x4);
dx_deta = (1/4)*((xsi-1)*x1 + (-xsi-1)*x2 + (xsi+1)*x3 + (1-xsi)*x4);

dy_dxsi = (1/4)*((eta-1)*y1 + (1-eta)*y2 +  (eta+1)*y3 + (-1-eta)*y4);
dy_deta = (1/4)*((xsi-1)*y1 + (-xsi-1)*y2 + (xsi+1)*y3 + (1-xsi)*y4);

%Jacobian and its inverse
J = [dx_dxsi, dx_deta;
     dy_dxsi, dy_deta];
Jinv = inv(J);

%---------------------------------------------------
% Derivatives, 5.25
dN1 = Jinv'*[(1/4)*(eta-1);
             (1/4)*(xsi-1)];

dN2 = Jinv'*[(1/4)*(1-eta);
             (1/4)*(-xsi-1)];

dN3 = Jinv'*[(1/4)*(eta+1);
             (1/4)*(xsi+1)];

dN4 = Jinv'*[(1/4)*(-eta-1);
             (1/4)*(1-xsi)];
%---------------------------------------------------
H0 = [dN1(1),  0  ,dN2(1),  0  ,dN3(1),  0  ,dN4(1),  0  ;              % 4x8
      dN1(2),  0  ,dN2(2),  0  ,dN3(2),  0  ,dN4(2),  0  ;
        0  ,dN1(1),  0  ,dN2(1),  0  ,dN3(1),  0  ,dN4(1);
        0  ,dN1(2),  0  ,dN2(2),  0  ,dN3(2),  0  ,dN4(2)];

AES = H0*ed';                                   % 4x1

defgrad{n} = [1+AES(1) AES(2) AES(3) 1+AES(4)]';
end
end

