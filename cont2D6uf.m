function fe = cont2D6uf(ec,t,stress)
x1 = ec(1,1);
x2 = ec(1,2);
x3 = ec(1,3);
x4 = ec(1,4);
x5 = ec(1,5);
x6 = ec(1,6);
x = [x1,x2,x3,x4,x5,x6];

y1 = ec(2,1);
y2 = ec(2,2);
y3 = ec(2,3);
y4 = ec(2,4);
y5 = ec(2,5);
y6 = ec(2,6);
y = [y1,y2,y3,y4,y5,y6];

zmatrix = 1/3*[2    1   1;
               1    2   1;
               1    1   2];
fe = zeros(12,1);
for i = 1:3
    
    z1 = zmatrix(i,1);
    z2 = zmatrix(i,2);
    z3 = zmatrix(i,3);
    
    N1 = z1*(2*z1-1);
    N2 = z2*(2*z2-1);
    N3 = z3*(2*z3-1);
    N4 = 4*z1*z2;
    N5 = 4*z2*z3;
    N6 = 4*z3*z1;
    
    dNa_dz1 = [4*z1-1 0 0 4*z2 0 4*z3];
    dNa_dz2 = [0 4*z2-1 0 4*z1 4*z3 0];
    dNa_dz3 = [0 0 4*z3-1 0 4*z2 4*z1];
    
    
    %Components for jacobian
    dx_dz1 = dNa_dz1*x';
    dx_dz2 = dNa_dz2*x';
    dx_dz3 = dNa_dz3*x';
    
    dy_dz1 = dNa_dz1*y';
    dy_dz2 = dNa_dz2*y';
    dy_dz3 = dNa_dz3*y';
    
    
    
    
    %Jacobian and its invers
    J = [   1  ,   1   ,   1   ;
         dx_dz1, dx_dz2, dx_dz3;
         dy_dz1, dy_dz2, dy_dz3];
    
    Jinv = inv(J);
    
    P = [0,1,0;0,0,1]*Jinv';
    
    dN = P*[4*z1-1 0 0 4*z2 0 4*z3;
            0 4*z2-1 0 4*z1 4*z3 0;
            0 0 4*z3-1 0 4*z2 4*z1];
    
        
    B0 = [dN(1,1),    0   ,dN(1,2),  0     ,dN(1,3),   0    ,dN(1,4),    0   ,dN(1,5),    0   ,dN(1,6),    0   ;            % 3x8
               0    ,dN(2,1),   0    ,dN(2,2),     0  ,dN(2,3),  0     ,dN(2,4),    0   ,dN(2,5),    0   ,dN(2,6);
            dN(2,1),dN(1,1),dN(2,2),dN(1,2),dN(2,3),dN(1,3),dN(2,4),dN(1,4),dN(2,5),dN(1,5),dN(2,6),dN(1,6)];    
        
    
    fe = fe + B0'*stress{i,1}*1/2*det(J)*t*1/3;
    
end

end