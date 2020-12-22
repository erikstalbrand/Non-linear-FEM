function Ke = contAxi3te(ec, t, D, ed, stress)
r1=ec(1,1);
r2=ec(1,2);
r3=ec(1,3);
z1=ec(2,1);
z2=ec(2,2);
z3=ec(2,3);

r0 = (r1+r2+r3)/3;
z0 = (z1+z2+z3)/3;

A = [ 1, 1, 1;
     r1,r2,r3;
     z1,z2,z3];
    
    
Ae = 0.5*det(A);

R0=[stress(1),stress(4),  0  ,  0  ,    0     ;
    stress(4),stress(2),  0  ,  0  ,    0     ;
      0  ,  0  ,stress(1),stress(4),    0     ;
      0  ,  0  ,stress(4),stress(2),    0     ;
      0  ,  0  ,    0    ,    0    ,stress(3)];

N1 = 1/(2*Ae)*(r2*z3-r3*z2+(z2-z3)*r0+(r3-r2)*z0);

N2 = 1/(2*Ae)*(r3*z1-r1*z3+(z3-z1)*r0+(r1-r3)*z0);

N3 = 1/(2*Ae)*(r1*z2-r2*z1+(z1-z2)*r0+(r2-r1)*z0);

B0le=1/(2*Ae)*[z2-z3     ,  0  ,z3-z1     ,  0  ,z1-z2     ,  0  ;
                 0       ,r3-r2,  0       ,r1-r3,  0       ,r2-r1;
               2*Ae*N1/r0,  0  ,2*Ae*N2/r0,  0  ,2*Ae*N3/r0,  0  ;
               r3-r2     ,z2-z3,r1-r3     ,z3-z1,r2-r1     ,z1-z2];

H0=1/(2*Ae)*[z2-z3       ,  0  ,z3-z1     ,  0  ,z1-z2     ,  0  ; 
             r3-r2       ,  0  ,r1-r3     ,  0  ,r2-r1     ,  0  ;
               0         ,z2-z3,  0       ,z3-z1,  0       ,z1-z2;
               0         ,r3-r2,  0       ,r1-r3,  0       ,r2-r1;
               2*Ae*N1/r0,  0  ,2*Ae*N2/r0,  0  ,2*Ae*N3/r0,  0  ];

AES=H0*ed';


%HUR PLACERAR MAN UT AES HÄR??
AE=[AES(1),   0  ,AES(3),   0  ,  0   ;
      0   ,AES(2),  0   ,AES(4),  0   ;
      0   ,   0  ,  0   ,   0  ,AES(5);
    AES(2),AES(1),AES(4),AES(3),  0   ];
    

B0=B0le+AE*H0;
          
Ke = 2*pi*r0*(B0'*D*B0+H0'*R0*H0)*Ae;
end
