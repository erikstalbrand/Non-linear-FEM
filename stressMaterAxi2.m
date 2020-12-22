function [stress] = stressMaterAxi2(flag, mpara, defgrad)

E=mpara(1);
v=mpara(2);

F = [defgrad(1), defgrad(2),      0    ;
     defgrad(3), defgrad(4),      0    ;
     0         ,     0     , defgrad(5)];
 
%Jacobian
J=det(F);   

%Shear modulus
k = E/(3*(1-2*v));

%Shear modulus
mu = E/(2*(1+v));

%Greens strain
C=F'*F;

% E=0.5*(C-eye(3));

I1c = C(1,1)+C(2,2)+C(3,3);

%Calculate stress tensors
S=k/2*(J^2-1)*inv(C)+mu*J^-(2/3)*(eye(3)-I1c/3*inv(C)); %Second Piola-Kirchhoff

sigma=1/J*F*S*F';    %Cauchy

tau=J*sigma;       %Kirchoff

if flag == 1
    
    stress = [S(1,1);
              S(2,2);
              S(3,3);
              S(1,2)];
end        
if flag == 3
        
    stress = [sigma(1,1);
              sigma(2,2);
              sigma(3,3);
              sigma(1,2)];  
end       
if flag == 2
    stress = [tau(1,1);
              tau(2,2);
              tau(3,3);
              tau(1,2)];  
end
end