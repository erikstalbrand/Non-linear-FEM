function [D] = dMater2D2(flag,mpara,defgrad) %Neo-Hookean Model

E=mpara(1);
v=mpara(2);

F = [defgrad(1), defgrad(2), 0;
     defgrad(3), defgrad(4), 0;
     0         ,     0     , 1];

 J=det(F);   %Jacobian

%Shear modulus
k = E/(3*(1-2*v));

%Shear modulus
mu = E/(2*(1+v));

%Greens strain
C = F'*F;
Cinv=inv(C);
trC = C(1,1)+C(2,2)+C(3,3);

%Invariants a1,a2,a3

a1 = k*J^2+2*mu/9*J^-(2/3)*trC;
a2 = 2*mu/3*J^-(2/3);
a3 = mu/3*J^-(2/3)*trC-k/2*(J^2-1);

%Kroneckers delta
d=[1,0;
   0,1];

D_func_tl=@(I,J,K,L) 0;

D_func_tl = @(I,J,K,L) D_func_tl(I,J,K,L) + a1*Cinv(I,J)*Cinv(K,L)-a2*(d(I,J)*Cinv(K,L)+Cinv(I,J)*d(K,L))+a3*(Cinv(I,K)*Cinv(J,L)+Cinv(I,L)*Cinv(J,K));

if flag == 1

  D = [D_func_tl(1,1,1,1),D_func_tl(1,1,2,2),D_func_tl(1,1,1,2);
       D_func_tl(2,2,1,1),D_func_tl(2,2,2,2),D_func_tl(2,2,1,2);
       D_func_tl(1,2,1,1),D_func_tl(1,2,2,2),D_func_tl(1,2,1,2)];
    
end
if flag==2
    F = [defgrad(1),defgrad(2),0;
         defgrad(3),defgrad(4),0;
             0     ,     0    ,1];

    D_ul_func= @(i,j,k,l) 0;

    for I=1:2
        for J=1:2
            for K=1:2
                for L=1:2
                    D_ul_func = @(i,j,k,l) D_ul_func(i,j,k,l) + (1/det(F))*F(i,I)*F(j,J)*D_func_tl(I,J,K,L)*F(k,K)*F(l,L);
                end
            end
        end
    end
    D  =  [D_ul_func(1,1,1,1),D_ul_func(1,1,2,2),D_ul_func(1,1,1,2);
           D_ul_func(2,2,1,1),D_ul_func(2,2,2,2),D_ul_func(2,2,1,2);
           D_ul_func(1,2,1,1),D_ul_func(1,2,2,2),D_ul_func(1,2,1,2)];
end
end

