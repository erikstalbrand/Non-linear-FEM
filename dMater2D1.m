function [D] = dMater2D1(flag, mpara, defgrad)

mu = mpara(1)/(2*(1+mpara(2)));
d = eye(2);

D_func_tl = @(I,J,K,L) 0;
D_func_tl = @(I,J,K,L) D_func_tl(I,J,K,L) + 2*mu*(1/2*(d(I,K)*d(J,L)+d(I,L)*d(J,K))+mpara(2)/(1-2*mpara(2))*d(I,J)*d(K,L));

if flag == 1
D =[D_func_tl(1,1,1,1), D_func_tl(1,1,2,2), D_func_tl(1,1,1,2);
    D_func_tl(2,2,1,1), D_func_tl(2,2,2,2), D_func_tl(2,2,1,2);
    D_func_tl(1,2,1,1), D_func_tl(1,2,2,2), D_func_tl(1,2,1,2)];
end

if flag == 2
F = [defgrad(1) defgrad(2) 0;
     defgrad(3) defgrad(4)  0;
     0          0           1];

Jac = det(F);
D_func_ul = @(i, j, k, l) 0;

    for I = 1:2
        for J = 1:2
            for K = 1:2
                for L = 1:2                
                  D_func_ul = @(i,j,k,l) D_func_ul(i,j,k,l) + (1/Jac)*(F(i,I)*F(j,J)*D_func_tl(I, J, K , L)*F(k,K)*F(l,L));
                end
            end
        end
    end

D = [D_func_ul(1,1,1,1), D_func_ul(1,1,2,2), D_func_ul(1,1,1,2);
     D_func_ul(2,2,1,1), D_func_ul(2,2,2,2), D_func_ul(2,2,1,2);
     D_func_ul(1,2,1,1), D_func_ul(1,2,2,2), D_func_ul(1,2,1,2)];
end
end

