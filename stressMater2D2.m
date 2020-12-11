function [stress] = stressMater2D2(flag, mpara, defgrad)   % Neo-Hookean
I = eye(3);
F = [defgrad(1) defgrad(2) 0;
     defgrad(3) defgrad(4) 0;
        0          0       1];                              % 1 is (1,1), 4 is (2,2), 2 is (1,2), 3 is (2,1)

C = F'*F;
mu = mpara(1)/(2*(1+mpara(2)));
I_1C = C(1, 1) + C(2, 2) + C(3, 3);

k = mpara(1)/(3*(1-2*mpara(2)));
J = det(F);

S = k/2*(J^2-1)*inv(C)+mu*J^(-2/3)*(I-I_1C/3*inv(C));
tau = F*S*F';

if flag == 1                                                % Second Piola Kirchoff
    stress = [S(1, 1) S(2, 2) S(1, 2)]';
    
else if flag == 2                                           % Kirchoff stress, tau
        stress = [tau(1, 1) tau(2, 2) tau(1, 2)]';
        
    else                                                    % Cauchy stress, sigma
        stress = [tau(1, 1) tau(2, 2) tau(1, 2)]'/J;
    end
end
end