function [stress] = stressMater2D1(flag,mpara,defgrad)      % St. Venant-Kirchoff
I = eye(3);
F = [defgrad(1) defgrad(2) 0;
     defgrad(3) defgrad(4) 0;
     0              0      1];                              % 1 is (1,1), 4 is (2,2), 2 is (1,2), 3 is (2,1)

C = F'*F;
E = 1/2*(C-I);
I_1E = E(1, 1) + E(2, 2) + E(3, 3);
mu = mpara(1)/(2*(1+mpara(2)));

S = 2*mu*(E + mpara(2)/(1-2*mpara(2))*I_1E*I);              
tau = F*S*F';                                               
J = det(F);
if flag == 1                                                % Second Piola Kirchoff

    stress = [S(1, 1) S(2, 2) S(1, 2)]';
    
else if flag == 2                                           % Kirchoff stress, tau
        stress = [tau(1, 1) tau(2, 2) tau(1, 2)]';
    else                                                    % Cauchy stress, sigma
        stress = [tau(1, 1) tau(2, 2) tau(1, 2)]'/J;
    end
end
end