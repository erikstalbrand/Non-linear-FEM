function [D] = dMater2D1(flag,mpara,defgrad)

D_1111 = 1 - mpara(2);
D_1122 = mpara(2);
D_1112 = 0;
D_2211 = mpara(2);
D_2222 = 1 - mpara(2);
D_2212 = 0;
D_1211 = 0;
D_1222 = 0;
D_1212 = 1/2*(1-2*mpara(2));

D = mpara(1)/((1+mpara(2))*(1-2*mpara(2)))*[D_1111 D_1122 D_1112;
     D_2211 D_2222 D_2212;
     D_1211 D_1222 D_1212];
end

