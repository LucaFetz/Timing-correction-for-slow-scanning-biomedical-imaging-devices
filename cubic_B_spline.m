function u = cubic_B_spline(x)
tmp1 = x;
tmp2 = x;
tmp1(abs(tmp1)<1)= 2/3 - (1/2).*tmp1(abs(tmp1)<1).^2.*(2-abs(tmp1(abs(tmp1)<1))); %tmp1(abs(tmp1)<1) is used to take only values smaller than 1. In the right side it is used to modify them in function of themselves
tmp1(abs(tmp1)>=1) = 0;%will work as result of precedent line can't be >=1.

tmp2(abs(tmp2)<1) = 0;
tmp2(abs(tmp2)>=2) = 0;
tmp2(tmp2 ~=0) = 1/6.*(2-abs(tmp2(tmp2 ~=0))).^3;
u = tmp1 + tmp2;
end