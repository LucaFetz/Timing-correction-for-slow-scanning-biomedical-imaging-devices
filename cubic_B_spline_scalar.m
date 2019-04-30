function u = cubic_B_spline_scalar(x)
x = abs(x);
if(x<1)
    u = 2/3 - (1/2)*x^2*(2-x);
elseif(x<2)
    u = 1/6*(2-x)^3;
else
    u = 0;
end
end