function u = B_spline(x)
u = 1-abs(x);
u(u<0)=0;
end