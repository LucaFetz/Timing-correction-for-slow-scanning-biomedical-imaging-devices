% clear previous variables if they exist
clear all
close all
% define variable and functions
f1 = @(y,z) gaussmf(y, [2 5]).*gaussmf(z, [2 5]);
syms x t
f0(x,t) = gaussmf(x, [2 5])*gaussmf(t, [2 5]);
dx = (0:1:10);
dt = (0:1:10);
f0(dx,dt)
% plot function
fsurf(f0,[0 10])

%way of building measurements from H 
tmp = df2';
ym = tmp(:);
ym = H*tmp(:);