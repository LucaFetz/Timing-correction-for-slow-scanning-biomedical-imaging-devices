%% clear workspace and add library
clc
clear all
close all
path = genpath('GlobalBioIm-release');
addpath(path);

%% define parameters
param.spline_order = 3; %1 or 3. any other number will result in cubic splines (order 3)
param.noise = 0; % 0 or 1
param.noise_snr = 10; %signal to noise ratio
param.lambda = 0*3e-2; %thune regulation term. 0 for no regulation.

%% Create 2D ground truth
%define GT parameters
param.GT.sigma_x = 2;
param.GT.centering_x = 3;
param.GT.sigma_y = 2;
param.GT.centering_y = 4;
param.GT.sigma_t = 2;
param.GT.centering_t = 5;
%create GT
[measurements, f0, df0] = create_2d_ground_truth(param.GT.sigma_x,param.GT.sigma_y,param.GT.sigma_t,param.GT.centering_x,param.GT.centering_y, param.GT.centering_t, param.noise, param.noise_snr);

%% plot measurements with animation in time
[Nx, Ny, Nt] = size(df0);
Nt = Nt/(Nx*Ny);
test = permute(reshape(measurements,[Nx,Ny,Nt]),[2 1 3]);
for i = 1:Nt
    imagesc(test(:,:,i));colorbar;caxis([0 1]);
    pause(0.2);
end

%% Create forward model
H = create_2d_forward_model(Nx,Ny,Nt,param.spline_order);

%% plot H
figure, imagesc(H);

%% Find C with inverse problem

%% Reconstruct good frames

%% visualize results
